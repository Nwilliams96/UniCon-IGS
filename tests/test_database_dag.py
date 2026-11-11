"""Isolated dry-runs of the real workflow; no reference downloads or QIIME runs.

The workflow copy receives synthetic classifiers and a minimal generated-config
fixture. Tests are skipped when Snakemake is not installed in this interpreter.
"""

import importlib.util
import itertools
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
class DatabaseDagTests(unittest.TestCase):
    def setup_project(self, root):
        root = Path(root)
        shutil.copytree(ROOT / "workflow", root / "workflow")
        config_dir = root / "config"
        (config_dir / "schemas").mkdir(parents=True)
        # Configuration-builder tests exercise its full public schema. Here
        # only valid sample identifiers are needed to parse the actual DAG.
        schema_version = '$schema: "https://json-schema.org/draft/2020-12/schema"\n'
        (config_dir / "schemas/config.schema.yml").write_text(schema_version + "type: object\n")
        (config_dir / "schemas/samples.schema.yml").write_text(
            schema_version + "type: object\nproperties:\n  sample:\n    type: string\nrequired: [sample]\n"
        )
        (config_dir / "samples.tsv").write_text("sample\n001\n")
        (config_dir / "prok_and_euk_SSU_amplicon_molarities.tsv").write_text(
            "sample_type\tamount_pM\n16S\t1000\n18S\t1000\n"
        )
        shared = root / "shared databases"
        (shared / "bbsplit-db/EUK-PROK-bbsplit-db").mkdir(parents=True)
        return dict(
            samplesheet="config/samples.tsv",
            projectName="DAG-fixture", database_dir=str(shared),
            rawdatadir=str(root / "raw"), R1file_ending=".R1.fastq.gz",
            R2file_ending=".R2.fastq.gz", use_internal_standards=False,
            use_preexisting_bbsplit_database=True,
            trunclens=dict(truncR1=220, truncR2=180),
            fwdPrimer="GTGYCAGCMGCCGCGGTAA", revPrimer="CCGYCAATTYMTTTRAGTTT",
            qiime2version="../envs/rachis-qiime2-linux-64-2026.7.yml",
            SILVAversion="138.2", PR2version="5.1.1",
        )

    def set_mode(self, root, config, family, mode):
        lower = family.lower()
        config["use_preexisting_" + lower + "_database"] = mode == "existing"
        if mode == "existing":
            artifact = Path(root) / "provided" / (family + ".qza")
            artifact.parent.mkdir(exist_ok=True)
            with zipfile.ZipFile(artifact, "w") as archive:
                archive.writestr("fixture/metadata.yaml", "type: TaxonomicClassifier\n")
                archive.writestr("fixture/data/classifier.pkl", b"fixture")
                archive.writestr(
                    "fixture/provenance/artifacts/ancestor/metadata.yaml",
                    "type: FeatureData[Sequence]\n",
                )
            config[lower + "_classifier"] = str(artifact)
        elif mode == "download":
            config[family + "_download_link"] = f"https://references.example/{family}.qza"

    def dry_run(self, root, config, targets=("ensure_silva_classifier", "ensure_pr2_classifier")):
        (Path(root) / "config/config.yml").write_text(json.dumps(config))
        process = subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", "workflow/Snakefile",
             "--dry-run", "--forceall", "--nolock", "--shared-fs-usage", "input-output",
             "persistence", "--cores", "8", "--printshellcmds", *targets],
            cwd=root, text=True, capture_output=True, timeout=60,
        )
        output = process.stdout + process.stderr
        self.assertEqual(process.returncode, 0, output)
        return output

    def test_all_nine_choices_schedule_only_the_selected_preparation_route(self):
        with tempfile.TemporaryDirectory() as tmp:
            for silva_mode, pr2_mode in itertools.product(("existing", "download", "train"), repeat=2):
                with self.subTest(SILVA=silva_mode, PR2=pr2_mode):
                    root = Path(tmp) / (silva_mode + "-" + pr2_mode)
                    config = self.setup_project(root)
                    self.set_mode(root, config, "SILVA", silva_mode)
                    self.set_mode(root, config, "PR2", pr2_mode)
                    output = self.dry_run(root, config)
                    for family, mode in (("SILVA", silva_mode), ("PR2", pr2_mode)):
                        download = "rule download_" + family + "_classifier:"
                        train = "rule train_classifier:" if family == "SILVA" else "rule train_classifier_pr2:"
                        extract = "rule extract_primers:" if family == "SILVA" else "rule extract_primers_pr2:"
                        self.assertEqual(download in output, mode == "download", output)
                        self.assertEqual(train in output, mode == "train", output)
                        self.assertEqual(extract in output, mode == "train", output)

    def test_named_qiime_environment_parses_the_real_workflow_without_a_yaml_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            for family in ("SILVA", "PR2"):
                self.set_mode(tmp, config, family, "existing")
            prefix = Path(tmp) / "envs/qiime-fixture"
            (prefix / "conda-meta").mkdir(parents=True)
            conda = Path(tmp) / "fixture-conda"
            conda.write_text(
                "#!" + sys.executable + "\nimport json, sys\n"
                + "print(json.dumps(" + repr({"envs": [str(prefix)]}) + ") if sys.argv[1] == 'info' else "
                + "json.dumps([{'name': 'qiime2', 'version': '2026.7'}, "
                + "{'name': 'scikit-learn', 'version': '1.7.1'}]))\n"
            )
            conda.chmod(0o700)
            config["qiime2version"] = "qiime-fixture"
            with mock.patch.dict(os.environ, {"CONDA_EXE": str(conda)}):
                output = self.dry_run(tmp, config)
            self.assertIn("rule ensure_silva_classifier:", output)
            self.assertIn("rule ensure_pr2_classifier:", output)

    def test_explicit_legacy_yaml_sets_real_workflow_classifier_cache_identity(self):
        for filename in ("qiime2-amplicon-ubuntu-2024-5-conda.yml",
                         "qiime2-amplicon-ubuntu-2025-7-conda.yml"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                config = self.setup_project(tmp)
                definition = Path(tmp) / filename
                # Deliberately distinct fixture content: no environment is
                # installed or solved during this workflow dry-run.
                definition.write_text("name: explicit-legacy-fixture\ndependencies: []\n")
                config["qiime2version"] = str(definition)
                snakefile = Path(tmp) / "workflow/Snakefile"
                snakefile.write_text(snakefile.read_text() + '''
rule inspect_selected_environment_fixture:
    input:
        rules.ensure_silva_classifier.output,
        rules.ensure_pr2_classifier.output
    output:
        "selected-environment-fixture.txt"
    params:
        requested=config["qiime2_environment"]["requested"],
        classifier=PR2_CLASSIFIER
    shell:
        "echo {params.requested:q} {params.classifier:q} > {output:q}"
''')
                output = self.dry_run(tmp, config, ("inspect_selected_environment_fixture",))
                identity = hashlib.sha256(definition.read_bytes()).hexdigest()[:12]
                self.assertIn(str(definition), output)
                self.assertIn("-" + identity + ".qza", output)
                self.assertIn("rule train_classifier:", output)
                self.assertIn("rule train_classifier_pr2:", output)

    def test_current_internal_standard_lists_connect_the_full_results_dag(self):
        for ids in (["BP"], ["BP", "DR"], ["Custom-A", "Custom_B", "SpikeIn3", "SpikeIn4"]):
            with self.subTest(ids=ids), tempfile.TemporaryDirectory() as tmp:
                config = self.setup_project(tmp)
                config.update(use_internal_standards=True, intstds=ids)
                for family in ("SILVA", "PR2"):
                    self.set_mode(tmp, config, family, "existing")
                (Path(tmp) / "config/samples.tsv").write_text(
                    "sample\t" + "\t".join(id + "_ng" for id in ids) + "\n"
                    + "001\t" + "\t".join("1" for id in ids) + "\n"
                )
                (Path(tmp) / "config/internal_stds.tsv").write_text(
                    "internal_std_ID\trRNA_copy_number\tgenome_len_bp\tfull_SSU_sequence\n"
                    + "".join(f"{id}\t1\t100\tACGT\n" for id in ids)
                )
                raw = Path(config["rawdatadir"])
                raw.mkdir()
                for direction in ("R1", "R2"):
                    (raw / ("001." + direction + ".fastq.gz")).touch()
                output = self.dry_run(tmp, config, ("all",))
                for rule in ("prepare_internal_standard_fastas", "intstd_correct_data",
                             "generate_pipeline_report", "export_results_for_download"):
                    self.assertIn("rule " + rule + ":", output)
                self.assertIn("DAG-fixture." + "_".join(ids) + ".ISD_corrected_long_data.tsv", output)
                for id in ids:
                    self.assertIn("config/intstd_fastas/" + id + ".fasta", output)

    def test_full_length_training_bypasses_primer_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            config.update(silva_reference_mode="full_length", pr2_reference_mode="full_length")
            output = self.dry_run(tmp, config)
            self.assertIn("rule train_classifier:", output)
            self.assertIn("rule train_classifier_pr2:", output)
            self.assertNotIn("rule extract_primers:", output)
            self.assertNotIn("rule extract_primers_pr2:", output)
            self.assertNotIn("rule dereplicated_sliced_data:", output)
            self.assertNotIn("rule dereplicate_extracted_pr2_reads:", output)

    def test_native_version_capability_check_precedes_raw_read_processing(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            config["SILVAversion"] = "144"
            self.set_mode(tmp, config, "PR2", "existing")
            raw = Path(config["rawdatadir"])
            raw.mkdir()
            for direction in ("R1", "R2"):
                (raw / ("001." + direction + ".fastq.gz")).touch()
            output = self.dry_run(tmp, config, ("results/00-trimmed/001.qc.txt",))
            self.assertIn("rule check_silva_reference_capability:", output)
            self.assertIn("rule cutadapt:", output)
            cutadapt = output.split("rule cutadapt:", 1)[1]
            self.assertRegex(cutadapt.split("output:", 1)[0],
                             r"database-build/SILVA/[a-f0-9]+/capability.json")

    def test_completed_compatible_downloads_and_training_are_reused(self):
        spec = importlib.util.spec_from_file_location(
            "database_config", ROOT / "workflow/scripts/database_config.py"
        )
        resolver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(resolver)
        for mode in ("download", "train"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                config = self.setup_project(tmp)
                for family in ("SILVA", "PR2"):
                    self.set_mode(tmp, config, family, mode)
                env = Path(tmp) / "workflow/envs/rachis-qiime2-linux-64-2026.7.yml"
                settings = dict(config, _database_environment_id=hashlib.sha256(env.read_bytes()).hexdigest()[:12])
                for family in ("SILVA", "PR2"):
                    path = Path(resolver.resolve_database(settings, family)["classifier"])
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with zipfile.ZipFile(path, "w") as archive:
                        archive.writestr("fixture/metadata.yaml", "type: TaxonomicClassifier\n")
                        archive.writestr("fixture/data/classifier.pkl", b"completed fixture")
                output = self.dry_run(tmp, config)
                self.assertIn("rule ensure_silva_classifier:", output)
                self.assertIn("rule ensure_pr2_classifier:", output)
                for forbidden in ("download_SILVA_classifier", "download_PR2_classifier",
                                  "train_classifier", "train_classifier_pr2",
                                  "check_silva_reference_capability", "check_pr2_reference_capability"):
                    self.assertNotIn("rule " + forbidden + ":", output)

    def test_full_results_dag_includes_normalized_pr2_inputs_and_matching_species_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.setup_project(tmp)
            for family in ("SILVA", "PR2"):
                self.set_mode(tmp, config, family, "existing")
            raw = Path(config["rawdatadir"])
            raw.mkdir()
            for direction in ("R1", "R2"):
                (raw / ("001." + direction + ".fastq.gz")).touch()
            output = self.dry_run(tmp, config, ("all",))
            for rule in ("normalize_canonical_16S_PR2_rank_prefixes",
                         "normalize_unmerged_16S_PR2_rank_prefixes",
                         "normalize_18S_PR2_rank_prefixes",
                         "import_normalized_18S_PR2_taxonomy",
                         "resolve_PR2_plastids", "split_resolved_prok_categories",
                         "generate_pipeline_report", "export_results_for_download"):
                self.assertIn("rule " + rule + ":", output)
            self.assertNotIn("rule download_SILVA_144_species_reference:", output)
            self.assertNotIn("rule assign_SILVA_144_species:", output)
            resolver = output.split("rule resolve_PR2_plastids:", 1)[1].split("output:", 1)[0]
            self.assertIn(".snakemake/classification-exports/canonical-16S/PR2/taxonomy.tsv", resolver)
            supplement = output.split("rule export_unmerged_16S_supplement:", 1)[1].split("output:", 1)[0]
            self.assertIn(".snakemake/classification-exports/unmerged-16S/PR2/taxonomy.tsv", supplement)
            self.assertIn(
                ".snakemake/classification-exports/unmerged-16S/sequences/linked-dna-sequences.fasta",
                supplement,
            )
            self.assertNotIn("/sequences/dna-sequences.fasta", supplement)
            report = output.split("rule generate_pipeline_report:", 1)[1].split("output:", 1)[0]
            self.assertIn("DAG-fixture.unmerged_linked_16S_summary.tsv", report)
            export = output.split("rule export_results_for_download:", 1)[1].split("log:", 1)[0]
            for filename in ("DAG-fixture.unmerged_linked_16S.tsv",
                             "DAG-fixture.unmerged_linked_16S_summary.tsv"):
                self.assertIn("results/02-proks/11-unmerged-16S-rescue/" + filename, export)
                self.assertIn("DAG-fixture-Results-Export/" + filename, export)
            euk_export = output.split("rule export_tax_convert_biom_euk:", 1)[1].split("output:", 1)[0]
            self.assertIn("DAG-fixture_PR2.normalized.classified.qza", euk_export)
            self.assertIn("rule assign_SILVA_species:", output)

    def test_report_only_rerun_uses_custom_samples_and_cannot_schedule_sequence_processing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.setup_project(root)
            for database in ("SILVA", "PR2"):
                self.set_mode(root, config, database, "existing")
            config["samplesheet"] = "config/custom-samples.tsv"
            (root / "config/samples.tsv").rename(root / config["samplesheet"])
            (root / "config/config.yml").write_text(json.dumps(config))
            project = config["projectName"]
            files = [
                f"results/{project}.eukfrac-per-sample.tsv",
                f"results/02-proks/04-DADA2d-plaintext-exports/{project}.16S.latest_stats.tsv",
                f"results/02-euks/09-DADA2d-plaintext-exports/{project}.18S.latest_stats.tsv",
                f"results/03-merged/{project}.16S_18S_correction_factors.tsv",
                f"results/02-proks/10-exports/{project}.SILVA-138.2.species-assignment-summary.tsv",
                f"results/02-proks/09-subsetting/tax-merged/{project}.PR2-plastid-routing-audit.tsv",
                f"results/02-proks/09-subsetting/tax-merged/{project}.PR2-plastid-routing-summary.tsv",
                f"results/00-trimmed/001.qc.txt",
                f"results/04-formatted/{project}.long_data.tsv",
                f"results/04-formatted/{project}.asv_sequences.tsv",
                f"results/02-proks/11-unmerged-16S-rescue/{project}.unmerged_linked_16S.tsv",
                f"results/02-proks/11-unmerged-16S-rescue/{project}.unmerged_linked_16S_summary.tsv",
            ]
            for relative in files:
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.touch()
            for relative in ("results/02-proks/02-quality-plots-R1-R2",
                             "results/02-euks/02-quality-plots-R1-R2",
                             "results/02-euks/07-quality-plots-concat"):
                (root / relative).mkdir(parents=True)
            process = subprocess.run(
                [sys.executable, "-m", "snakemake", "export_results_for_download",
                 "--snakefile", "workflow/Snakefile", "--cores", "1", "--dry-run",
                 "--shared-fs-usage", "input-output", "persistence",
                 "--forcerun", "generate_pipeline_report", "export_results_for_download",
                 "--allowed-rules", "generate_pipeline_report", "export_results_for_download"],
                cwd=root, text=True, capture_output=True, timeout=60,
            )
            output = process.stdout + process.stderr
            self.assertEqual(process.returncode, 0, output)
            self.assertIn("config/custom-samples.tsv", output)
            self.assertNotIn("config/samples.tsv", output)
            scheduled = re.findall(r"(?:local)?rule (\w+):", output)
            self.assertEqual(sorted(scheduled), ["export_results_for_download", "generate_pipeline_report"])
            self.assertFalse(Path(config["rawdatadir"]).exists())


if __name__ == "__main__":
    unittest.main()
