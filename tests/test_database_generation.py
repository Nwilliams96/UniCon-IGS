from contextlib import chdir
import importlib.util
import hashlib
import os
import re
import tempfile
import unittest
import zipfile
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFORMAT_SCRIPT = ROOT / "workflow/scripts/reformat_pr2_reference.py"
SPEC = importlib.util.spec_from_file_location("reformat_pr2_reference", REFORMAT_SCRIPT)
REFORMAT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REFORMAT)

DOWNLOAD_SCRIPT = ROOT / "workflow/scripts/download_verified_reference.py"
DOWNLOAD_SPEC = importlib.util.spec_from_file_location(
    "download_verified_reference", DOWNLOAD_SCRIPT
)
DOWNLOAD = importlib.util.module_from_spec(DOWNLOAD_SPEC)
DOWNLOAD_SPEC.loader.exec_module(DOWNLOAD)

SAVE_SCRIPT = ROOT / "workflow/scripts/save_classifier.py"
SAVE_SPEC = importlib.util.spec_from_file_location(
    "save_classifier", SAVE_SCRIPT
)
SAVE = importlib.util.module_from_spec(SAVE_SPEC)
SAVE_SPEC.loader.exec_module(SAVE)


class DatabaseGenerationTests(unittest.TestCase):
    @staticmethod
    def write_classifier(path, semantic_type="TaxonomicClassifier"):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("fixture/metadata.yaml", "type: " + semantic_type + "\n")
            archive.writestr("fixture/data/classifier.pkl", b"test classifier")
            archive.writestr(
                "fixture/provenance/artifacts/ancestor/metadata.yaml",
                "type: FeatureData[Sequence]\n",
            )
        return path

    def test_database_rules_are_self_contained_and_use_supported_parallelism(self):
        classification = (ROOT / "workflow/rules/prepare_classification_dbs.smk").read_text(
            encoding="utf-8"
        )
        common = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        downloads = (ROOT / "workflow/rules/download_databases.smk").read_text(
            encoding="utf-8"
        )
        bbsplit = (ROOT / "workflow/rules/prepare_bbsplit_db.smk").read_text(
            encoding="utf-8"
        )
        classifier_saver = (ROOT / "workflow/scripts/save_classifier.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("tax-classifier-construction", classification)
        for rule_name in (
            "download_pr2_reference",
            "normalize_pr2_reference_taxonomy",
            "cull_pr2_seqs",
            "derep_seqs_taxonomy",
            "extract_primers_pr2",
            "dereplicate_extracted_pr2_reads",
            "train_classifier_pr2",
            "ensure_pr2_classifier",
        ):
            self.assertIn(f"rule {rule_name}:", classification)
        for rule_name in (
            "download_SILVA", "reverse_transcribe", "qc_seqs_cull",
            "qc_seqs_filter", "qc_seqs_dereplicate", "extract_primers",
            "dereplicated_sliced_data", "train_classifier", "ensure_silva_classifier",
            "download_SILVA_classifier", "download_PR2_classifier",
        ):
            self.assertIn(f"rule {rule_name}:", classification)
        for removed_rule in (
            "decompress_pr2_reference", "clean_pr2_fasta_extract_headers",
            "import_pr2_fasta", "import_pr2_taxonomy",
        ):
            self.assertNotIn(f"rule {removed_rule}:", classification)
        self.assertGreaterEqual(classification.count("threads: 8"), 1)
        self.assertEqual(classification.count("fit-classifier-naive-bayes"), 2)
        self.assertIn("--o-read-extraction-stats", classification)
        self.assertIn("--p-rank-handles disable", classification)
        self.assertIn("get-pr2-data", classification)
        self.assertIn("get-silva-data", classification)
        self.assertIn("--p-labels Bacteria Archaea Eukaryota --p-min-lens 1200 900 1400", classification)
        self.assertIn("fcntl.flock", classifier_saver)
        self.assertIn("os.replace(temporary, destination)", classifier_saver)
        self.assertIn('["qiime", "tools", "validate", str(path)]', classifier_saver)
        self.assertFalse((ROOT / "workflow/scripts/ensure_pr2_classifier.py").exists())
        self.assertIn("rule download_SILVA_144_species_reference:", classification)
        self.assertIn('config.get("SILVAversion", "144")', common)
        self.assertIn("silva_v144_assignSpecies.fa.gz", common)
        self.assertIn("444de7c0cce0b66addda7a3f8b38e012", common)
        self.assertIn("rachis-qiime2-linux-64-2026.7.yml", common)
        self.assertIn("marker=PR2_READY_MARKER", classification)
        self.assertIn("marker=SILVA_READY_MARKER", classification)
        self.assertIn("threads: 8", bbsplit)
        self.assertIn("bbsplit.sh build=1 threads={threads}", bbsplit)
        self.assertIn("rule initialize_database_directories:", downloads)
        self.assertIn('DATABASE_PREFIX + "classification/SILVA"', downloads)
        self.assertNotIn("rule download_pr2:", downloads)
        self.assertNotIn('"v7.2.0/utils/aria2c"', downloads)
        self.assertIn('"logs/bbsplit_prok_db_download.log"', downloads)
        self.assertIn('"logs/bbsplit_euk_db_download.log"', downloads)
        self.assertIn("46683398d8859aeef635fff49112c4d6453d4ebd70bf8b9fee006cb5b84273b1", downloads)
        self.assertIn("c4a396035271db516a0054f1d750d3081ca5a11831b44756e91590a79ed0f109", downloads)
        self.assertRegex(classification, r'script:\s+"\.\./scripts/download_verified_reference\.py"')
        self.assertIn("USE_PREEXISTING_BBSPLIT_DATABASE", common)
        self.assertNotIn("USE_PREEXISTING_SILVA_DATABASE", common)
        self.assertNotIn("USE_PREEXISTING_PR2_DATABASE", common)
        snakefile = (ROOT / "workflow/Snakefile").read_text(encoding="utf-8")
        self.assertIn("if not USE_PREEXISTING_BBSPLIT_DATABASE:", snakefile)

        species_script = (ROOT / "workflow/scripts/assign_silva_species.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("allowMultiple = FALSE", species_script)
        self.assertIn("tryRC = TRUE", species_script)
        self.assertIn("Genus_Agrees", species_script)
        self.assertIn("assign_species_in_chunks", species_script)
        self.assertIn("stage_reference_locally", species_script)
        self.assertIn("Using node-local SILVA species reference", species_script)
        self.assertIn("gc(full = TRUE)", species_script)
        species_rule = (
            ROOT / "workflow/rules/02-denoise-and-export-prok.smk"
        ).read_text(encoding="utf-8")
        species_rule = species_rule.split("rule assign_SILVA_species:", 1)[1].split(
            "\nrule ", 1
        )[0]
        self.assertIn("chunk_size=500", species_rule)
        self.assertIn("mem_mb=110000", species_rule)
        qiime_env = (
            ROOT / "workflow/envs/rachis-qiime2-linux-64-2026.7.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("bioconductor-dada2=", qiime_env)
        self.assertIn("q2-dada2=2026.7.0", qiime_env)
        self.assertIn("q2-feature-classifier=2026.7.0", qiime_env)
        self.assertIn("rescript=2026.7.0", qiime_env)
        self.assertIn("scikit-learn=1.7.1", qiime_env)
        self.assertNotIn("\n- deblur=", qiime_env)
        self.assertNotIn("\n- q2-deblur=", qiime_env)
        self.assertNotIn("\n- sortmerna=", qiime_env)

        self.assertNotIn(
            "input:\n        temp(",
            downloads,
        )

    def test_shared_bbsplit_reference_outputs_are_persistent(self):
        downloads = (ROOT / "workflow/rules/download_databases.smk").read_text()
        for rule, filename in (
            ("download_prok_db", "SILVA_132_PROK.cdhit95pc.fasta"),
            ("download_euk_db", "SILVA_132_and_PR2_EUK.cdhit95pc.fasta"),
        ):
            with self.subTest(rule=rule):
                body = downloads.split(f"rule {rule}:", 1)[1].split("\nrule ", 1)[0]
                output = body.split("output:", 1)[1].split("log:", 1)[0]
                self.assertNotRegex(output, r"\btemp\s*\(")
                self.assertIn('DATABASE_PREFIX + "bbsplit-db/' + filename + '"', output)

    def test_reference_download_is_atomic_and_checksum_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.qza"
            destination = root / "classification" / "classifier.qza"
            source.write_bytes(b"signed SILVA classifier fixture")
            checksum = hashlib.md5(source.read_bytes()).hexdigest()

            DOWNLOAD.download_verified(source.as_uri(), destination, checksum)
            self.assertEqual(destination.read_bytes(), source.read_bytes())
            self.assertFalse(destination.with_name(destination.name + ".part").exists())

            destination.unlink()
            with self.assertRaises(ValueError):
                DOWNLOAD.download_verified(source.as_uri(), destination, "0" * 32)
            self.assertFalse(destination.exists())
            self.assertFalse(destination.with_name(destination.name + ".part").exists())

            sha_destination = root / "absolute-style" / "reference.fasta.gz"
            sha_checksum = hashlib.sha256(source.read_bytes()).hexdigest()
            DOWNLOAD.download_verified(
                source.as_uri(), sha_destination, sha_checksum, "sha256"
            )
            self.assertEqual(sha_destination.read_bytes(), source.read_bytes())

    def test_shared_reference_download_uses_lock_and_unique_staging(self):
        source = DOWNLOAD_SCRIPT.read_text(encoding="utf-8")
        self.assertIn("fcntl.flock", source)
        self.assertIn("destination_lock(destination)", source)
        self.assertIn("tempfile.NamedTemporaryFile", source)
        self.assertIn("os.replace(partial, destination)", source)
        self.assertNotIn('destination.with_name(destination.name + ".part")', source)

    def test_classifier_download_without_checksum_requires_valid_type_and_qiime(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.write_classifier(root / "source.qza")
            destination = root / "absolute destination" / "classifier.qza"
            log = root / "logs" / "download.log"
            with mock.patch.object(DOWNLOAD.subprocess, "run") as validate:
                DOWNLOAD.download_verified(
                    source.as_uri(), destination, validate_classifier=True, log_path=log
                )
            self.assertEqual(destination.read_bytes(), source.read_bytes())
            self.assertEqual(validate.call_count, 1)
            self.assertEqual(validate.call_args.args[0][:3], ["qiime", "tools", "validate"])
            self.assertIn("No publisher checksum supplied", log.read_text())
            self.assertFalse(list(destination.parent.glob("*.part")))
            with self.assertRaisesRegex(ValueError, "require a published checksum"):
                DOWNLOAD.download_verified(source.as_uri(), root / "raw.fasta")

    def test_classifier_download_rejects_nonclassifier_archive_and_preserves_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.write_classifier(root / "sequences.qza", "FeatureData[Sequence]")
            destination = root / "classifier.qza"
            destination.write_bytes(b"prior reference")
            with mock.patch.object(DOWNLOAD.subprocess, "run") as qiime:
                with self.assertRaisesRegex(ValueError, "TaxonomicClassifier"):
                    DOWNLOAD.download_verified(
                        source.as_uri(), destination, validate_classifier=True
                    )
                qiime.assert_not_called()
            self.assertEqual(destination.read_bytes(), b"prior reference")
            self.assertFalse(list(root.glob("*.part")))

    def test_valid_cached_classifier_is_revalidated_without_network_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = self.write_classifier(Path(tmp) / "classifier.qza")
            with mock.patch.object(DOWNLOAD.subprocess, "run") as qiime:
                with mock.patch.object(DOWNLOAD.urllib.request, "urlopen") as network:
                    DOWNLOAD.download_verified(
                        "https://references.example/classifier.qza", destination,
                        validate_classifier=True,
                    )
                    network.assert_not_called()
                    qiime.assert_called_once()

    def test_pr2_saving_is_locked_atomic_and_reusable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "private" / "classifier.qza"
            destination = root / "shared" / "classifier.qza"
            marker = root / "ready" / "classifier.ready"
            lock = root / "shared" / ".build.lock"
            log = root / "logs" / "save.log"
            source.parent.mkdir(parents=True)
            source.write_bytes(b"validated fixture")

            def valid_if_nonempty(path, _log_handle):
                path = Path(path)
                return path.is_file() and path.stat().st_size > 0

            with mock.patch.object(SAVE, "classifier_is_valid", valid_if_nonempty):
                SAVE.save_classifier(source, destination, lock, marker, log)
                self.assertEqual(destination.read_bytes(), b"validated fixture")
                self.assertTrue(marker.exists())

                source.write_bytes(b"new staged bytes should not replace valid shared")
                SAVE.save_classifier(source, destination, lock, marker, log)
                self.assertEqual(destination.read_bytes(), b"validated fixture")
                self.assertFalse(list(destination.parent.glob("*.saving.*")))

    def test_existing_classifier_validation_never_writes_beside_the_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.write_classifier(root / "read-only source" / "provided.qza")
            before = source.read_bytes()
            marker = root / "project" / "classifier.ready"
            with mock.patch.object(SAVE.subprocess, "run") as qiime:
                qiime.return_value.returncode = 0
                SAVE.save_classifier(
                    source, source, root / "project" / "save.lock",
                    marker, root / "project" / "save.log",
                    read_only=True, label="SILVA",
                )
            self.assertTrue(marker.exists())
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(list(source.parent.iterdir()), [source])
            self.assertIn(str(source.resolve()), marker.read_text())

    def test_saving_rejects_wrong_artifact_type_before_running_qiime(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.write_classifier(root / "sequences.qza", "FeatureData[Sequence]")
            with mock.patch.object(SAVE.subprocess, "run") as qiime:
                with self.assertRaisesRegex(RuntimeError, "failed QIIME validation"):
                    SAVE.save_classifier(
                        source, root / "classifier.qza", root / "save.lock",
                        root / "ready", root / "save.log", label="SILVA",
                    )
                qiime.assert_not_called()
            self.assertFalse((root / "classifier.qza").exists())
            self.assertFalse((root / "ready").exists())

    def test_expected_download_checksum_prevents_reusing_different_valid_classifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self.write_classifier(root / "private" / "classifier.qza")
            destination = self.write_classifier(root / "shared" / "classifier.qza")
            with zipfile.ZipFile(destination, "a") as archive:
                archive.writestr("fixture/data/old.txt", "old valid classifier")
            expected = hashlib.md5(source.read_bytes()).hexdigest()
            with mock.patch.object(SAVE.subprocess, "run") as qiime:
                qiime.return_value.returncode = 0
                SAVE.save_classifier(
                    source, destination, root / "shared" / "save.lock",
                    root / "ready", root / "save.log", expected_md5=expected,
                )
            self.assertEqual(destination.read_bytes(), source.read_bytes())

    def test_long_data_parser_is_rank_named_and_silva_144_safe(self):
        source = (ROOT / "workflow/scripts/long-data-preparation.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('Kingdom = extract_prefixed_rank(Taxonomy, "k__")', source)
        self.assertIn("parse_prefixed_taxonomy", source)
        self.assertNotIn('separate(Taxonomy, c("Domain","Phylum"', source)
        self.assertIn('"Domain", "Kingdom", "Supergroup"', source)
        self.assertIn(
            '"Taxonomy", "Domain", "Kingdom", "Supergroup"', source
        )

        subset_script = (
            ROOT / "workflow/scripts/P09b-split-resolved-prok-categories.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("p__Cyanobacteria,p__Cyanobacteriota", subset_script)
        self.assertIn(":plas,o__Chloroplast", subset_script)

        prok_rules = (
            ROOT / "workflow/rules/02-denoise-and-export-prok.smk"
        ).read_text(encoding="utf-8")
        self.assertIn("rule classify_all_16S_with_PR2:", prok_rules)
        self.assertIn("rule resolve_PR2_plastids:", prok_rules)
        self.assertIn("rule split_resolved_prok_categories:", prok_rules)
        self.assertIn("all_16S_ASVs_PR2.classified.qza", prok_rules)
        self.assertIn(
            "resolved_taxonomy=rules.import_resolved_PR2_plastid_taxonomy.output.taxonomy",
            prok_rules,
        )
        self.assertNotIn("rule splitchloroplasts:", prok_rules)

        for rule_name in (
            "export_SILVA_taxonomy_for_plastid_resolution",
            "export_PR2_taxonomy_for_plastid_resolution",
            "resolve_PR2_plastids",
            "import_resolved_PR2_plastid_taxonomy",
            "split_resolved_prok_categories",
        ):
            self.assertIn(f"rule {rule_name}:", prok_rules)

        self.assertIn(
            "silva_taxonomy=rules.export_SILVA_taxonomy_for_plastid_resolution.output",
            prok_rules,
        )
        self.assertIn(
            "pr2_taxonomy=rules.normalize_canonical_16S_PR2_rank_prefixes.output.normalized",
            prok_rules,
        )
        self.assertIn(
            "rules.resolve_PR2_plastids.output.taxonomy_tsv",
            prok_rules,
        )

        resolver = (ROOT / "workflow/scripts/resolve_pr2_plastids.py").read_text(
            encoding="utf-8"
        )
        supplement = (ROOT / "workflow/scripts/export_unmerged_16S.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("subprocess", resolver)
        self.assertNotIn("qiime tools", resolver.lower())
        self.assertNotIn("subprocess", supplement)
        self.assertNotIn("tempfile", supplement)
        self.assertNotIn("qiime tools", supplement.lower())
        self.assertNotIn("biom convert", supplement.lower())
        for rule_name in (
            "export_unmerged_16S_retained_table",
            "export_unmerged_16S_standard_table",
            "export_unmerged_16S_sequences",
            "export_unmerged_16S_SILVA_taxonomy",
            "export_unmerged_16S_PR2_taxonomy",
            "convert_unmerged_16S_retained_biom",
            "convert_unmerged_16S_standard_biom",
            "export_unmerged_16S_supplement",
        ):
            self.assertIn(f"rule {rule_name}:", prok_rules)

        self.assertIn('str_detect(Taxonomy, regex(":plas|', source)
        self.assertIn('plastid_16S_rRNA == "yes" ~ "Chloroplast_16S"', source)

    def test_provenance_is_documented_without_inline_ai_boilerplate(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        ledger_url = (
            "https://github.com/Nwilliams96/515FY-926R-snakemake-NW-edits/"
            "blob/codex/code-provenance/CODE_PROVENANCE.md"
        )
        self.assertIn(ledger_url, readme)
        self.assertFalse((ROOT / "CODE_PROVENANCE.md").exists())
        runtime_files = [
            ROOT / "workflow/rules/common.smk",
            ROOT / "workflow/rules/prepare_classification_dbs.smk",
            ROOT / "workflow/rules/02-denoise-and-export-prok.smk",
            ROOT / "workflow/rules/05-internal-standard-correction.smk",
            ROOT / "workflow/scripts/save_classifier.py",
            ROOT / "workflow/scripts/resolve_pr2_plastids.py",
            ROOT / "workflow/scripts/export_unmerged_16S.py",
            ROOT / "workflow/scripts/assign_silva_species.R",
            ROOT / "workflow/scripts/correct-w-isds-snakemake-v5.R",
            ROOT / "workflow/scripts/correct_16S_18S_ASV-snakemake-v1.R",
            ROOT / "workflow/scripts/generate_pipeline_report.py",
            ROOT / "workflow/scripts/P03-DADA2.sh",
            ROOT / "workflow/scripts/E08-DADA2.sh",
            ROOT / "workflow/scripts/check_reference_capability.py",
            ROOT / "workflow/scripts/database_config.py",
            ROOT / "workflow/scripts/download_verified_reference.py",
            ROOT / "workflow/scripts/normalize_pr2_taxonomy.py",
            ROOT / "run_snakemake_USC_CARC_only.sh",
        ]
        for path in runtime_files:
            with self.subTest(path=path.relative_to(ROOT)):
                source = path.read_text(encoding="utf-8")
                self.assertNotIn("AI-assisted change", source)
                self.assertNotIn("Human review status:", source)
                self.assertNotIn("provenance ledger linked from README.md", source)

    def test_18s_is_classified_once_with_clearly_labelled_pr2_outputs(self):
        rules = (ROOT / "workflow/rules/02-denoise-and-export-euk.smk").read_text(
            encoding="utf-8"
        )
        classifier = (
            ROOT / "workflow/scripts/E10-classify-18S-with-PR2.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("rule classify_18S_with_PR2:", rules)
        self.assertNotIn("rule classify_ASVs_euk:", rules)
        self.assertNotIn("rule euk_PR2_reclassify:", rules)
        self.assertIn("_PR2.classified.qza", rules)
        self.assertIn("_PR2.normalized.classified.qza", rules)
        self.assertIn("12-PR2-barplots", rules)
        self.assertNotIn("SILVA-only", rules)
        self.assertEqual(classifier.count("feature-classifier classify-sklearn"), 1)
        self.assertFalse((ROOT / "workflow/scripts/E14a-PR2-alternative-class.sh").exists())

    def test_new_amplicon_concentration_name_and_isd_ids_reach_outputs(self):
        common = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        merge_rule = (ROOT / "workflow/rules/03-merge-16S-18S.smk").read_text(
            encoding="utf-8"
        )
        correction_rule = (
            ROOT / "workflow/rules/05-internal-standard-correction.smk"
        ).read_text(encoding="utf-8")
        self.assertIn("prok_and_euk_SSU_amplicon_molarities.tsv", common)
        for removed in ("NEW_AMPLICON_CONCENTRATIONS_FILE",
                        "PREVIOUS_AMPLICON_CONCENTRATIONS_FILE",
                        "LEGACY_AMPLICON_CONCENTRATIONS_FILE"):
            self.assertNotIn(removed, common)
        self.assertIn("amplicon_concentrations=AMPLICON_CONCENTRATIONS_FILE", merge_rule)
        self.assertIn("correction_factors=", merge_rule)
        self.assertIn('f"{standard_id}_recovery_ratio"', common)
        self.assertIn('f"mean_{first}_and_{second}_recovery_ratio"', common)
        self.assertIn("corrected=ISD_CORRECTED_LONG_TABLE", correction_rule)
        self.assertIn("filtered=ISD_FILTERED_LONG_TABLE", correction_rule)
        self.assertIn("RESULTS_LONG_DATA = (", common)
        self.assertIn("ISD_FILTERED_LONG_TABLE", common)

    def test_current_molarity_filename_is_selected_even_when_legacy_files_exist(self):
        source = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        start = source.index("\nAMPLICON_CONCENTRATIONS_FILE =") + 1
        block = source[start:source.index("# Keep shared supporting resources", start)]
        with tempfile.TemporaryDirectory() as tmp, chdir(tmp):
            config_dir = Path("config")
            config_dir.mkdir()
            for filename in ("bioanalyzer.tsv", "prok_and_euk_SSU_amplicon_concentrations.tsv"):
                (config_dir / filename).write_text("legacy fixture\n", encoding="utf-8")
            for current_exists in (False, True):
                with self.subTest(current_exists=current_exists):
                    if current_exists:
                        (config_dir / "prok_and_euk_SSU_amplicon_molarities.tsv").touch()
                    scope = {"os": os}
                    exec(compile(block, "common.smk:molarity-file", "exec"), scope)
                    self.assertEqual(
                        scope["AMPLICON_CONCENTRATIONS_FILE"],
                        "config/prok_and_euk_SSU_amplicon_molarities.tsv",
                    )

    def test_reformats_pr2_fasta_without_appending_stale_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "pr2.fasta"
            clean = root / "clean.fasta"
            taxonomy = root / "taxonomy.tsv"
            source.write_text(">taxon one\nACGT\n>taxon two\nTGCA\n", encoding="utf-8")
            clean.write_text("stale\n", encoding="utf-8")
            taxonomy.write_text("stale\n", encoding="utf-8")

            REFORMAT.reformat_pr2_reference(source, clean, taxonomy)

            self.assertEqual(clean.read_text(encoding="utf-8"), ">feature_1\nACGT\n>feature_2\nTGCA\n")
            self.assertEqual(
                taxonomy.read_text(encoding="utf-8"),
                "feature_1\ttaxon one\nfeature_2\ttaxon two\n",
            )

    def test_pr2_prefix_normalization_preserves_nine_ranks_and_plastid_markers(self):
        # RESCRIPt 2026.7 get_pr2.py prepends these handles and joins fields
        # with semicolons without spaces. Exercise the actual rule expression
        # against that serialization rather than merely checking its text.
        rules = (ROOT / "workflow/rules/prepare_classification_dbs.smk").read_text()
        normalizer = rules.split("rule normalize_pr2_reference_taxonomy:", 1)[1].split("rule cull_pr2_seqs:", 1)[0]
        search = re.search(r"--p-search-strings '([^']+)'", normalizer).group(1)
        replacement = re.search(r"--p-replacement-strings '([^']+)'", normalizer).group(1)
        fields = ["Eukaryota", "Archaeplastida", "Chlorophyta:plas", "Chlorophyta_X",
                  "Mamiellophyceae", "Mamiellales", "Mamiellaceae", "Micromonas", "Micromonas_commoda"]
        handles = ["d__", "sgr__", "dv__", "dvs__", "c__", "o__", "f__", "g__", "s__"]
        serialized = ";".join(handle + name for handle, name in zip(handles, fields))
        normalized = re.sub(search, replacement, serialized)
        self.assertEqual(normalized.split(";"), fields)
        self.assertEqual(normalized.count(";"), 8)
        self.assertIn(":plas", normalized)


if __name__ == "__main__":
    unittest.main()
