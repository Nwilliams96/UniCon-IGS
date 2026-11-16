"""PR2 rank notation must not make its lineages look like SILVA downstream."""

import csv
import importlib.util
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
PROK_RULES = (ROOT / "workflow/rules/02-denoise-and-export-prok.smk").read_text()
EUK_RULES = (ROOT / "workflow/rules/02-denoise-and-export-euk.smk").read_text()
SCRIPT = ROOT / "workflow/scripts/normalize_pr2_taxonomy.py"
SPEC = importlib.util.spec_from_file_location("normalize_pr2_taxonomy", SCRIPT)
NORMALIZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NORMALIZER)
RESOLVER_SPEC = importlib.util.spec_from_file_location(
    "resolve_pr2_plastids_for_normalization", ROOT / "workflow/scripts/resolve_pr2_plastids.py"
)
RESOLVER = importlib.util.module_from_spec(RESOLVER_SPEC)
RESOLVER_SPEC.loader.exec_module(RESOLVER)


def rule_body(source, name):
    return source.split(f"rule {name}:", 1)[1].split("\nrule ", 1)[0]


class Pr2RankNormalizationTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
    def test_real_snakemake_entrypoint_receives_one_path_in_all_three_branches(self):
        for source, normalization, exported, raw_path in (
            (PROK_RULES, "normalize_unmerged_16S_PR2_rank_prefixes", "export_unmerged_16S_PR2_taxonomy",
             ".snakemake/classification-exports/unmerged-16S/PR2-raw/taxonomy.tsv"),
            (PROK_RULES, "normalize_canonical_16S_PR2_rank_prefixes", "export_PR2_taxonomy_for_plastid_resolution",
             ".snakemake/classification-exports/canonical-16S/PR2-raw/taxonomy.tsv"),
            (EUK_RULES, "normalize_18S_PR2_rank_prefixes", "export_18S_PR2_taxonomy_for_normalization",
             ".snakemake/classification-exports/18S/PR2-raw/taxonomy.tsv"),
        ):
            with self.subTest(branch=normalization), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                raw = root / raw_path
                raw.parent.mkdir(parents=True)
                raw.write_text(
                    "Feature ID\tTaxon\tConfidence\n"
                    "001\td__Eukaryota;sgr__Archaeplastida:plas\t0.8000\n"
                )
                body = rule_body(source, normalization).replace(
                    '"../scripts/normalize_pr2_taxonomy.py"', repr(str(SCRIPT))
                )
                target = re.search(r'normalized=temp\("([^"]+)"\)', body)[1]
                # The export fixture is already complete. Only the actual
                # normalization rule and its real script may execute here.
                (root / "Snakefile").write_text(
                    # This isolated rule fixture does not include common.smk.
                    'LOG_DIR = "fixture-logs"\n\n'
                    f'rule {exported}:\n    output:\n        {raw_path!r}\n'
                    '    shell:\n        "exit 99"\n\n'
                    f'rule {normalization}:' + body
                )
                process = subprocess.run(
                    [sys.executable, "-m", "snakemake", "--snakefile", "Snakefile",
                     "--cores", "1", "--nolock", "--notemp", "--shared-fs-usage",
                     "input-output", "persistence", "--", target],
                    cwd=root, text=True, capture_output=True, timeout=60,
                )
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                self.assertEqual((root / target).read_text(),
                                 "Feature ID\tTaxon\tConfidence\n"
                                 "001\tEukaryota;Archaeplastida:plas\t0.8000\n")

    def test_normalization_preserves_confidence_ids_and_extra_columns(self):
        # Nine-rank fixture includes an organelle marker and an internal d__
        # substring that is part of a name, not an actual rank prefix.
        plain = (
            "Eukaryota; Archaeplastida; Chlorophyta:plas; Chlorophyta_X; "
            "Class_d__name; Order_name; Family_name; Genus_name; Species_name"
        )
        prefixed = (
            "d__Eukaryota; sgr__Archaeplastida; dv__Chlorophyta:plas; "
            "dvs__Chlorophyta_X; c__Class_d__name; o__Order_name; "
            "f__Family_name; g__Genus_name; s__Species_name"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, target, repeated = (root / name for name in ("raw.tsv", "normalized.tsv", "twice.tsv"))
            header = ["Feature ID", "Taxon", "Confidence", "Additional column"]
            rows = [
                ["001", prefixed, "0.8", "verbatim"],
                ["002", plain, "0.8000", ""],
                ["003", "Unassigned", "0.000", "unchanged"],
            ]
            with source.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)
            self.assertEqual(NORMALIZER.normalize_pr2_taxonomy(source, target), 1)
            with target.open(newline="") as handle:
                observed = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(observed[0], header)
            expected = [row[:] for row in rows]
            expected[0][1] = plain
            self.assertEqual(observed[1:], expected)
            self.assertEqual(len(observed[1][1].split(";")), 9)
            self.assertEqual(observed[1][2], "0.8")
            self.assertEqual(observed[2][2], "0.8000")
            self.assertEqual(NORMALIZER.normalize_pr2_taxonomy(target, repeated), 0)
            self.assertEqual(target.read_bytes(), repeated.read_bytes())

    def test_invalid_taxonomy_header_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "raw.tsv", Path(tmp) / "out.tsv"
            source.write_text("Feature ID\tConfidence\n001\t0.8\n")
            with self.assertRaisesRegex(ValueError, "Taxon or taxonomy"):
                NORMALIZER.normalize_pr2_taxonomy(source, target)

    def test_confidence_survives_normalization_for_plastid_decisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, normalized = root / "raw.tsv", root / "normalized.tsv"
            source.write_text(
                "Feature ID\tTaxon\tConfidence\n"
                "001\td__Eukaryota;sgr__Archaeplastida:plas\t0.8\n"
                "002\td__Eukaryota;sgr__Archaeplastida:plas\t0.69\n"
            )
            NORMALIZER.normalize_pr2_taxonomy(source, normalized)
            pr2_rows = RESOLVER.read_taxonomy(normalized)
            silva = [{
                "feature_id": feature_id, "taxonomy": "d__Bacteria;p__Test",
                "confidence": 0.95, "confidence_text": "0.95",
            } for feature_id in ("001", "002")]
            final, audit, summary = RESOLVER.resolve_taxonomies(silva, pr2_rows, 0.7)
            self.assertEqual(final[0]["Taxon"], "Eukaryota;Archaeplastida:plas")
            self.assertEqual(final[0]["Confidence"], "0.8")
            self.assertEqual(final[1]["Taxon"], "d__Bacteria;p__Test")
            self.assertEqual(audit[1]["PR2_confidence"], "0.69")
            self.assertEqual(summary["PR2_rescued_ASVs"], 1)
            self.assertEqual(summary["PR2_low_confidence_not_accepted_ASVs"], 1)

    def test_all_consumers_receive_normalized_taxonomy(self):
        for source, consumer, normalized in (
            (PROK_RULES, "resolve_PR2_plastids", "normalize_canonical_16S_PR2_rank_prefixes"),
            (PROK_RULES, "export_unmerged_16S_supplement", "normalize_unmerged_16S_PR2_rank_prefixes"),
            (EUK_RULES, "make_PR2_euk_barplots", "import_normalized_18S_PR2_taxonomy"),
            (EUK_RULES, "export_tax_convert_biom_euk", "import_normalized_18S_PR2_taxonomy"),
        ):
            with self.subTest(consumer=consumer):
                self.assertIn(
                    f"rules.{normalized}.output.normalized",
                    rule_body(source, consumer),
                )
        # Source acquisition has no effect on these rules: the raw classification
        # always feeds normalization, so all three acquisition methods are covered.
        for source, normalization, exported in (
            (PROK_RULES, "normalize_canonical_16S_PR2_rank_prefixes", "export_PR2_taxonomy_for_plastid_resolution"),
            (PROK_RULES, "normalize_unmerged_16S_PR2_rank_prefixes", "export_unmerged_16S_PR2_taxonomy"),
            (EUK_RULES, "normalize_18S_PR2_rank_prefixes", "export_18S_PR2_taxonomy_for_normalization"),
        ):
            body = rule_body(source, normalization)
            self.assertIn(f"taxonomy=rules.{exported}.output[0]", body)
            self.assertIn('"../scripts/normalize_pr2_taxonomy.py"', body)
            self.assertNotIn("qiime rescript edit-taxonomy", body)
        self.assertIn("rules.classify_all_16S_with_PR2.output.classified", rule_body(PROK_RULES, "export_PR2_taxonomy_for_plastid_resolution"))
        self.assertIn("rules.classify_unmerged_16S_with_PR2.output.classified", rule_body(PROK_RULES, "export_unmerged_16S_PR2_taxonomy"))
        self.assertIn("rules.classify_18S_with_PR2.output.normalized", rule_body(EUK_RULES, "export_18S_PR2_taxonomy_for_normalization"))
        self.assertIn("rules.normalize_18S_PR2_rank_prefixes.output.normalized", rule_body(EUK_RULES, "import_normalized_18S_PR2_taxonomy"))
        self.assertIn(
            '_PR2.normalized.classified.qza',
            rule_body(EUK_RULES, "import_normalized_18S_PR2_taxonomy"),
        )


if __name__ == "__main__":
    unittest.main()
