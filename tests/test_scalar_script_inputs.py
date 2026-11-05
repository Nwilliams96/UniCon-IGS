"""Execute the real pure-Python consumer rules with completed exporter fixtures."""
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RULES = (ROOT / "workflow/rules/02-denoise-and-export-prok.smk").read_text()


def body(name):
    return RULES.split(f"rule {name}:", 1)[1].split("\nrule ", 1)[0]


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
class ScalarScriptInputTests(unittest.TestCase):
    def run_consumer(self, root, name, script, fixtures, target, settings=""):
        declarations = [settings]
        for rule, path, contents, output_name in fixtures:
            file = root / path
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_text(contents)
            output = (output_name + "=" if output_name else "") + repr(path)
            declarations.append(f'rule {rule}:\n    output:\n        {output}\n    shell:\n        "exit 99"\n')
        consumer = body(name).replace('"../scripts/' + script + '"', repr(str(ROOT / "workflow/scripts" / script)))
        declarations.append(f"rule {name}:" + consumer)
        (root / "Snakefile").write_text("\n".join(declarations))
        process = subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", "Snakefile",
             "--cores", "1", "--nolock", "--notemp", "--shared-fs-usage",
             "input-output", "persistence", "--", target],
            cwd=root, text=True, capture_output=True, timeout=60,
        )
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)

    def rows(self, path):
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_plastid_resolver_receives_scalar_silva_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = "results/02-proks/09-subsetting/tax-merged/fixture.PR2-plastid-resolved-taxonomy.tsv"
            self.run_consumer(root, "resolve_PR2_plastids", "resolve_pr2_plastids.py", [
                ("export_SILVA_taxonomy_for_plastid_resolution", "fixture/SILVA.tsv",
                 "Feature ID\tTaxon\tConfidence\n001\td__Bacteria;p__Test\t0.95\n", ""),
                ("normalize_canonical_16S_PR2_rank_prefixes", "fixture/PR2.tsv",
                 "Feature ID\tTaxon\tConfidence\n001\tEukaryota;Archaeplastida:plas\t0.8000\n", "normalized"),
            ], target, "config = {'projectName': 'fixture'}\nCHLOROPLAST_PR2_MIN_CONFIDENCE = 0.7\n")
            result = self.rows(root / target)
            self.assertEqual(result, [{"Feature ID": "001", "Taxon": "Eukaryota;Archaeplastida:plas", "Confidence": "0.8000"}])
            summary = self.rows(root / "results/02-proks/09-subsetting/tax-merged/fixture.PR2-plastid-routing-summary.tsv")
            self.assertEqual(summary[0]["PR2_rescued_ASVs"], "1")

    def test_supplement_receives_scalar_counts_fasta_and_silva_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.run_consumer(root, "export_unmerged_16S_supplement", "export_unmerged_16S.py", [
                ("convert_unmerged_16S_retained_biom", "fixture/retained.tsv", "#OTU ID\t001\ncanonical\t10\nlinked\t5\n", ""),
                ("convert_unmerged_16S_standard_biom", "fixture/standard.tsv", "#OTU ID\t001\ncanonical\t10\n", ""),
                ("export_unmerged_16S_sequences", "fixture/linked-dna-sequences.fasta", ">linked\nACGTACGT\n", ""),
                ("export_unmerged_16S_SILVA_taxonomy", "fixture/SILVA.tsv", "Feature ID\tTaxon\tConfidence\nlinked\td__Bacteria;p__Test\t0.95\n", ""),
                ("normalize_unmerged_16S_PR2_rank_prefixes", "fixture/PR2.tsv", "Feature ID\tTaxon\tConfidence\nlinked\tEukaryota;Archaeplastida:plas\t0.8\n", "normalized"),
            ], "results/supplement.tsv", 'UNMERGED_16S_TABLE="results/supplement.tsv"\nUNMERGED_16S_SUMMARY="results/summary.tsv"\nCHLOROPLAST_PR2_MIN_CONFIDENCE=0.7\n')
            result = self.rows(root / "results/supplement.tsv")
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0]["Feature_ID"], "linked")
            self.assertEqual(result[0]["Linked_Sequence"], "ACGTACGT")
            self.assertEqual(result[0]["001"], "5")
            summary = self.rows(root / "results/summary.tsv")[0]
            self.assertEqual(summary["standard_non_chimeric_16S_reads"], "10")
            self.assertEqual(summary["supplemental_unmerged_16S_reads"], "5")
            self.assertEqual(summary["supplemental_chloroplast_reads"], "5")


if __name__ == "__main__":
    unittest.main()
