"""The generated projectName is the sole naming input across workflow consumers."""
from pathlib import Path
import runpy
import subprocess
import tempfile
from types import SimpleNamespace
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ProjectNamingTests(unittest.TestCase):
    def test_workflow_has_no_legacy_name_or_helper_references(self):
        paths = [*ROOT.joinpath("workflow/rules").glob("*.smk"),
                 *ROOT.joinpath("workflow/scripts").glob("*"),
                 ROOT / "run_snakemake_USC_CARC_only.sh"]
        for path in paths:
            if not path.is_file():
                continue
            with self.subTest(path=path.name):
                text = path.read_text()
                self.assertNotIn("studyName", text)
                self.assertNotIn("publish_pr2_classifier.py", text)
                self.assertNotIn("publish_classifier(", text)

    def test_carc_name_reader_accepts_only_the_project_name(self):
        source = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text()
        function = source.split("read_project_name() {", 1)[1].split("\n}", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            config = Path(tmp) / "config"
            config.mkdir()
            for name in ("001", "AMT29-DNA", "Test_Project.1"):
                (config / "config.yml").write_text(f'projectName: "{name}"\n')
                process = subprocess.run(
                    ["bash", "-c", "read_project_name() {" + function + "\n}\nread_project_name"],
                    cwd=tmp, capture_output=True, text=True, check=True,
                )
                self.assertEqual(process.stdout.strip(), name)

    def test_export_script_needs_only_project_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv"
            source.write_text("SampleID\tPhylum\tOrder\tCorrected_Sequence_Counts\n001\tP\tO\t2\n")
            export = root / "Example-Results-Export"
            destination = export / "Example.long_data.tsv"
            readme = export / "README.txt"
            snakemake = SimpleNamespace(
                config={"projectName": "Example"},
                input=SimpleNamespace(data_files=[source]),
                output=SimpleNamespace(data_files=[destination], readme=readme, summaries=[]),
            )
            runpy.run_path(str(ROOT / "workflow/scripts/export_results.py"),
                           init_globals={"snakemake": snakemake})
            self.assertEqual(source.read_bytes(), destination.read_bytes())
            self.assertIn("Results export for Example", readme.read_text())


if __name__ == "__main__":
    unittest.main()
