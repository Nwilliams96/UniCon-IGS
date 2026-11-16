"""Project-named rule logs and safe CARC submission without a real Slurm job."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "run_snakemake_USC_CARC_only.sh"


class ProjectLogPathTests(unittest.TestCase):
    def test_rule_logs_share_one_project_directory_and_are_ignored_by_git(self):
        common = (ROOT / "workflow/rules/common.smk").read_text()
        self.assertIn('LOG_DIR = PROJECT_NAME + "-logs"', common)
        self.assertLess(common.index('PROJECT_NAME = config["projectName"]'), common.index('LOG_DIR ='))
        paths = []
        for path in (ROOT / "workflow/rules").glob("*.smk"):
            source = path.read_text()
            self.assertNotRegex(source, r'''["']logs/''', str(path))
            paths.extend(re.findall(r'LOG_DIR \+ "(/[^"\n]+)"', source))
        self.assertGreaterEqual(len(paths), 61)
        self.assertIn("/06-report/generate_pipeline_report.log", paths)
        self.assertIn("*-logs/", (ROOT / ".gitignore").read_text().splitlines())

    def submit_fixture(self, root, project, returncode=0):
        root = Path(root)
        (root / "config").mkdir()
        (root / "config/config.yml").write_text(f'projectName: "{project}"\n')
        binaries = root / "bin"
        binaries.mkdir()
        capture = root / "captured.json"
        stub = binaries / "sbatch"
        stub.write_text(
            "#!" + sys.executable + "\nimport json, os, sys\nfrom pathlib import Path\n"
            "args = sys.argv[1:]\n"
            "folder = next(a.split('=',1)[1] for a in args if a.startswith('--output=')).split('/')[0]\n"
            "Path(os.environ['LOG_TEST_CAPTURE']).write_text(json.dumps({'args': args, 'directory_exists': Path(folder).is_dir()}))\n"
            "raise SystemExit(int(os.environ['LOG_TEST_EXIT_CODE']))\n"
        )
        stub.chmod(0o700)
        env = dict(os.environ, PATH=str(binaries) + os.pathsep + os.environ.get("PATH", ""),
                   LOG_TEST_CAPTURE=str(capture), LOG_TEST_EXIT_CODE=str(returncode))
        env.pop("SLURM_JOB_ID", None)  # Always test the submission branch, not a live allocation.
        target = "results/a path/sample-metadata.tsv"
        process = subprocess.run(["bash", str(RUNNER), target], cwd=root, env=env,
                                 capture_output=True, text=True, timeout=15)
        return process, json.loads(capture.read_text()), target

    def test_bash_submission_creates_named_logs_first_and_preserves_target_arguments(self):
        for project in ("AMT29", "Test_Project.1"):
            with self.subTest(project=project), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                old = root / "logs/previous.err"
                old.parent.mkdir()
                old.write_text("Historical fixture log\n")
                before = (old.stat().st_mtime_ns, old.read_bytes())
                process, captured, target = self.submit_fixture(root, project)
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                self.assertEqual(captured["args"], [f"--job-name={project}",
                                 f"--output={project}-logs/%x_%j.out", f"--error={project}-logs/%x_%j.err",
                                 str(RUNNER), target])
                self.assertTrue(captured["directory_exists"])
                self.assertEqual((old.stat().st_mtime_ns, old.read_bytes()), before)
                self.assertNotIn("Started at:", process.stdout)  # No compute job ran.

    def test_submission_failure_is_returned_without_starting_the_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            process, captured, _ = self.submit_fixture(tmp, "Failed-submission", returncode=27)
            self.assertEqual(process.returncode, 27)
            self.assertTrue(captured["directory_exists"])
            self.assertNotIn("Started at:", process.stdout)


if __name__ == "__main__":
    unittest.main()
