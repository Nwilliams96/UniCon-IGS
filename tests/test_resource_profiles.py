"""Resource policy and harmless launcher/scheduler fixtures; no Slurm or Conda jobs."""

import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "workflow/profiles/resources"
CLASSIFIERS = {
    "classify_ASVs", "classify_unmerged_16S_with_SILVA", "classify_all_16S_with_PR2",
    "classify_unmerged_16S_with_PR2", "classify_18S_with_PR2",
}
HEAVY_RULES = CLASSIFIERS | {
    "prepare_bbsplit_db", "bbsplit_prok_euk", "denoise_prok_dada2",
    "denoise_prok_dada2_with_unmerged", "denoise_euk_dada2", "assign_SILVA_species",
    "train_classifier", "train_classifier_pr2",
}
MEMORY = {
    "download_prok_db": 1024, "download_euk_db": 1024,
    "prepare_bbsplit_db": 64000, "train_classifier": 110000,
    "train_classifier_pr2": 32000, "ensure_pr2_classifier": 2000,
    "bbsplit_prok_euk": 64000, "deinterleave_split_reads_euk": 4000,
    "deinterleave_split_reads_prok": 4000, "bbduk_cut_reads": 2500,
    "fuse_trimmed_euk_seqs": 2500, "denoise_prok_dada2": 96000,
    "denoise_prok_dada2_with_unmerged": 96000, "denoise_euk_dada2": 96000,
    "classify_ASVs": 64000, "classify_unmerged_16S_with_SILVA": 64000,
    "classify_all_16S_with_PR2": 110000, "classify_unmerged_16S_with_PR2": 110000,
    "classify_18S_with_PR2": 110000, "assign_SILVA_species": 110000,
}


class ResourceProfileTests(unittest.TestCase):
    def test_profile_keeps_explicit_estimates_and_caps_classifier_workers(self):
        profile = yaml.safe_load((PROFILE / "config.yaml").read_text())
        self.assertEqual(profile["resources"], {"heavy_step": 1})
        self.assertEqual(profile["default-resources"], {"mem_mb": 2000})
        self.assertEqual(profile["set-threads"], {name: 1 for name in CLASSIFIERS})
        resources = profile["set-resources"]
        self.assertEqual({name: values["mem_mb"] for name, values in resources.items()}, MEMORY)
        self.assertEqual({name for name, values in resources.items() if values.get("heavy_step") == 1},
                         HEAVY_RULES)
        self.assertEqual({name: values["runtime"] for name, values in resources.items() if "runtime" in values},
                         {"download_prok_db": 30, "download_euk_db": 30,
                          "prepare_bbsplit_db": 240, "train_classifier": 720,
                          "train_classifier_pr2": 360, "ensure_pr2_classifier": 30,
                          "assign_SILVA_species": 360})
        source = "\n".join(p.read_text() for p in (ROOT / "workflow/rules").glob("*.smk"))
        available = set(re.findall(r"^\s*rule\s+(\w+):", source, re.M))
        self.assertLessEqual(set(resources) | set(profile["set-threads"]), available)
        # Keep the public profile visible without exposing downloaded resources.
        ignore = (ROOT / ".gitignore").read_text().splitlines()
        self.assertIn("resources/", ignore)
        self.assertIn("!/workflow/profiles/resources/", ignore)
        self.assertIn("/workflow/profiles/resources/*", ignore)
        self.assertIn("!/workflow/profiles/resources/config.yaml", ignore)

    def run_launcher_fixture(self, root, *, local=False, allocation=None, custom_profile=None,
                             arguments=()):
        """Run the real launcher against executable stubs confined to the fixture."""
        root = Path(root)
        (root / "config").mkdir()
        shared = root / "shared rule envs"
        (root / "config/config.yml").write_text(
            f'projectName: "Resource_Test.1"\nconda_envs_dir: "{shared}"\n'
        )
        binaries = root / "bin"
        binaries.mkdir()
        capture = root / "captured.json"
        (binaries / "conda").write_text(
            "#!/usr/bin/env bash\n"
            "if [[ \"$1\" == shell.bash && \"${2:-}\" == hook ]]; then\n"
            "  printf 'conda() { :; }\\n'\nfi\n"
        )
        (binaries / "scontrol").write_text("#!/usr/bin/env bash\nexit 0\n")
        (binaries / "snakemake").write_text(
            "#!" + sys.executable + "\nimport json, os, sys\nfrom pathlib import Path\n"
            "if sys.argv[1:] == ['--version']:\n    print('fixture')\n    raise SystemExit(0)\n"
            "Path(os.environ['RESOURCE_TEST_CAPTURE']).write_text(json.dumps({\n"
            "    'args': sys.argv[1:], 'RQCMEM': os.environ.get('RQCMEM'),\n"
            "    'CONDA_PKGS_DIRS': os.environ.get('CONDA_PKGS_DIRS')}))\n"
        )
        for binary in binaries.iterdir():
            binary.chmod(0o700)
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        for key in ("SNAKEMAKE_CONDA_PREFIX", "EASV_RESOURCE_PROFILE", "RQCMEM"):
            env.pop(key, None)
        tmpdir = root / "node tmp"
        tmpdir.mkdir()
        env.update(PATH=str(binaries) + os.pathsep + env.get("PATH", ""),
                   TMPDIR=str(tmpdir), RESOURCE_TEST_CAPTURE=str(capture))
        if not local:
            env.update(SLURM_JOB_ID="fixture-job", SLURM_CPUS_PER_TASK="4")
            env.update(allocation or {})
        if custom_profile is not None:
            env["EASV_RESOURCE_PROFILE"] = str(custom_profile)
        launcher = ROOT / ("run_snakemake.sh" if local else "run_snakemake_USC_CARC_only.sh")
        process = subprocess.run(["bash", str(launcher), *arguments], cwd=root, env=env,
                                 capture_output=True, text=True, timeout=15)
        return process, json.loads(capture.read_text()) if capture.exists() else None

    def test_carc_budget_uses_allocated_memory_and_keeps_target_intact(self):
        cases = [
            ({"SLURM_MEM_PER_NODE": "131072"}, 131072),
            ({"SLURM_MEM_PER_NODE": "65536", "SLURM_CPUS_PER_TASK": "2"}, 65536),
            ({"SLURM_MEM_PER_CPU": "8192"}, 8192 * 4),
            ({"SLURM_MEM_PER_NODE": "65536", "SLURM_MEM_PER_CPU": "1024"}, 65536),
            ({"SLURM_MEM_PER_NODE": "0", "SLURM_MEM_PER_CPU": "2048"}, 2048 * 4),
            ({"SLURM_MEM_PER_NODE": "invalid", "SLURM_MEM_PER_CPU": "2048"}, 2048 * 4),
        ]
        for allocation, memory_mib in cases:
            with self.subTest(allocation=allocation), tempfile.TemporaryDirectory() as tmp:
                target = "results/a path/sample-metadata.tsv"
                process, captured = self.run_launcher_fixture(tmp, allocation=allocation,
                                                              arguments=(target,))
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                args = captured["args"]
                budget = (memory_mib * 1048576 // 1000000) * 90 // 100
                self.assertEqual(args[args.index("--cores") + 1], allocation.get("SLURM_CPUS_PER_TASK", "4"))
                self.assertEqual(args[args.index("--resources") + 1:args.index("--resources") + 3],
                                 [f"mem_mb={budget}", "heavy_step=1"])
                self.assertEqual(args[args.index("--workflow-profile") + 1], "workflow/profiles/resources")
                self.assertEqual(captured["RQCMEM"], str(budget * 1000000 // 1048576))
                self.assertEqual(args[-1], target)
                self.assertTrue(Path(captured["CONDA_PKGS_DIRS"]).is_relative_to(Path(tmp)))

    def test_carc_rejects_missing_zero_or_invalid_allocation_memory(self):
        for allocation in ({}, {"SLURM_MEM_PER_NODE": "0"}, {"SLURM_MEM_PER_NODE": "-1"},
                           {"SLURM_MEM_PER_NODE": "64G"}, {"SLURM_MEM_PER_CPU": "0"},
                           {"SLURM_MEM_PER_NODE": "bad", "SLURM_MEM_PER_CPU": "bad"}):
            with self.subTest(allocation=allocation), tempfile.TemporaryDirectory() as tmp:
                process, captured = self.run_launcher_fixture(tmp, allocation=allocation)
                self.assertEqual(process.returncode, 2, process.stdout + process.stderr)
                self.assertIsNone(captured)
                self.assertIn("Specify a positive Slurm --mem or --mem-per-cpu allocation.", process.stderr)

    def test_both_launchers_preserve_custom_profile_paths_with_spaces(self):
        for local in (False, True):
            with self.subTest(local=local), tempfile.TemporaryDirectory() as tmp:
                profile = Path(tmp) / "custom resource profile"
                process, captured = self.run_launcher_fixture(
                    tmp, local=local, allocation={"SLURM_MEM_PER_NODE": "65536"},
                    custom_profile=profile,
                )
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                args = captured["args"]
                self.assertEqual(args[args.index("--workflow-profile") + 1], str(profile))

    def test_local_launcher_retains_one_core_default_and_all_user_arguments(self):
        with tempfile.TemporaryDirectory() as tmp:
            arguments = ["--dry-run", "--cores", "2", "--resources", "mem_mb=16000", "heavy_step=1",
                         "results/a path/report.html"]
            process, captured = self.run_launcher_fixture(tmp, local=True, arguments=arguments)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            args = captured["args"]
            self.assertEqual(args[args.index("--cores") + 1], "1")
            self.assertEqual(args[-len(arguments):], arguments)
            self.assertEqual(args[args.index("--workflow-profile") + 1], "workflow/profiles/resources")
            self.assertIsNone(captured["RQCMEM"])  # No new local hardware inference was approved.


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
class NativeResourceSchedulingTests(unittest.TestCase):
    def test_profile_serializes_heavy_rules_but_permits_light_work(self):
        """Execute only sleeping synthetic jobs, with a two-heavy positive control."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(PROFILE, root / "profile with spaces")
            worker = root / "worker.py"
            worker.write_text(
                "import json, sys, time\nfrom pathlib import Path\n"
                "start = time.monotonic_ns()\n"
                "time.sleep(1.5 if sys.argv[1] != 'light' else 3.0)\n"
                "Path(sys.argv[2]).write_text(json.dumps({'start': start, 'end': time.monotonic_ns(), 'threads': int(sys.argv[3])}))\n"
            )
            jobs = {"classify_ASVs": "first.json", "classify_unmerged_16S_with_SILVA": "second.json",
                    "light": "light.json"}
            snakefile = "rule all:\n    input: " + repr(list(jobs.values())) + "\n"
            for name, output in jobs.items():
                snakefile += (f"\nrule {name}:\n    output: {output!r}\n"
                              f"    threads: {1 if name == 'light' else 2}\n"
                              "    resources: mem_mb=2000\n"
                              f"    shell: {str(sys.executable + ' worker.py ' + name + ' {output:q} {threads}')!r}\n")
            (root / "Snakefile").write_text(snakefile)
            intervals = {}
            for capacity in (1, 2):
                process = subprocess.run(
                    [sys.executable, "-m", "snakemake", "--snakefile", "Snakefile", "--cores", "3",
                     "--workflow-profile", str(root / "profile with spaces"), "--resources", "mem_mb=150000",
                     f"heavy_step={capacity}", "--forceall", "--nolock", "--shared-fs-usage",
                     "input-output", "persistence"], cwd=root,
                    capture_output=True, text=True, timeout=45,
                )
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                intervals[capacity] = {name: json.loads((root / output).read_text())
                                       for name, output in jobs.items()}
            def overlap(first, second):
                return max(first["start"], second["start"]) < min(first["end"], second["end"])
            limited = intervals[1]
            self.assertEqual(limited["classify_ASVs"]["threads"], 1)
            self.assertEqual(limited["classify_unmerged_16S_with_SILVA"]["threads"], 1)
            self.assertFalse(overlap(limited["classify_ASVs"], limited["classify_unmerged_16S_with_SILVA"]))
            self.assertTrue(any(overlap(limited["light"], limited[name]) for name in CLASSIFIERS if name in limited))
            control = intervals[2]
            self.assertTrue(overlap(control["classify_ASVs"], control["classify_unmerged_16S_with_SILVA"]))

    def test_actual_workflow_dry_run_uses_profile_without_changing_scientific_rules(self):
        # Reuse the existing synthetic database/config builder without duplicating its formats.
        specification = importlib.util.spec_from_file_location("profile_dag_fixture", ROOT / "tests/test_database_dag.py")
        helpers = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(helpers)
        fixture = helpers.DatabaseDagTests(methodName="runTest")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = fixture.setup_project(root)
            for database in ("SILVA", "PR2"):
                fixture.set_mode(root, config, database, "existing")
            raw = Path(config["rawdatadir"])
            raw.mkdir()
            for mate in ("R1", "R2"):
                (raw / f"001.{mate}.fastq.gz").touch()
            (root / "config/config.yml").write_text(json.dumps(config))
            process = subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", "workflow/Snakefile", "--dry-run",
                 "--forceall", "--nolock", "--shared-fs-usage", "input-output", "persistence",
                 "--cores", "4", "--workflow-profile", "workflow/profiles/resources",
                 "--resources", "mem_mb=60000", "heavy_step=1", "--printshellcmds", "all"],
                cwd=root, capture_output=True, text=True, timeout=60,
            )
            output = process.stdout + process.stderr
            self.assertEqual(process.returncode, 0, output)
            for name in HEAVY_RULES - {"prepare_bbsplit_db", "train_classifier", "train_classifier_pr2"}:
                match = re.search(r"(?:local)?rule " + re.escape(name) + r":(.*?)(?=\n\[[^\n]+\]\n|\Z)",
                                  output, re.S)
                self.assertIsNotNone(match, output)
                self.assertRegex(match.group(1), r"resources: [^\n]*heavy_step=1")
                if name in CLASSIFIERS:
                    # Snakemake omits the thread field when the effective count is one.
                    self.assertNotRegex(match.group(1), r"(?m)^\s+threads: (?:[2-9]|[1-9][0-9]+)\s*$")
                if name.startswith("denoise_"):
                    self.assertRegex(match.group(1), r"(?m)^\s+threads: 4\s*$")


if __name__ == "__main__":
    unittest.main()
