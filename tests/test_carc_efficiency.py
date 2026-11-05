import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class CarcEfficiencyTests(unittest.TestCase):
    def test_runner_uses_slurm_cpu_and_memory_allocations(self):
        runner = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text(encoding="utf-8")
        self.assertIn("#SBATCH --cpus-per-task=8", runner)
        self.assertIn('PIPELINE_CORES="${SLURM_CPUS_PER_TASK:-1}"', runner)
        self.assertIn('--cores "${PIPELINE_CORES}"', runner)
        self.assertIn("--resources mem_mb=120000", runner)

    def test_runner_labels_jobs_from_project_name(self):
        runner = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text(encoding="utf-8")
        self.assertIn("read_project_name()", runner)
        self.assertIn('JOB_LABEL="${PROJECT_NAME}"', runner)
        self.assertIn('sbatch --job-name="${JOB_LABEL}" "$0" "$@"', runner)
        self.assertIn("#SBATCH --output=logs/%x_%j.out", runner)
        self.assertIn("#SBATCH --error=logs/%x_%j.err", runner)

    def test_runner_repairs_only_future_dated_config_files(self):
        runner = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text(encoding="utf-8")
        self.assertIn('find config -type f -newer "${TIMESTAMP_MARKER}" -print0', runner)
        self.assertIn('touch -r "${TIMESTAMP_MARKER}" "${future_config_files[@]}"', runner)

    def test_bbsplit_receives_threads_and_realistic_memory(self):
        rules = (ROOT / "workflow/rules/01-split-16S-18S.smk").read_text(
            encoding="utf-8"
        )
        self.assertIn("rule bbsplit_prok_euk:", rules)
        self.assertIn("threads: 8", rules)
        self.assertIn("mem_mb=64000", rules)
        self.assertIn("bbsplit.sh threads={threads}", rules)

    def test_split_counts_use_fastq_records_not_header_like_quality_lines(self):
        script = (ROOT / "workflow/scripts/count-seqs.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("awk 'END { print NR / 4 }'", script)
        self.assertNotIn('zgrep -c "^@"', script)

    def test_dada2_and_classifiers_receive_rule_threads(self):
        for relative_path in (
            "workflow/scripts/P03-DADA2.sh",
            "workflow/scripts/E08-DADA2.sh",
        ):
            script = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("--p-n-threads ${snakemake[threads]}", script)

        for relative_path in (
            "workflow/scripts/P05-classify-eASVs.sh",
            "workflow/scripts/E10-classify-18S-with-PR2.sh",
        ):
            script = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("--p-n-jobs", script)
            self.assertIn("${snakemake[threads]}", script)

    def test_large_pr2_classifications_limit_worker_memory(self):
        rules = (ROOT / "workflow/rules/02-denoise-and-export-prok.smk").read_text(
            encoding="utf-8"
        )
        for rule_name in (
            "classify_unmerged_16S_with_PR2",
            "classify_all_16S_with_PR2",
        ):
            rule_block = rules.split(f"rule {rule_name}:", 1)[1].split("\nrule ", 1)[0]
            self.assertIn("threads: 2", rule_block)
            self.assertIn("mem_mb=110000", rule_block)


if __name__ == "__main__":
    unittest.main()
