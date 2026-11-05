#!/bin/bash
#
#SBATCH --account=fuhrman_1138
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --job-name=pipeline
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err


set -euo pipefail

read_project_name() {
  awk '
    /^[[:space:]]*projectName[[:space:]]*:/ {
      sub(/^[^:]*:[[:space:]]*/, "")
      gsub(/^["'\'' ]+|["'\'' ]+$/, "")
      print
      exit
    }
  ' config/config.yml
}

read_config_value() {
  local key="$1"
  awk -v key="${key}" '
    $0 ~ "^[[:space:]]*" key "[[:space:]]*:" {
      sub(/^[^:]*:[[:space:]]*/, "")
      gsub(/^["\047 ]+|["\047 ]+$/, "")
      print
      exit
    }
  ' config/config.yml
}

PROJECT_NAME="$(read_project_name)"
if [[ -z "${PROJECT_NAME}" ]]; then
  PROJECT_NAME="pipeline"
fi
JOB_LABEL="${PROJECT_NAME}"
JOB_LABEL="$(printf '%s' "${JOB_LABEL}" | sed 's/[^A-Za-z0-9._-]/_/g')"

# Slurm reads #SBATCH directives before running this file, so it cannot obtain
# projectName from config.yml at that stage. When launched with bash, submit this
# same script with the detected project name. Inside the allocation, continue to
# the workflow normally.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  mkdir -p logs
  echo "Submitting CARC job ${JOB_LABEL}"
  sbatch --job-name="${JOB_LABEL}" "$0" "$@"
  exit $?
fi

# Also rename the queue entry if someone used `sbatch script.sh` directly.
# In that case the log filename retains the generic name because Slurm opened
# it before this update; using `bash script.sh` gives both named jobs and logs.
scontrol update JobId="${SLURM_JOB_ID}" JobName="${JOB_LABEL}" || true

echo "Started at: $(date) on $(hostname)"
echo "Project: ${PROJECT_NAME}; Slurm job: ${JOB_LABEL} (${SLURM_JOB_ID})"
mkdir -p logs

# ZIP files transferred from another computer can preserve timestamps that are
# ahead of the CARC compute-node clock. Snakemake rejects outputs created before
# those future-dated inputs. Reset only affected config files to the current
# compute-node time; normally this block finds nothing and changes nothing.
TIMESTAMP_MARKER="$(mktemp)"
future_config_files=()
while IFS= read -r -d '' config_file; do
  future_config_files+=("${config_file}")
done < <(find config -type f -newer "${TIMESTAMP_MARKER}" -print0)

if (( ${#future_config_files[@]} > 0 )); then
  echo "Correcting future timestamps on ${#future_config_files[@]} config file(s):"
  printf '  %s\n' "${future_config_files[@]}"
  touch -r "${TIMESTAMP_MARKER}" "${future_config_files[@]}"
fi
rm -f "${TIMESTAMP_MARKER}"

# Prevent libraries from taking cores outside the per-rule limits set by
# Snakemake. Multithreaded tools receive their thread counts explicitly.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# CARC assigns each compute job a unique node-local TMPDIR. Keep Conda's
# downloaded-package and repository-metadata cache there so simultaneous jobs
# on different nodes never read or modify the same cache in ~/.conda/pkgs.
# Completed Snakemake rule environments remain in the shared conda_envs_dir.
CONDA_JOB_TMPDIR="${TMPDIR:-/tmp/SLURM_${SLURM_JOB_ID}}"
CONDA_JOB_PKGS_DIR="${CONDA_JOB_TMPDIR%/}/eASV-conda-pkgs"
mkdir -p "${CONDA_JOB_PKGS_DIR}"
export CONDA_PKGS_DIRS="${CONDA_JOB_PKGS_DIR}"
echo "Job-local Conda package cache: ${CONDA_PKGS_DIRS}"

# Ensure conda is available in this shell and activate your snakemake env
# (this uses the conda shims so it works with modern conda installs)
eval "$(conda shell.bash hook)"
conda activate snakemake

echo "Using snakemake: $(which snakemake) ; version: $(snakemake --version || true)"

# Reuse rule-specific environments across project clones. Snakemake identifies
# environments from their definitions and creates only missing/changed ones.
# SNAKEMAKE_CONDA_PREFIX can override the config value for a particular system.
CONDA_ENV_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$(read_config_value conda_envs_dir)}"
if [[ -z "${CONDA_ENV_PREFIX}" ]]; then
  if [[ -d .snakemake/conda ]]; then
    CONDA_ENV_PREFIX=".snakemake"
  else
    CONDA_ENV_PREFIX="../eASV-conda-envs"
  fi
fi
mkdir -p "${CONDA_ENV_PREFIX}"
echo "Rule-specific Conda environments: ${CONDA_ENV_PREFIX}"

# Optional target passed to the script becomes the Snakemake target.
# Example: bash run_snakemake_USC_CARC_only.sh results/02-proks/sample-metadata.tsv
TARGET="${1:-}"
PIPELINE_CORES="${SLURM_CPUS_PER_TASK:-1}"

SNAKEMAKE_ARGS=(
  --snakefile workflow/Snakefile
  --cores "${PIPELINE_CORES}"
  --resources mem_mb=120000
  --use-conda
  --conda-prefix "${CONDA_ENV_PREFIX}"
  --rerun-incomplete
  --latency-wait 60
  --printshellcmds
)

if [[ -n "${TARGET}" ]]; then
  SNAKEMAKE_ARGS+=("${TARGET}")
fi

# Snakemake invocation
# - --cores: use every CPU allocated by Slurm
# - --resources: prevent concurrent rules from exceeding the node's memory
# - --use-conda: create/use conda envs declared in workflow
# - --conda-prefix: reuse matching rule environments across project clones
# - --rerun-incomplete: pick up partial outputs
# - --latency-wait: wait for N seconds when a file appears missing on shared filesystems
if snakemake "${SNAKEMAKE_ARGS[@]}"; then
  rc=0
else
  rc=$?
fi

echo "Snakemake finished with exit code ${rc} at: $(date)"
exit ${rc}
