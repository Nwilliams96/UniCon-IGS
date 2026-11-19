#!/usr/bin/env bash

set -euo pipefail

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

# Keep rule environments outside individual analysis clones so Snakemake can
# reuse an existing environment whenever its environment definition is unchanged.
# SNAKEMAKE_CONDA_PREFIX can override the config value for a particular system.
CONDA_ENV_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$(read_config_value conda_envs_dir)}"
if [[ -z "${CONDA_ENV_PREFIX}" ]]; then
  if [[ -d .snakemake/conda ]]; then
    # Preserve reuse for projects created before conda_envs_dir was introduced.
    CONDA_ENV_PREFIX=".snakemake"
  else
    CONDA_ENV_PREFIX="../eASV-conda-envs"
  fi
fi

mkdir -p "${CONDA_ENV_PREFIX}"
echo "Rule-specific Conda environments: ${CONDA_ENV_PREFIX}"

# Select the same editable resource policy for local execution.
RESOURCE_PROFILE="${EASV_RESOURCE_PROFILE:-workflow/profiles/resources}"

snakemake \
  --snakefile workflow/Snakefile \
  --workflow-profile "${RESOURCE_PROFILE}" \
  --use-conda \
  --conda-prefix "${CONDA_ENV_PREFIX}" \
  --cores 1 \
  "$@"
