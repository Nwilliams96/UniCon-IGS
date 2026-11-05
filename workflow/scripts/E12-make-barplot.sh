#!/usr/bin/env bash

set -euo pipefail

timestamp=$(date +"%y%m%d-%H%M")

qiime taxa barplot \
  --i-table "${snakemake_input[euktable]}" \
  --i-taxonomy "${snakemake_input[euktax]}" \
  --m-metadata-file "${snakemake_input[eukmetadata]}" \
  --output-dir "${snakemake_output[0]}"

cp \
  "${snakemake_output[0]}/visualization.qzv" \
  "${snakemake_output[0]}/${timestamp}.${snakemake_params[projectName]}.18S.PR2.barplot.qzv"
