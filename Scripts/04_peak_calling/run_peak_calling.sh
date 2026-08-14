#!/usr/bin/env bash
set -euo pipefail

# Usage:
# run_macs2.sh SAMPLE CHIP_BAM INPUT_BAM

sample="$1"
chip_bam="$2"
input_bam="$3"

output_dir="peak_calling_results/${sample}"

mkdir -p "${output_dir}"

macs2 callpeak \
  --treatment "${chip_bam}" \
  --control "${input_bam}" \
  --format BAMPE \
  --gsize dm \
  --name "${sample}" \
  --pvalue 1e-4 \
  --outdir "${output_dir}"