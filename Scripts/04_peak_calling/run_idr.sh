#!/usr/bin/env bash
set -euo pipefail

# Usage:
# run_idr.sh REP1_PEAKS REP2_PEAKS OUTPUT_PREFIX

rep1_peaks="$1"
rep2_peaks="$2"
output_prefix="$3"

output_dir="idr_results"

mkdir -p "${output_dir}"

idr \
  --samples "${rep1_peaks}" "${rep2_peaks}" \
  --input-file-type narrowPeak \
  --output-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --soft-idr-threshold 0.05 \
  --output-file "${output_dir}/${output_prefix}_IDR0.05.narrowPeak" \
  --plot \
  2> "${output_dir}/${output_prefix}_IDR.log"