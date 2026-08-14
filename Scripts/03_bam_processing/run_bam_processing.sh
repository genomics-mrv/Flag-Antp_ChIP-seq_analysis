#!/usr/bin/env bash
set -euo pipefail

# Usage:
# run_bam_processing.sh SAMPLE INPUT_NAMESORT_BAM

sample="$1"
input_bam="$2"

threads=8
output_dir="bam_results"

mkdir -p "${output_dir}"

samtools fixmate \
  --threads "${threads}" \
  -m \
  "${input_bam}" \
  "${output_dir}/${sample}.fixmate.bam"

samtools sort \
  --threads "${threads}" \
  -o "${output_dir}/${sample}.coordsort.bam" \
  "${output_dir}/${sample}.fixmate.bam"

samtools markdup \
  --threads "${threads}" \
  -r \
  "${output_dir}/${sample}.coordsort.bam" \
  "${output_dir}/${sample}.bam"

samtools index \
  --threads "${threads}" \
  "${output_dir}/${sample}.bam"

samtools flagstat \
  --threads "${threads}" \
  "${output_dir}/${sample}.bam" \
  > "${output_dir}/${sample}.flagstat.txt"