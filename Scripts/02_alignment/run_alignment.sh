#!/usr/bin/env bash
set -euo pipefail

# Usage:
# run_alignment.sh SAMPLE R1.fastq.gz R2.fastq.gz BOWTIE2_INDEX

sample="$1"
r1="$2"
r2="$3"
index="$4"

threads=8
output_dir="alignment_results"

mkdir -p "${output_dir}"

bowtie2 \
  --very-sensitive \
  --mm \
  --no-mixed \
  --no-discordant \
  --threads "${threads}" \
  -x "${index}" \
  -1 "${r1}" \
  -2 "${r2}" \
  2> "${output_dir}/${sample}.bowtie2.log" \
| samtools view \
    --threads "${threads}" \
    -b \
    -F 4 \
    -q 30 \
    - \
| samtools sort \
    --threads "${threads}" \
    -n \
    -o "${output_dir}/${sample}.namesort.bam" \
    -