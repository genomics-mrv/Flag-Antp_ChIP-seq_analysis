#!/usr/bin/env bash
set -euo pipefail

mkdir -p fastqc_results

fastqc \
  --threads 8 \
  --outdir fastqc_results \
  *.fastq.gz