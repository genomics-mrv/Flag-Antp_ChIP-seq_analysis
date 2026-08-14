#!/usr/bin/env bash
set -euo pipefail

# Generate CPM-normalized bigWig tracks and calculate the Pearson
# correlation between the two Flag-Antp biological replicates.
#
# Usage:
# run_signal_tracks_qc.sh REP1_BAM REP2_BAM
#
# Example:
# scripts/06_signal_tracks_qc/run_signal_tracks_qc.sh \
#   bam_results/Flag-Antp_rep1.bam \
#   bam_results/Flag-Antp_rep2.bam

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 REP1_BAM REP2_BAM" >&2
  exit 1
fi

rep1_bam="$1"
rep2_bam="$2"

threads=8
output_dir="signal_tracks_qc_results"

rep1_bigwig="${output_dir}/Flag-Antp_rep1_CPMnorm.bw"
rep2_bigwig="${output_dir}/Flag-Antp_rep2_CPMnorm.bw"

mkdir -p "${output_dir}"

# Check that the input BAM files exist.

if [[ ! -f "${rep1_bam}" ]]; then
  echo "Error: BAM file not found: ${rep1_bam}" >&2
  exit 1
fi

if [[ ! -f "${rep2_bam}" ]]; then
  echo "Error: BAM file not found: ${rep2_bam}" >&2
  exit 1
fi

# 1. Generate a CPM-normalized bigWig track for replicate 1.

bamCoverage \
  --bam "${rep1_bam}" \
  --outFileName "${rep1_bigwig}" \
  --outFileFormat bigwig \
  --binSize 25 \
  --normalizeUsing CPM \
  --minMappingQuality 30 \
  --smoothLength 75 \
  --skipNonCoveredRegions \
  --ignoreDuplicates \
  --numberOfProcessors "${threads}" \
  2> "${output_dir}/Flag-Antp_rep1_bamCoverage.log"

# 2. Generate a CPM-normalized bigWig track for replicate 2.

bamCoverage \
  --bam "${rep2_bam}" \
  --outFileName "${rep2_bigwig}" \
  --outFileFormat bigwig \
  --binSize 25 \
  --normalizeUsing CPM \
  --minMappingQuality 30 \
  --smoothLength 75 \
  --skipNonCoveredRegions \
  --ignoreDuplicates \
  --numberOfProcessors "${threads}" \
  2> "${output_dir}/Flag-Antp_rep2_bamCoverage.log"

# 3. Summarize the bigWig signal in 10-kb genomic bins.

multiBigwigSummary bins \
  --bwfiles \
    "${rep1_bigwig}" \
    "${rep2_bigwig}" \
  --labels \
    Flag-Antp_rep1 \
    Flag-Antp_rep2 \
  --binSize 10000 \
  --numberOfProcessors "${threads}" \
  --outFileName "${output_dir}/correlation_results.npz" \
  --outRawCounts "${output_dir}/correlation_bin_values.tsv" \
  2> "${output_dir}/multiBigwigSummary.log"

# 4. Calculate the Pearson correlation and generate a scatterplot.

plotCorrelation \
  --corData "${output_dir}/correlation_results.npz" \
  --corMethod pearson \
  --whatToPlot scatterplot \
  --labels \
    Flag-Antp_rep1 \
    Flag-Antp_rep2 \
  --plotTitle "Flag-Antp replicate correlation" \
  --plotFile "${output_dir}/replicate_correlation.svg" \
  --outFileCorMatrix "${output_dir}/Pearson_correlation_matrix.tsv" \
  2> "${output_dir}/plotCorrelation.log"

echo "Signal-track generation and replicate-correlation analysis completed."
echo "Results were written to: ${output_dir}"