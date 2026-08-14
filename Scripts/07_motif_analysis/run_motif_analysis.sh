#!/usr/bin/env bash
set -euo pipefail

# Extract 200-bp sequences centered on IDR peak summits and scan
# them for Antp motifs using FIMO.
#
# Usage:
# run_motif_analysis.sh IDR_PEAKS GENOME_FASTA CHROM_SIZES MOTIFS
#
# Example:
# scripts/07_motif_analysis/run_motif_analysis.sh \
#   idr_results/Flag-Antp_IDR0.05.narrowPeak \
#   reference/dm6.fa \
#   reference/dm6.chrom.sizes \
#   scripts/07_motif_analysis/motifs.meme

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 IDR_PEAKS GENOME_FASTA CHROM_SIZES MOTIFS" >&2
  exit 1
fi

idr_peaks="$1"
genome_fasta="$2"
chrom_sizes="$3"
motif_file="$4"

output_dir="motif_analysis_results"
summits="${output_dir}/Flag-Antp_summits.bed"
peak_windows="${output_dir}/Flag-Antp_peakcenter200.bed"
peak_fasta="${output_dir}/Flag-Antp_peakcenter200.fa"
fimo_output="${output_dir}/fimo_peakcenter200"

mkdir -p "${output_dir}"

# Confirm that the required input files exist.

for input_file in \
  "${idr_peaks}" \
  "${genome_fasta}" \
  "${chrom_sizes}" \
  "${motif_file}"
do
  if [[ ! -f "${input_file}" ]]; then
    echo "Error: input file not found: ${input_file}" >&2
    exit 1
  fi
done

# 1. Convert the summit offsets in the narrowPeak file into
#    one-base genomic summit intervals.
#
# narrowPeak column 10 contains the summit position relative
# to the start of each peak.

awk 'BEGIN {OFS="\t"}
{
  summit = $2 + $10
  print $1, summit, summit + 1, $4
}' "${idr_peaks}" > "${summits}"

# 2. Generate exactly 200-bp intervals centered on each summit.
#
# The one-base summit interval is extended by 100 bp upstream
# and 99 bp downstream, giving a total interval length of 200 bp.
# bedtools slop also prevents intervals from extending beyond
# chromosome boundaries.

bedtools slop \
  -i "${summits}" \
  -g "${chrom_sizes}" \
  -l 100 \
  -r 99 \
  > "${peak_windows}"

# 3. Extract the corresponding dm6 genomic sequences.
#
# The default FASTA headers retain genomic coordinates in the
# format chromosome:start-end.

bedtools getfasta \
  -fi "${genome_fasta}" \
  -bed "${peak_windows}" \
  -fo "${peak_fasta}"

# 4. Scan the peak-centered sequences for motif occurrences.

fimo \
  --oc "${fimo_output}" \
  --thresh 1e-3 \
  "${motif_file}" \
  "${peak_fasta}" \
  2> "${output_dir}/fimo.log"

echo "Motif occurrence analysis completed."
echo "Peak-centered BED: ${peak_windows}"
echo "Peak-centered FASTA: ${peak_fasta}"
echo "FIMO results: ${fimo_output}"