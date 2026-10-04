#!/usr/bin/env bash

# ============================================================
# Flag-Antp ChIP-seq: Quantitative ChIP/input enrichment
#
# Extract mean ChIP/input ratios for:
#   1. Motif-containing peaks (n = 561)
#   2. Motif-free peaks       (n = 627)
#   3. Shared peaks          (n = 739)
#   4. Non-shared peaks      (n = 449)
#
# Requires:
#   - deepTools (multiBigwigSummary)
#   - FlagAntp_fold_enrichment.bw
#   - BED files for the four peak categories
# ============================================================

set -euo pipefail

# 1. Input files
# ============================================================

BIGWIG="FlagAntp_fold_enrichment.bw"

MOTIF_BED="FlagAntp_Motif-Containing_IDR_peaks.bed"
MOTIF_FREE_BED="FlagAntp_Motif-Free_IDR_peaks.bed"

SHARED_BED="Shared_vs_Unique/Flag-Antp_shared.bed"
UNIQUE_BED="Shared_vs_Unique/Flag-Antp_unique.bed"

THREADS=8

# 2. Output directory
# ============================================================

mkdir -p data

# 3. Check input files
# ============================================================

for file in \
  "$BIGWIG" \
  "$MOTIF_BED" \
  "$MOTIF_FREE_BED" \
  "$SHARED_BED" \
  "$UNIQUE_BED"
do
  if [[ ! -f "$file" ]]; then
    echo "ERROR: Missing input file: $file" >&2
    exit 1
  fi
done

# 4. Extract peak-level enrichment
# ============================================================

echo "Extracting motif-containing peak enrichment..."

multiBigwigSummary BED-file \
  --bwfiles "$BIGWIG" \
  --BED "$MOTIF_BED" \
  --outRawCounts data/Motif_Containing_enrichment.tsv \
  -o data/Motif_Containing_enrichment.npz \
  -p "$THREADS"


echo "Extracting motif-free peak enrichment..."

multiBigwigSummary BED-file \
  --bwfiles "$BIGWIG" \
  --BED "$MOTIF_FREE_BED" \
  --outRawCounts data/Motif_Free_enrichment.tsv \
  -o data/Motif_Free_enrichment.npz \
  -p "$THREADS"


echo "Extracting shared peak enrichment..."

multiBigwigSummary BED-file \
  --bwfiles "$BIGWIG" \
  --BED "$SHARED_BED" \
  --outRawCounts data/Shared_enrichment.tsv \
  -o data/Shared_enrichment.npz \
  -p "$THREADS"


echo "Extracting non-shared peak enrichment..."

multiBigwigSummary BED-file \
  --bwfiles "$BIGWIG" \
  --BED "$UNIQUE_BED" \
  --outRawCounts data/Unique_enrichment.tsv \
  -o data/Unique_enrichment.npz \
  -p "$THREADS"


# 5. Verify output
# ============================================================

echo ""
echo "Peak counts in the output files:"
echo "--------------------------------"

for file in \
  data/Motif_Containing_enrichment.tsv \
  data/Motif_Free_enrichment.tsv \
  data/Shared_enrichment.tsv \
  data/Unique_enrichment.tsv
do
  # Subtract the header line
  count=$(( $(wc -l < "$file") - 1 ))

  echo "$(basename "$file"): $count peaks"
done

echo ""
echo "Enrichment extraction completed successfully."