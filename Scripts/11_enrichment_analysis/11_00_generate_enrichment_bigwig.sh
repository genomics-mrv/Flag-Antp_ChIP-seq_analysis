#!/usr/bin/env bash

# ============================================================
# Flag-Antp ChIP-seq: Generate ChIP/input enrichment bigWig
#
# Steps:
#   1. Merge the two Flag-Antp ChIP-seq replicates
#   2. Index the merged BAM file
#   3. Calculate the normalized ChIP/input ratio
#
# Requirements:
#   - SAMtools
#   - deepTools (bamCompare)
#
# Output:
#   FlagAntp_Embryos_S5-11_MERGED.bam
#   FlagAntp_Embryos_S5-11_MERGED.bam.bai
#   FlagAntp_fold_enrichment.bw
# ============================================================

set -euo pipefail

# 1. Working directory and parameters
# ============================================================

WORKDIR="${1:-.}"
THREADS=8

REP1="${WORKDIR}/FlagAntp_Embryos_S5-11.5_Rep1.bam"
REP2="${WORKDIR}/FlagAntp_Embryos_S5-11.5_Rep2.bam"
INPUT="${WORKDIR}/FlagAntp_Embryos_S5-11.5_Input.bam"

MERGED="${WORKDIR}/FlagAntp_Embryos_S5-11_MERGED.bam"
BIGWIG="${WORKDIR}/FlagAntp_fold_enrichment.bw"

# 2. Check dependencies and input files
# ============================================================

for program in samtools bamCompare; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "ERROR: $program is not installed or is not in PATH." >&2
    exit 1
  fi
done

for file in "$REP1" "$REP2" "$INPUT"; do
  if [[ ! -f "$file" ]]; then
    echo "ERROR: Missing input file: $file" >&2
    exit 1
  fi
done

# Prevent accidental overwriting of existing results
for file in "$MERGED" "$BIGWIG"; do
  if [[ -e "$file" ]]; then
    echo "ERROR: Output already exists: $file" >&2
    echo "Move or remove it before rerunning the script." >&2
    exit 1
  fi
done

# 3. Merge ChIP-seq biological replicates
# ============================================================

echo "Merging Flag-Antp ChIP-seq replicates..."

samtools merge \
  -@ "$THREADS" \
  "$MERGED" \
  "$REP1" \
  "$REP2"

# 4. Index the merged BAM
# ============================================================

echo "Indexing merged BAM..."

samtools index \
  -@ "$THREADS" \
  "$MERGED"

# 5. Generate normalized ChIP/input ratio bigWig
# ============================================================

echo "Calculating normalized ChIP/input ratios..."

bamCompare \
  -b1 "$MERGED" \
  -b2 "$INPUT" \
  --operation ratio \
  --scaleFactorsMethod None \
  --normalizeUsing CPM \
  --binSize 10 \
  --pseudocount 1 \
  -p "$THREADS" \
  -o "$BIGWIG"

# 6. Verify outputs
# ============================================================

echo ""
echo "Generated files:"
echo "--------------------------------"

for file in \
  "$MERGED" \
  "${MERGED}.bai" \
  "$BIGWIG"
do
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty output: $file" >&2
    exit 1
  fi

  echo "$file"
done

echo ""
echo "ChIP/input enrichment calculation completed successfully."
