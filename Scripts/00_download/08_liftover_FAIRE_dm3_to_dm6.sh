#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Convert published FAIRE-seq peak coordinates from dm3 to dm6 using UCSC liftOver
#
# Expected inputs:
#   - dm3 FAIRE BED.gz files downloaded by 04_download_FAIRE.sh
#   - dm3ToDm6.over.chain.gz downloaded by 07_download_liftover_chain.sh
#
# Usage:
#   bash 08_liftover_FAIRE_dm3_to_dm6.sh \
#       data/external/FAIRE/dm3 \
#       data/reference/dm3ToDm6.over.chain.gz \
#       data/external/FAIRE/dm6
#
# Default paths:
#   input directory:  data/external/FAIRE/dm3
#   chain file:       data/reference/dm3ToDm6.over.chain.gz
#   output directory: data/external/FAIRE/dm6
#
# Requirements:
#   - UCSC/Kent liftOver
#   - gzip
#   - awk
#   - sort
#
# Outputs:
#   dm6-converted BED files
#   unmapped BED files
#   FAIRE_liftover_summary.tsv
# ============================================================================

INDIR="${1:-data/external/FAIRE/dm3}"
CHAIN="${2:-data/reference/dm3ToDm6.over.chain.gz}"
OUTDIR="${3:-data/external/FAIRE/dm6}"

mkdir -p "$OUTDIR"

# ----------------------------- REQUIREMENTS ----------------------------------

for cmd in liftOver gzip awk sort; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required program not found in PATH: $cmd" >&2
        exit 1
    }
done

[[ -s "$CHAIN" ]] || {
    echo "ERROR: liftOver chain file not found or empty: $CHAIN" >&2
    exit 1
}

gzip -t "$CHAIN" || {
    echo "ERROR: chain file is not a valid gzip archive: $CHAIN" >&2
    exit 1
}

echo "Software:"
liftOver 2>&1 | head -n 1 || true
gzip --version | head -n 1
echo

# ------------------------------- INPUT FILES ---------------------------------

declare -A FILES=(
    ["FAIRE_2-4h"]="GSM948712_e2-4hr_FAIRE_peaks.bed.gz"
    ["FAIRE_6-8h"]="GSM948713_e6-8hr_FAIRE_peaks.bed.gz"
    ["FAIRE_16-18h"]="GSM948714_e16-18hr_FAIRE_peaks.bed.gz"
    ["FAIRE_wing_disc"]="GSM948717_wing_disc_FAIRE_peaks.bed.gz"
)

DATASET_ORDER=(
    "FAIRE_2-4h"
    "FAIRE_6-8h"
    "FAIRE_16-18h"
    "FAIRE_wing_disc"
)

SUMMARY="$OUTDIR/FAIRE_liftover_summary.tsv"
echo -e "dataset_id\tinput_dm3_intervals\tmapped_dm6_intervals\tunmapped_intervals\tpercent_mapped\toutput_bed" > "$SUMMARY"

# ------------------------------- CONVERSION ----------------------------------

for dataset in "${DATASET_ORDER[@]}"; do
    infile_gz="$INDIR/${FILES[$dataset]}"

    [[ -s "$infile_gz" ]] || {
        echo "ERROR: required FAIRE file not found or empty: $infile_gz" >&2
        exit 1
    }

    gzip -t "$infile_gz" || {
        echo "ERROR: invalid gzip file: $infile_gz" >&2
        exit 1
    }

    tmp_dm3="$OUTDIR/.${dataset}.dm3.tmp.bed"
    mapped="$OUTDIR/${dataset}.dm6.bed"
    unmapped="$OUTDIR/${dataset}.unmapped.bed"

    echo "============================================================"
    echo "Dataset: $dataset"
    echo "Input:   $infile_gz"

    # Decompress and normalize to the first three BED columns.
    # liftOver only requires chromosome, start, and end for these peak intervals.
    gzip -dc "$infile_gz" \
        | awk 'BEGIN{OFS="\t"}
               $0 !~ /^#/ &&
               NF >= 3 &&
               $2 ~ /^[0-9]+$/ &&
               $3 ~ /^[0-9]+$/ &&
               $3 > $2 {
                   print $1,$2,$3
               }' \
        | LC_ALL=C sort -k1,1 -k2,2n \
        > "$tmp_dm3"

    input_count=$(wc -l < "$tmp_dm3" | tr -d ' ')

    if [[ "$input_count" -eq 0 ]]; then
        echo "ERROR: no valid BED intervals found in $infile_gz" >&2
        rm -f "$tmp_dm3"
        exit 1
    fi

    echo "Running liftOver dm3 -> dm6..."

    liftOver \
        "$tmp_dm3" \
        "$CHAIN" \
        "$mapped" \
        "$unmapped"

    # Sort mapped intervals for consistent downstream use.
    if [[ -s "$mapped" ]]; then
        LC_ALL=C sort -k1,1 -k2,2n "$mapped" -o "$mapped"
    fi

    mapped_count=$(wc -l < "$mapped" | tr -d ' ')

    # UCSC liftOver's unmapped file contains explanatory comment lines beginning
    # with '#'. Count only BED-like data rows as unmapped intervals.
    unmapped_count=$(awk '
        $0 !~ /^#/ &&
        NF >= 3 &&
        $2 ~ /^[0-9]+$/ &&
        $3 ~ /^[0-9]+$/ {
            n++
        }
        END {print n+0}
    ' "$unmapped")

    percent_mapped=$(awk -v m="$mapped_count" -v n="$input_count" \
        'BEGIN {printf "%.4f", (n > 0 ? 100*m/n : 0)}')

    echo -e "${dataset}\t${input_count}\t${mapped_count}\t${unmapped_count}\t${percent_mapped}\t${mapped}" \
        >> "$SUMMARY"

    rm -f "$tmp_dm3"

    echo "Input intervals:    $input_count"
    echo "Mapped intervals:   $mapped_count"
    echo "Unmapped intervals: $unmapped_count"
    echo "Mapped:              ${percent_mapped}%"
    echo
done

# ------------------------------- SUMMARY -------------------------------------

echo "FAIRE dm3 -> dm6 liftOver completed successfully."
echo
echo "Output directory:"
echo "  $OUTDIR"
echo
echo "Summary:"
echo "  $SUMMARY"
echo
echo "Converted BED files:"
find "$OUTDIR" -maxdepth 1 -type f -name "*.dm6.bed" -print | sort
