#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download Flag-Antp ChIP-seq reads from GEO/SRA
#
# Samples:
#   GSM9491492  Flag-Antp replicate 1
#   GSM9491493  Flag-Antp replicate 2
#   GSM9491494  Flag-Antp input
#
# The script resolves each GEO sample accession (GSM) to its underlying SRA
# run accession(s) using pysradb, then converts each run to FASTQ with
# fasterq-dump and compresses the result with gzip.
#
# Usage:
#   bash 01_download_FlagAntp_GEO.sh [output_directory]
#
# Default output directory:
#   data/raw/FlagAntp
#
# Requirements:
#   - pysradb
#   - SRA Toolkit (fasterq-dump)
#   - gzip
#   - awk, sort
# ============================================================================

OUTDIR="${1:-data/raw/FlagAntp}"
mkdir -p "$OUTDIR"

for cmd in pysradb fasterq-dump gzip awk sort; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required program not found in PATH: $cmd" >&2
        exit 1
    }
done

echo "Software:"
pysradb --version 2>/dev/null || true
fasterq-dump --version 2>/dev/null | head -n 1 || true
echo

declare -A SAMPLE_NAMES=(
    ["GSM9491492"]="FlagAntp_rep1"
    ["GSM9491493"]="FlagAntp_rep2"
    ["GSM9491494"]="FlagAntp_input"
)

GSM_ACCESSIONS=(
    "GSM9491492"
    "GSM9491493"
    "GSM9491494"
)

MAPPING_FILE="$OUTDIR/GEO_to_SRA_runs.tsv"
echo -e "GEO_accession\tsample_name\tSRA_run" > "$MAPPING_FILE"

for gsm in "${GSM_ACCESSIONS[@]}"; do
    sample="${SAMPLE_NAMES[$gsm]}"

    echo "============================================================"
    echo "Sample: $sample"
    echo "GEO accession: $gsm"
    echo "Resolving SRA run accession(s)..."

    mapfile -t runs < <(
        pysradb gsm-to-srr "$gsm" \
        | awk 'NR > 1 {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^SRR[0-9]+$/) print $i
            }
        }' \
        | sort -u
    )

    if [[ ${#runs[@]} -eq 0 ]]; then
        echo "ERROR: no SRR accession could be resolved for $gsm." >&2
        exit 1
    fi

    echo "Resolved run(s): ${runs[*]}"

    run_number=0
    for srr in "${runs[@]}"; do
        run_number=$((run_number + 1))
        echo -e "${gsm}\t${sample}\t${srr}" >> "$MAPPING_FILE"

        if [[ ${#runs[@]} -eq 1 ]]; then
            prefix="$sample"
        else
            prefix="${sample}_run${run_number}_${srr}"
        fi

        tmpdir="$OUTDIR/.tmp_${srr}"
        mkdir -p "$tmpdir"

        echo "Downloading/converting $srr with fasterq-dump..."

        fasterq-dump \
            --outdir "$tmpdir" \
            --threads 4 \
            "$srr"

        if [[ -s "$tmpdir/${srr}.fastq" ]]; then
            mv "$tmpdir/${srr}.fastq" "$OUTDIR/${prefix}.fastq"
            gzip -f "$OUTDIR/${prefix}.fastq"

        elif [[ -s "$tmpdir/${srr}_1.fastq" && -s "$tmpdir/${srr}_2.fastq" ]]; then
            mv "$tmpdir/${srr}_1.fastq" "$OUTDIR/${prefix}_R1.fastq"
            mv "$tmpdir/${srr}_2.fastq" "$OUTDIR/${prefix}_R2.fastq"
            gzip -f "$OUTDIR/${prefix}_R1.fastq"
            gzip -f "$OUTDIR/${prefix}_R2.fastq"

        else
            echo "ERROR: expected FASTQ output was not produced for $srr." >&2
            rm -rf "$tmpdir"
            exit 1
        fi

        rm -rf "$tmpdir"
        echo "Finished: $srr"
        echo
    done
done

echo "All Flag-Antp GEO samples downloaded successfully."
echo
echo "Output directory:"
echo "  $OUTDIR"
echo
echo "GEO-to-SRA mapping:"
echo "  $MAPPING_FILE"
echo
echo "FASTQ files:"
find "$OUTDIR" -maxdepth 1 -type f -name "*.fastq.gz" -print | sort
