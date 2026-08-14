#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download GFP-Antp wing-disc ChIP-seq reads from GEO/SRA
#
# Samples:
#   GSM3578082  Antp-GFP ChIP, wing disc
#   GSM3578083  Antp-GFP input, wing disc
#
# The script resolves each GEO sample accession (GSM) to its underlying SRA
# run accession(s) using pysradb, then converts each run to FASTQ with
# fasterq-dump and compresses the result with gzip.
#
# Usage:
#   bash 03_download_GFP_Antp_wing_disc.sh [output_directory]
#
# Default output directory:
#   data/raw/GFP_Antp_wing_disc
#
# Requirements:
#   - pysradb
#   - SRA Toolkit (fasterq-dump)
#   - gzip
#   - awk, sort
# ============================================================================

OUTDIR="${1:-data/raw/GFP_Antp_wing_disc}"
mkdir -p "$OUTDIR"

# ----------------------------- REQUIREMENTS ----------------------------------

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

# ------------------------------- SAMPLES -------------------------------------

declare -A SAMPLE_NAMES=(
    ["GSM3578082"]="GFP_Antp_wing_ChIP"
    ["GSM3578083"]="GFP_Antp_wing_input"
)

GSM_ACCESSIONS=(
    "GSM3578082"
    "GSM3578083"
)

# Save the exact GEO-to-SRA mapping used during the download.
MAPPING_FILE="$OUTDIR/GEO_to_SRA_runs.tsv"
echo -e "GEO_accession\tsample_name\tSRA_run" > "$MAPPING_FILE"

# ------------------------------- DOWNLOAD ------------------------------------

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

        # Preserve individual runs if one GEO sample contains multiple SRR files.
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

        # Handle either single-end or paired-end output without assuming the
        # library layout from the GEO accession alone.
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

# ------------------------------- SUMMARY -------------------------------------

echo "GFP-Antp wing-disc datasets downloaded successfully."
echo
echo "Output directory:"
echo "  $OUTDIR"
echo
echo "GEO-to-SRA mapping:"
echo "  $MAPPING_FILE"
echo
echo "FASTQ files:"
find "$OUTDIR" -maxdepth 1 -type f -name "*.fastq.gz" -print | sort
