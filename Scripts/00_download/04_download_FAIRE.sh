#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download published FAIRE-seq peak BED files from GEO
#
# Datasets:
#   GSM948712  Embryos 2–4 h
#   GSM948713  Embryos 6–8 h
#   GSM948714  Embryos 16–18 h
#   GSM948717  Wing discs
#
# These are processed FAIRE peak BED files provided by GEO.
# IMPORTANT: the downloaded coordinates are in the dm3 genome assembly.
# Conversion to dm6 is performed separately by:
#
#   08_liftover_FAIRE_dm3_to_dm6.sh
#
# Usage:
#   bash 04_download_FAIRE.sh [output_directory]
#
# Default output directory:
#   data/external/FAIRE/dm3
#
# Requirements:
#   - curl
#   - gzip
# ============================================================================

OUTDIR="${1:-data/external/FAIRE/dm3}"
mkdir -p "$OUTDIR"

# ----------------------------- REQUIREMENTS ----------------------------------

for cmd in curl gzip; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required program not found in PATH: $cmd" >&2
        exit 1
    }
done

echo "Software:"
curl --version | head -n 1
gzip --version | head -n 1
echo

# ------------------------------- DATASETS ------------------------------------

declare -A ACCESSIONS=(
    ["FAIRE_2-4h"]="GSM948712"
    ["FAIRE_6-8h"]="GSM948713"
    ["FAIRE_16-18h"]="GSM948714"
    ["FAIRE_wing_disc"]="GSM948717"
)

declare -A FILENAMES=(
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

# Manifest records the exact accession, original GEO filename, genome assembly,
# and local path used for each downloaded file.
MANIFEST="$OUTDIR/FAIRE_download_manifest.tsv"
echo -e "dataset_id\tGEO_accession\tassembly\tfilename\tlocal_path" > "$MANIFEST"

# ------------------------------- DOWNLOAD ------------------------------------

for dataset in "${DATASET_ORDER[@]}"; do
    accession="${ACCESSIONS[$dataset]}"
    filename="${FILENAMES[$dataset]}"

    # GEO supplementary-file endpoint.
    # The filename is URL-encoded only where required by GEO.
    case "$dataset" in
        "FAIRE_2-4h")
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM948712&format=file&file=GSM948712%5Fe2%2D4hr%5FFAIRE%5Fpeaks%2Ebed%2Egz"
            ;;
        "FAIRE_6-8h")
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM948713&format=file&file=GSM948713%5Fe6%2D8hr%5FFAIRE%5Fpeaks%2Ebed%2Egz"
            ;;
        "FAIRE_16-18h")
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM948714&format=file&file=GSM948714%5Fe16%2D18hr%5FFAIRE%5Fpeaks%2Ebed%2Egz"
            ;;
        "FAIRE_wing_disc")
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM948717&format=file&file=GSM948717%5Fwing%5Fdisc%5FFAIRE%5Fpeaks%2Ebed%2Egz"
            ;;
        *)
            echo "ERROR: unknown dataset ID: $dataset" >&2
            exit 1
            ;;
    esac

    outfile="$OUTDIR/$filename"

    echo "============================================================"
    echo "Dataset: $dataset"
    echo "GEO accession: $accession"
    echo "Assembly: dm3"
    echo "Downloading: $filename"

    curl -fL "$url" -o "$outfile"

    [[ -s "$outfile" ]] || {
        echo "ERROR: downloaded file is missing or empty: $outfile" >&2
        exit 1
    }

    # Validate gzip integrity without decompressing the file.
    gzip -t "$outfile" || {
        echo "ERROR: gzip integrity check failed: $outfile" >&2
        exit 1
    }

    echo -e "${dataset}\t${accession}\tdm3\t${filename}\t${outfile}" >> "$MANIFEST"

    echo "Downloaded successfully."
    echo
done

# ------------------------------- SUMMARY -------------------------------------

echo "All FAIRE datasets downloaded successfully."
echo
echo "IMPORTANT:"
echo "  These BED files use dm3 coordinates."
echo "  Convert them to dm6 before downstream comparison/coverage analyses."
echo
echo "Output directory:"
echo "  $OUTDIR"
echo
echo "Manifest:"
echo "  $MANIFEST"
echo
echo "Files:"
find "$OUTDIR" -maxdepth 1 -type f -name "*.bed.gz" -print | sort
