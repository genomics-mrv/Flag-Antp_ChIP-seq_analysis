#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download the UCSC liftOver chain file for dm3 -> dm6 conversion
#
# File:
#   dm3ToDm6.over.chain.gz
#
# This chain file is used to convert the published FAIRE-seq peak coordinates
# from the Drosophila melanogaster dm3 genome assembly to dm6.
#
# Usage:
#   bash 07_download_liftover_chain.sh [output_directory]
#
# Default output directory:
#   data/reference
#
# Requirements:
#   - curl
#   - gzip
#
# Output:
#   dm3ToDm6.over.chain.gz
#   liftover_chain_manifest.tsv
#
# The actual coordinate conversion is performed separately by:
#   08_liftover_FAIRE_dm3_to_dm6.sh
# ============================================================================

OUTDIR="${1:-data/reference}"
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

# ------------------------------- SOURCE --------------------------------------

CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz"
CHAIN_FILE="$OUTDIR/dm3ToDm6.over.chain.gz"

# ------------------------------- DOWNLOAD ------------------------------------

echo "Downloading UCSC dm3 -> dm6 liftOver chain..."
echo "Source:"
echo "  $CHAIN_URL"

curl -fL \
    "$CHAIN_URL" \
    -o "$CHAIN_FILE"

[[ -s "$CHAIN_FILE" ]] || {
    echo "ERROR: downloaded chain file is missing or empty: $CHAIN_FILE" >&2
    exit 1
}

# Validate the compressed chain file.
gzip -t "$CHAIN_FILE" || {
    echo "ERROR: gzip integrity check failed: $CHAIN_FILE" >&2
    exit 1
}

# ------------------------------- MANIFEST ------------------------------------

MANIFEST="$OUTDIR/liftover_chain_manifest.tsv"

{
    echo -e "resource\tsource_assembly\ttarget_assembly\tsource_url\tlocal_file"
    echo -e "dm3ToDm6_chain\tdm3\tdm6\t${CHAIN_URL}\t${CHAIN_FILE}"
} > "$MANIFEST"

# ------------------------------- SUMMARY -------------------------------------

echo
echo "liftOver chain downloaded successfully."
echo
echo "File:"
echo "  $CHAIN_FILE"
echo
echo "Manifest:"
echo "  $MANIFEST"
echo
echo "Use this chain with:"
echo "  08_liftover_FAIRE_dm3_to_dm6.sh"
