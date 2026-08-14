#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download the Drosophila melanogaster dm6 blacklist
#
# Source:
#   Boyle Lab Blacklist repository
#   dm6-blacklist.v2.bed.gz
#
# Usage:
#   bash 06_download_dm6_blacklist.sh [output_directory]
#
# Default output directory:
#   data/reference
#
# Requirements:
#   - curl
#   - gzip
#
# Output:
#   dm6-blacklist.v2.bed.gz
#   dm6-blacklist.v2.bed
#   dm6_blacklist_manifest.tsv
#
# The compressed source file is retained, and an uncompressed BED copy is
# generated for use in downstream read filtering.
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

BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz"

BLACKLIST_GZ="$OUTDIR/dm6-blacklist.v2.bed.gz"
BLACKLIST_BED="$OUTDIR/dm6-blacklist.v2.bed"

# ------------------------------- DOWNLOAD ------------------------------------

echo "Downloading dm6 blacklist..."
echo "Source:"
echo "  $BLACKLIST_URL"

curl -fL \
    "$BLACKLIST_URL" \
    -o "$BLACKLIST_GZ"

[[ -s "$BLACKLIST_GZ" ]] || {
    echo "ERROR: downloaded blacklist is missing or empty: $BLACKLIST_GZ" >&2
    exit 1
}

# Validate the compressed source file before using it.
gzip -t "$BLACKLIST_GZ" || {
    echo "ERROR: gzip integrity check failed: $BLACKLIST_GZ" >&2
    exit 1
}

# ------------------------------ DECOMPRESS -----------------------------------

echo "Creating uncompressed BED copy..."

gzip -dc "$BLACKLIST_GZ" > "$BLACKLIST_BED"

[[ -s "$BLACKLIST_BED" ]] || {
    echo "ERROR: decompressed blacklist is missing or empty: $BLACKLIST_BED" >&2
    exit 1
}

# ------------------------------- VALIDATION ----------------------------------

# Basic BED sanity check: require at least one line containing chromosome,
# start, and end coordinates with numeric start/end values.
if ! awk '
    $0 !~ /^#/ &&
    NF >= 3 &&
    $2 ~ /^[0-9]+$/ &&
    $3 ~ /^[0-9]+$/ &&
    $3 > $2 {
        valid=1
        exit
    }
    END {exit !valid}
' "$BLACKLIST_BED"; then
    echo "ERROR: blacklist does not appear to contain valid BED intervals." >&2
    exit 1
fi

# ------------------------------- MANIFEST ------------------------------------

MANIFEST="$OUTDIR/dm6_blacklist_manifest.tsv"

{
    echo -e "resource\tassembly\tsource\tlocal_file"
    echo -e "dm6_blacklist_v2\tdm6\t${BLACKLIST_URL}\t${BLACKLIST_GZ}"
    echo -e "dm6_blacklist_v2_uncompressed\tdm6\tderived_from_dm6-blacklist.v2.bed.gz\t${BLACKLIST_BED}"
} > "$MANIFEST"

# ------------------------------- SUMMARY -------------------------------------

echo
echo "dm6 blacklist downloaded successfully."
echo
echo "Files:"
echo "  $BLACKLIST_GZ"
echo "  $BLACKLIST_BED"
echo
echo "Manifest:"
echo "  $MANIFEST"
