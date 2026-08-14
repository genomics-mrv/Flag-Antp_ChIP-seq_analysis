#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download Drosophila melanogaster dm6 reference files from UCSC
#
# Files:
#   dm6.fa.gz
#   dm6.chrom.sizes
#
# Usage:
#   bash 05_download_dm6_reference.sh [output_directory]
#
# Default output directory:
#   data/reference
#
# Requirements:
#   - curl
#   - gzip
#
# Notes:
#   - The FASTA is downloaded compressed and then decompressed to dm6.fa.
#   - Bowtie2 index construction is NOT performed here; that belongs in the
#     alignment module.
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

# ------------------------------- SOURCES -------------------------------------

FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
CHROMSIZES_URL="https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes"

FASTA_GZ="$OUTDIR/dm6.fa.gz"
FASTA="$OUTDIR/dm6.fa"
CHROMSIZES="$OUTDIR/dm6.chrom.sizes"

# ------------------------------- DOWNLOAD ------------------------------------

echo "Downloading dm6 reference genome FASTA..."
curl -fL "$FASTA_URL" -o "$FASTA_GZ"

[[ -s "$FASTA_GZ" ]] || {
    echo "ERROR: downloaded FASTA archive is missing or empty: $FASTA_GZ" >&2
    exit 1
}

gzip -t "$FASTA_GZ" || {
    echo "ERROR: gzip integrity check failed: $FASTA_GZ" >&2
    exit 1
}

echo "Decompressing dm6 FASTA..."
gzip -dc "$FASTA_GZ" > "$FASTA"

[[ -s "$FASTA" ]] || {
    echo "ERROR: decompressed FASTA is missing or empty: $FASTA" >&2
    exit 1
}

echo "Downloading dm6 chromosome sizes..."
curl -fL "$CHROMSIZES_URL" -o "$CHROMSIZES"

[[ -s "$CHROMSIZES" ]] || {
    echo "ERROR: chromosome sizes file is missing or empty: $CHROMSIZES" >&2
    exit 1
}

# ------------------------------- VALIDATION ----------------------------------

# Basic FASTA sanity check: file should contain at least one sequence header.
if ! grep -q '^>' "$FASTA"; then
    echo "ERROR: dm6 FASTA does not appear to contain valid sequence headers." >&2
    exit 1
fi

# Basic chromosome-size sanity check:
# first two columns must be chromosome name + positive integer length.
if ! awk 'NF >= 2 && $2 ~ /^[0-9]+$/ && $2 > 0 {valid=1} END{exit !valid}' "$CHROMSIZES"; then
    echo "ERROR: dm6.chrom.sizes does not appear to contain valid chromosome sizes." >&2
    exit 1
fi

# ------------------------------- MANIFEST ------------------------------------

MANIFEST="$OUTDIR/dm6_reference_manifest.tsv"

{
    echo -e "resource\tsource_url\tlocal_file"
    echo -e "dm6_fasta\t${FASTA_URL}\t${FASTA}"
    echo -e "dm6_chrom_sizes\t${CHROMSIZES_URL}\t${CHROMSIZES}"
} > "$MANIFEST"

# ------------------------------- SUMMARY -------------------------------------

echo
echo "dm6 reference files downloaded successfully."
echo
echo "Files:"
echo "  $FASTA"
echo "  $CHROMSIZES"
echo
echo "Compressed FASTA archive retained at:"
echo "  $FASTA_GZ"
echo
echo "Manifest:"
echo "  $MANIFEST"
