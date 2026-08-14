#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Download GFP-Antp embryonic ChIP-seq reads
#
# Samples:
#   ENCFF194GDT  GFP-Antp replicate 1
#   ENCFF716FHU  GFP-Antp replicate 2
#   SRR28101828  GFP-Antp input
#
# The two GFP-Antp ChIP files are single-end 50-bp (SE50) FASTQ files hosted
# by ENCODE and are downloaded directly using their ENCODE file accessions.
#
# The input sample is available as an SRA run and is converted to FASTQ with
# fasterq-dump.
#
# Usage:
#   bash 02_download_GFP_Antp_embryo.sh [output_directory]
#
# Default output directory:
#   data/raw/GFP_Antp_embryo
#
# Requirements:
#   - curl
#   - SRA Toolkit (fasterq-dump)
#   - gzip
# ============================================================================

OUTDIR="${1:-data/raw/GFP_Antp_embryo}"
mkdir -p "$OUTDIR"

# ----------------------------- REQUIREMENTS ----------------------------------

for cmd in curl fasterq-dump gzip; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required program not found in PATH: $cmd" >&2
        exit 1
    }
done

echo "Software:"
curl --version | head -n 1
fasterq-dump --version 2>/dev/null | head -n 1 || true
echo

# --------------------------- ENCODE REPLICATES -------------------------------

REP1_ACCESSION="ENCFF194GDT"
REP2_ACCESSION="ENCFF716FHU"

REP1_URL="https://www.encodeproject.org/files/${REP1_ACCESSION}/@@download/${REP1_ACCESSION}.fastq.gz"
REP2_URL="https://www.encodeproject.org/files/${REP2_ACCESSION}/@@download/${REP2_ACCESSION}.fastq.gz"

echo "Downloading GFP-Antp embryo replicate 1..."
echo "  ENCODE accession: $REP1_ACCESSION"

curl -fL     "$REP1_URL"     -o "$OUTDIR/GFP_Antp_emb_rep1.fastq.gz"

echo "Downloading GFP-Antp embryo replicate 2..."
echo "  ENCODE accession: $REP2_ACCESSION"

curl -fL     "$REP2_URL"     -o "$OUTDIR/GFP_Antp_emb_rep2.fastq.gz"

# Confirm that the downloaded ENCODE files are non-empty.
for file in     "$OUTDIR/GFP_Antp_emb_rep1.fastq.gz"     "$OUTDIR/GFP_Antp_emb_rep2.fastq.gz"
do
    [[ -s "$file" ]] || {
        echo "ERROR: downloaded file is missing or empty: $file" >&2
        exit 1
    }
done

# ------------------------------- INPUT ---------------------------------------

INPUT_SRR="SRR28101828"
TMPDIR="$OUTDIR/.tmp_${INPUT_SRR}"

mkdir -p "$TMPDIR"

echo "Downloading/converting GFP-Antp embryo input..."
echo "  SRA run: $INPUT_SRR"

fasterq-dump     --outdir "$TMPDIR"     --threads 4     "$INPUT_SRR"

# This dataset is expected to be single-end. The paired-end branch is retained
# as a safety check so the script does not silently mishandle unexpected SRA
# metadata.

if [[ -s "$TMPDIR/${INPUT_SRR}.fastq" ]]; then
    mv "$TMPDIR/${INPUT_SRR}.fastq" "$OUTDIR/GFP_Antp_emb_input.fastq"
    gzip -f "$OUTDIR/GFP_Antp_emb_input.fastq"

elif [[ -s "$TMPDIR/${INPUT_SRR}_1.fastq" && -s "$TMPDIR/${INPUT_SRR}_2.fastq" ]]; then
    echo "WARNING: $INPUT_SRR was returned as paired-end data." >&2
    mv "$TMPDIR/${INPUT_SRR}_1.fastq" "$OUTDIR/GFP_Antp_emb_input_R1.fastq"
    mv "$TMPDIR/${INPUT_SRR}_2.fastq" "$OUTDIR/GFP_Antp_emb_input_R2.fastq"
    gzip -f "$OUTDIR/GFP_Antp_emb_input_R1.fastq"
    gzip -f "$OUTDIR/GFP_Antp_emb_input_R2.fastq"

else
    echo "ERROR: expected FASTQ output was not produced for $INPUT_SRR." >&2
    rm -rf "$TMPDIR"
    exit 1
fi

rm -rf "$TMPDIR"

# ------------------------------- SUMMARY -------------------------------------

echo
echo "GFP-Antp embryonic datasets downloaded successfully."
echo
echo "Output directory:"
echo "  $OUTDIR"
echo
echo "Files:"
find "$OUTDIR" -maxdepth 1 -type f -name "*.fastq.gz" -print | sort
