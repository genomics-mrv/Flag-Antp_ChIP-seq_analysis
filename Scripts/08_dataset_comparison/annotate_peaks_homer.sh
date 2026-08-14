#!/usr/bin/env bash
set -euo pipefail

# Input
PEAKS="GFP-Antp_Emb_idr.bed"
GENOME="dm6"
OUT="GFP-Antp_Emb_idr_annotated.txt"

# Check HOMER
command -v annotatePeaks.pl >/dev/null 2>&1 || {
    echo "Error: annotatePeaks.pl not found in PATH"
    exit 1
}

# Check input
[[ -f "$PEAKS" ]] || {
    echo "Error: file not found: $PEAKS"
    exit 1
}

# Annotate
annotatePeaks.pl "$PEAKS" "$GENOME" > "$OUT"

echo "Done."
echo "Output: $OUT"