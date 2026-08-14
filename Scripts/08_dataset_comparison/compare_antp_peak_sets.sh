#!/usr/bin/env bash
set -euo pipefail

# ==========================================
# Input BED files
# ==========================================
FLAG="Flag-Antp_idr.bed"
GFP_EMB="GFP-Antp_Emb_idr.bed"
GFP_WD="GFP-Antp_WD.bed"

# ==========================================
# Checks
# ==========================================
for f in "$FLAG" "$GFP_EMB" "$GFP_WD"; do
    [[ -f "$f" ]] || { echo "Error: file not found: $f"; exit 1; }
done

command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found"; exit 1; }

# ==========================================
# Output directory
# ==========================================
OUTDIR="Antp_venn_regions"
mkdir -p "$OUTDIR"

# ==========================================
# Sort BED files
# ==========================================
sort -k1,1 -k2,2n "$FLAG"    > "$OUTDIR/Flag-Antp.sorted.bed"
sort -k1,1 -k2,2n "$GFP_EMB" > "$OUTDIR/GFP-Antp_Emb.sorted.bed"
sort -k1,1 -k2,2n "$GFP_WD"  > "$OUTDIR/GFP-Antp_WD.sorted.bed"

A="$OUTDIR/Flag-Antp.sorted.bed"
B="$OUTDIR/GFP-Antp_Emb.sorted.bed"
C="$OUTDIR/GFP-Antp_WD.sorted.bed"

# ==========================================
# 1) Shared by all three
# ==========================================
bedtools intersect -u -a "$A" -b "$B" | \
bedtools intersect -u -a - -b "$C" \
> "$OUTDIR/Flag-Antp__GFP-Antp-Emb__GFP-Antp-WD.shared_all3.bed"

# ==========================================
# 2) Pairwise-only regions
#    (shared by exactly two, excluding third)
# ==========================================

# Flag-Antp ∩ GFP-Antp embryos, but NOT GFP-Antp wing discs
bedtools intersect -u -a "$A" -b "$B" | \
bedtools intersect -v -a - -b "$C" \
> "$OUTDIR/Flag-Antp__GFP-Antp-Emb.only.bed"

# Flag-Antp ∩ GFP-Antp wing discs, but NOT GFP-Antp embryos
bedtools intersect -u -a "$A" -b "$C" | \
bedtools intersect -v -a - -b "$B" \
> "$OUTDIR/Flag-Antp__GFP-Antp-WD.only.bed"

# GFP-Antp embryos ∩ GFP-Antp wing discs, but NOT Flag-Antp
bedtools intersect -u -a "$B" -b "$C" | \
bedtools intersect -v -a - -b "$A" \
> "$OUTDIR/GFP-Antp-Emb__GFP-Antp-WD.only.bed"

# ==========================================
# 3) Unique regions
# ==========================================

# Only Flag-Antp
bedtools intersect -v -a "$A" -b "$B" "$C" \
> "$OUTDIR/Flag-Antp.only.bed"

# Only GFP-Antp embryos
bedtools intersect -v -a "$B" -b "$A" "$C" \
> "$OUTDIR/GFP-Antp-Emb.only.bed"

# Only GFP-Antp wing discs
bedtools intersect -v -a "$C" -b "$A" "$B" \
> "$OUTDIR/GFP-Antp-WD.only.bed"

# ==========================================
# Counts
# ==========================================
FLAG_ONLY=$(wc -l < "$OUTDIR/Flag-Antp.only.bed")
GFP_EMB_ONLY=$(wc -l < "$OUTDIR/GFP-Antp-Emb.only.bed")
GFP_WD_ONLY=$(wc -l < "$OUTDIR/GFP-Antp-WD.only.bed")

FLAG_GFP_EMB_ONLY=$(wc -l < "$OUTDIR/Flag-Antp__GFP-Antp-Emb.only.bed")
FLAG_GFP_WD_ONLY=$(wc -l < "$OUTDIR/Flag-Antp__GFP-Antp-WD.only.bed")
GFP_EMB_GFP_WD_ONLY=$(wc -l < "$OUTDIR/GFP-Antp-Emb__GFP-Antp-WD.only.bed")

ALL3=$(wc -l < "$OUTDIR/Flag-Antp__GFP-Antp-Emb__GFP-Antp-WD.shared_all3.bed")

# Totals
FLAG_TOTAL=$(wc -l < "$A")
GFP_EMB_TOTAL=$(wc -l < "$B")
GFP_WD_TOTAL=$(wc -l < "$C")

# ==========================================
# Summary table
# ==========================================
SUMMARY="$OUTDIR/venn_counts.tsv"

{
    echo -e "region\tcount"
    echo -e "Flag-Antp_only\t$FLAG_ONLY"
    echo -e "GFP-Antp_Emb_only\t$GFP_EMB_ONLY"
    echo -e "GFP-Antp_WD_only\t$GFP_WD_ONLY"
    echo -e "Flag-Antp_and_GFP-Antp_Emb_only\t$FLAG_GFP_EMB_ONLY"
    echo -e "Flag-Antp_and_GFP-Antp_WD_only\t$FLAG_GFP_WD_ONLY"
    echo -e "GFP-Antp_Emb_and_GFP-Antp_WD_only\t$GFP_EMB_GFP_WD_ONLY"
    echo -e "all_three\t$ALL3"
    echo -e "Flag-Antp_total\t$FLAG_TOTAL"
    echo -e "GFP-Antp_Emb_total\t$GFP_EMB_TOTAL"
    echo -e "GFP-Antp_WD_total\t$GFP_WD_TOTAL"
} > "$SUMMARY"

# ==========================================
# Print summary
# ==========================================
echo
echo "Venn regions calculated."
echo "Output directory: $OUTDIR"
echo
cat "$SUMMARY"
echo