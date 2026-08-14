#!/usr/bin/env bash

set -euo pipefail

# ==================================================================================================
# Calculate genomic representation of promoters, exons, introns, and intergenic regions
# ==================================================================================================
#
# Purpose:
#   This script calculates the percentage of the genome occupied by four mutually exclusive
#   annotation categories: promoters, exons, introns, and intergenic regions.
#
# Why this is needed:
#   Promoters, exons, and introns can overlap. For example, a promoter window around a TSS may
#   overlap an exon or intron. If raw BED files are counted independently, overlapping bases are
#   counted more than once and the percentages become inflated. To avoid this, the genome is
#   partitioned using the same priority hierarchy used for sequential ChIP-seq peak annotation:
#
#       promoter > exon > intron > intergenic
#
# Final categories:
#   1. Promoter   = merged promoter intervals
#   2. Exon       = merged exon intervals after removing promoter bases
#   3. Intron     = merged intron intervals after removing promoter and exon bases
#   4. Intergenic = all genomic bases not assigned to promoter, exon, or intron
#
# Requirements:
#   - bedtools
#   - awk
#   - sort
#
# Input files:
#   - Promoters.bed       BED file with promoter intervals, e.g. ±1 kb around TSS
#   - Exons.bed           BED file with exon intervals
#   - Introns.bed         BED file with intron intervals
#   - dm6.chrom.sizes     Chromosome sizes file: chr<TAB>length
#
# Usage:
#   bash calculate_genomic_feature_percentages.sh \
#        Promoters.bed Exons.bed Introns.bed dm6.chrom.sizes
#
# Output directory:
#   genomic_feature_percentages/
#
# Main output:
#   genomic_feature_percentages/genomic_feature_percentages.tsv
#
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# 1. Read input arguments
# --------------------------------------------------------------------------------------------------
# The script expects exactly four input files. The chromosome sizes file is essential because
# intergenic DNA is calculated as the complement of promoter/exon/intron regions across the genome.

PROMOTERS=${1:-Promoters.bed}
EXONS=${2:-Exons.bed}
INTRONS=${3:-Introns.bed}
CHROMSIZES=${4:-dm6.chrom.sizes}

OUTDIR="genomic_feature_percentages"
mkdir -p "$OUTDIR"

# --------------------------------------------------------------------------------------------------
# 2. Check that all required files and programs are available
# --------------------------------------------------------------------------------------------------

for f in "$PROMOTERS" "$EXONS" "$INTRONS" "$CHROMSIZES"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Required input file not found or empty: $f" >&2
        exit 1
    fi
done

if ! command -v bedtools >/dev/null 2>&1; then
    echo "ERROR: bedtools is not installed or not available in PATH." >&2
    exit 1
fi

# --------------------------------------------------------------------------------------------------
# 3. Sort and merge each annotation class independently
# --------------------------------------------------------------------------------------------------
# Merging removes redundant or overlapping intervals within the same category.
# This is critical for a base-pair analysis: each genomic base must be counted
# only once within a feature before the promoter > exon > intron hierarchy is applied.

bedtools sort -i "$PROMOTERS" | bedtools merge -i - > "$OUTDIR/promoters.merged.bed"
bedtools sort -i "$EXONS"     | bedtools merge -i - > "$OUTDIR/exons.merged.bed"
bedtools sort -i "$INTRONS"   | bedtools merge -i - > "$OUTDIR/introns.merged.bed"

# --------------------------------------------------------------------------------------------------
# 4. Keep promoters as the highest-priority category
# --------------------------------------------------------------------------------------------------
# Promoter intervals are not subtracted from anything at this step because promoters have
# highest priority. Any genomic base that belongs to both a promoter and another annotation
# is assigned to promoter, enforcing the hierarchy used throughout the manuscript.

cp "$OUTDIR/promoters.merged.bed" "$OUTDIR/promoters.final.bed"

# --------------------------------------------------------------------------------------------------
# 5. Define final exonic space after removing promoter-overlapping bases
# --------------------------------------------------------------------------------------------------
# Exons are the second-priority class. `bedtools subtract` removes only the exon portions
# overlapping promoters, so the remaining exon intervals are mutually exclusive from promoters.

bedtools subtract \
    -a "$OUTDIR/exons.merged.bed" \
    -b "$OUTDIR/promoters.final.bed" \
    > "$OUTDIR/exons.final.bed"

# --------------------------------------------------------------------------------------------------
# 6. Define final intronic space after removing promoter and exon bases
# --------------------------------------------------------------------------------------------------
# Introns are the third-priority class. Any intron base overlapping a promoter or exon is removed.
# This is important for cases such as promoter windows overlapping introns, nested genes, or exons
# from one transcript/gene overlapping introns from another.

cat "$OUTDIR/promoters.final.bed" "$OUTDIR/exons.final.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$OUTDIR/promoter_exon.final.merged.bed"

bedtools subtract \
    -a "$OUTDIR/introns.merged.bed" \
    -b "$OUTDIR/promoter_exon.final.merged.bed" \
    > "$OUTDIR/introns.final.bed"

# --------------------------------------------------------------------------------------------------
# 7. Define intergenic space as the complement of all assigned regions
# --------------------------------------------------------------------------------------------------
# Intergenic DNA is defined mathematically as the complement of all assigned promoter/exon/intron
# bases across the chromosome sizes file. This guarantees that the four final categories partition
# the genome and should sum to ~100%.

cat "$OUTDIR/promoters.final.bed" "$OUTDIR/exons.final.bed" "$OUTDIR/introns.final.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$OUTDIR/annotated_non_intergenic.final.bed"

bedtools complement \
    -i "$OUTDIR/annotated_non_intergenic.final.bed" \
    -g "$CHROMSIZES" \
    > "$OUTDIR/intergenic.final.bed"

# --------------------------------------------------------------------------------------------------
# 8. Calculate base-pair coverage for each final category
# --------------------------------------------------------------------------------------------------
# BED coordinates are 0-based and half-open. Therefore, the number of covered bases in an
# interval is exactly end - start; adding 1 here would systematically inflate all coverage values.

promoter_bp=$(awk '{sum += $3 - $2} END {print sum + 0}' "$OUTDIR/promoters.final.bed")
exon_bp=$(awk     '{sum += $3 - $2} END {print sum + 0}' "$OUTDIR/exons.final.bed")
intron_bp=$(awk   '{sum += $3 - $2} END {print sum + 0}' "$OUTDIR/introns.final.bed")
intergenic_bp=$(awk '{sum += $3 - $2} END {print sum + 0}' "$OUTDIR/intergenic.final.bed")

genome_bp=$(awk '{sum += $2} END {print sum + 0}' "$CHROMSIZES")
assigned_bp=$((promoter_bp + exon_bp + intron_bp + intergenic_bp))

# --------------------------------------------------------------------------------------------------
# 9. Write output table with bp counts and genome percentages
# --------------------------------------------------------------------------------------------------
# The percentages are calculated relative to the total genome size in the chromosome sizes file.
# The final categories should sum to approximately 100%, allowing for rounding.

{
    echo -e "category\tbp\tgenome_percent"
    awk -v bp="$promoter_bp"   -v genome="$genome_bp" 'BEGIN {printf "promoter\t%d\t%.6f\n", bp, 100*bp/genome}'
    awk -v bp="$exon_bp"       -v genome="$genome_bp" 'BEGIN {printf "exon\t%d\t%.6f\n", bp, 100*bp/genome}'
    awk -v bp="$intron_bp"     -v genome="$genome_bp" 'BEGIN {printf "intron\t%d\t%.6f\n", bp, 100*bp/genome}'
    awk -v bp="$intergenic_bp" -v genome="$genome_bp" 'BEGIN {printf "intergenic\t%d\t%.6f\n", bp, 100*bp/genome}'
    awk -v bp="$assigned_bp"   -v genome="$genome_bp" 'BEGIN {printf "total\t%d\t%.6f\n", bp, 100*bp/genome}'
} > "$OUTDIR/genomic_feature_percentages.tsv"

# --------------------------------------------------------------------------------------------------
# 10. Print a concise summary to the terminal
# --------------------------------------------------------------------------------------------------

echo "Done."
echo "Genome size used: $genome_bp bp"
echo "Assigned bp total: $assigned_bp bp"
echo "Output table: $OUTDIR/genomic_feature_percentages.tsv"
echo
echo "Results:"
column -t "$OUTDIR/genomic_feature_percentages.tsv" 2>/dev/null || cat "$OUTDIR/genomic_feature_percentages.tsv"

