#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Base-pair-weighted genomic annotation
# Priority: promoter > exon > intron > intergenic
#
# Outputs:
#   1) mutually exclusive annotation BED files
#   2) ChIP and FAIRE base-pair coverage by genomic feature
#   3) per-peak fractional ChIP annotation
#   4) combined genome / FAIRE / ChIP percentages
#   5) ChIP enrichment relative to genome and accessible DNA
# ============================================================================

# ------------------------------- INPUTS --------------------------------------
# All input paths are supplied on the command line rather than being hard-coded
# to a particular computer. This makes the analysis portable and reproducible.
#
# Usage:
#   bash calculate_feature_coverage.sh \
#       Promoters.bed \
#       Exons.bed \
#       Introns.bed \
#       ChIP_peaks.bed \
#       FAIRE_peaks.bed \
#       genomic_feature_percentages.tsv \
#       output_directory
#
# The ChIP and FAIRE BED files may contain more than three columns; only genomic
# coordinates (chromosome, start, end) are used for the coverage calculation.

if [[ $# -ne 7 ]]; then
    echo "Usage: $0 Promoters.bed Exons.bed Introns.bed ChIP_peaks.bed FAIRE_peaks.bed genomic_feature_percentages.tsv output_directory" >&2
    exit 1
fi

PROMOTER=$1
EXON=$2
INTRON=$3
CHIP=$4
FAIRE=$5
GENOME_PERCENTAGES=$6
OUTDIR=$7

# GENOME_PERCENTAGES must have been generated from the SAME promoter/exon/intron
# annotation files used here. Mixing annotation versions would make enrichment
# values internally inconsistent.
# Expected columns: category, bp, genome_percent

mkdir -p "$OUTDIR"/{annotation,chip,faire,tmp}

# ----------------------------- REQUIREMENTS ----------------------------------
command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools was not found in PATH." >&2
    exit 1
}

for file in "$EXON" "$INTRON" "$PROMOTER" "$CHIP" "$FAIRE" "$GENOME_PERCENTAGES"; do
    [[ -s "$file" ]] || {
        echo "ERROR: missing or empty input file: $file" >&2
        exit 1
    }
done

echo "Base-pair-weighted annotation"
echo "Feature priority: promoter > exon > intron > intergenic"
echo "Output directory: $OUTDIR"

# -------------------------- HELPER FUNCTIONS ---------------------------------
# Keep valid BED records, retain only genomic coordinates, sort them, and merge
# overlapping intervals. Merging is essential because coverage is measured in
# unique base pairs; otherwise overlapping peaks/regions would be counted twice.
normalize_bed3() {
    local input=$1
    local output=$2

    awk 'BEGIN{OFS="\t"}
         $0 !~ /^#/ && NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {
             print $1,$2,$3
         }' "$input" \
      | LC_ALL=C sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "$output"
}

# Sum the length of BED3 intervals.
bed_bp() {
    awk '{sum += $3-$2} END{print sum+0}' "$1"
}

# Summarize unique covered bp from a dataset against the mutually exclusive
# annotation. `bedtools intersect -wo` reports the exact number of overlapping
# base pairs, which is the key quantity used for the base-pair-weighted analysis.
summarize_coverage() {
    local dataset_name=$1
    local input=$2
    local workdir=$3

    local merged="$workdir/${dataset_name}.merged.bed"
    local overlap="$workdir/${dataset_name}.feature_overlap.tsv"
    local summary="$workdir/${dataset_name}_coverage_by_feature.tsv"

    normalize_bed3 "$input" "$merged"

    bedtools intersect \
        -a "$merged" \
        -b "$OUTDIR/annotation/features.exclusive.bed" \
        -wo \
        > "$overlap"

    local total_bp
    total_bp=$(bed_bp "$merged")

    awk -v total="$total_bp" 'BEGIN{OFS="\t"}
        {bp[$7] += $8}
        END {
            order[1]="promoter"; order[2]="exon"; order[3]="intron";
            assigned=0;
            for (i=1; i<=3; i++) assigned += bp[order[i]]+0;
            bp["intergenic"] = total-assigned;

            print "dataset","category","bp","percent";
            for (i=1; i<=4; i++) {
                if (i==4) category="intergenic"; else category=order[i];
                pct=(total>0 ? 100*bp[category]/total : 0);
                printf "%s\t%s\t%.0f\t%.6f\n", DATASET,category,bp[category]+0,pct;
            }
            printf "%s\ttotal\t%.0f\t100.000000\n", DATASET,total;
        }' DATASET="$dataset_name" "$overlap" > "$summary"

    echo "  $dataset_name: $total_bp unique bp"
}

# ------------------ BUILD MUTUALLY EXCLUSIVE FEATURES ------------------------
echo "Building mutually exclusive genomic features..."

normalize_bed3 "$PROMOTER" "$OUTDIR/annotation/promoter.merged.bed"
normalize_bed3 "$EXON"     "$OUTDIR/annotation/exon.raw.merged.bed"
normalize_bed3 "$INTRON"   "$OUTDIR/annotation/intron.raw.merged.bed"

# Promoter has first priority. Bases falling in both a promoter and another
# feature are therefore counted as promoter, matching the hierarchy described
# in the manuscript: promoter > exon > intron > intergenic.
cp "$OUTDIR/annotation/promoter.merged.bed" \
   "$OUTDIR/annotation/promoter.exclusive.bed"

# Remove promoter-overlapping bases from exons. This operation makes promoter
# and exon intervals mutually exclusive and prevents double-counting.
bedtools subtract \
    -a "$OUTDIR/annotation/exon.raw.merged.bed" \
    -b "$OUTDIR/annotation/promoter.exclusive.bed" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$OUTDIR/annotation/exon.exclusive.bed"

# Mask promoter + exclusive exon before constructing introns. Any intronic base
# overlapping either higher-priority category is removed before coverage is
# calculated.
cat "$OUTDIR/annotation/promoter.exclusive.bed" \
    "$OUTDIR/annotation/exon.exclusive.bed" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$OUTDIR/annotation/promoter_exon.mask.bed"

bedtools subtract \
    -a "$OUTDIR/annotation/intron.raw.merged.bed" \
    -b "$OUTDIR/annotation/promoter_exon.mask.bed" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$OUTDIR/annotation/intron.exclusive.bed"

# Add feature labels. Intergenic coverage is calculated later as the residual:
# total unique dataset bp minus bp assigned to promoter, exon, and intron.
# This lets peaks contribute proportionally to multiple categories while ensuring
# that each covered base is counted exactly once.
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"promoter"}' \
    "$OUTDIR/annotation/promoter.exclusive.bed" \
    > "$OUTDIR/annotation/features.exclusive.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"exon"}' \
    "$OUTDIR/annotation/exon.exclusive.bed" \
    >> "$OUTDIR/annotation/features.exclusive.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"intron"}' \
    "$OUTDIR/annotation/intron.exclusive.bed" \
    >> "$OUTDIR/annotation/features.exclusive.bed"

LC_ALL=C sort -k1,1 -k2,2n "$OUTDIR/annotation/features.exclusive.bed" \
    -o "$OUTDIR/annotation/features.exclusive.bed"

# Validate that the final promoter/exon/intron annotation contains no overlaps.
# If this check fails, the base-pair percentages would be inflated by duplicated
# genomic bases and the script stops rather than producing misleading results.
bedtools intersect \
    -a "$OUTDIR/annotation/features.exclusive.bed" \
    -b "$OUTDIR/annotation/features.exclusive.bed" \
    -wo \
  | awk '$4 != $8 || $2 != $6 || $3 != $7' \
  > "$OUTDIR/tmp/annotation_self_overlap_check.tsv"

if [[ -s "$OUTDIR/tmp/annotation_self_overlap_check.tsv" ]]; then
    echo "ERROR: exclusive annotation contains unexpected overlaps." >&2
    exit 1
fi

# --------------------- CHiP AND FAIRE COVERAGE -------------------------------
echo "Calculating ChIP and FAIRE base-pair distributions..."
summarize_coverage "ChIP"  "$CHIP"  "$OUTDIR/chip"
summarize_coverage "FAIRE" "$FAIRE" "$OUTDIR/faire"

# ------------------- PER-PEAK FRACTIONAL ANNOTATION --------------------------
echo "Calculating fractional annotation for each ChIP peak..."

# Assign a unique internal ID to every ChIP interval. narrowPeak/BED names are
# not guaranteed to be unique, so using generated IDs prevents different peaks
# from being accidentally combined during per-peak fractional annotation.
awk 'BEGIN{OFS="\t"}
     $0 !~ /^#/ && NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {
         original=(NF>=4 ? $4 : ".");
         print $1,$2,$3,"peak_"++n,original
     }' "$CHIP" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > "$OUTDIR/chip/ChIP.peaks.with_ids.bed"

bedtools intersect \
    -a "$OUTDIR/chip/ChIP.peaks.with_ids.bed" \
    -b "$OUTDIR/annotation/features.exclusive.bed" \
    -wo \
  > "$OUTDIR/chip/ChIP.per_peak.overlaps.tsv"

# Aggregate the exact overlap length for each peak × feature combination.
# A single peak can therefore contribute, for example, partly to promoter and
# partly to exon rather than being forced into one categorical annotation.
awk 'BEGIN{OFS="\t"}
     {key=$4 SUBSEP $9; bp[key]+=$10}
     END{
         for (key in bp) {
             split(key,a,SUBSEP);
             print a[1],a[2],bp[key]
         }
     }' "$OUTDIR/chip/ChIP.per_peak.overlaps.tsv" \
  > "$OUTDIR/tmp/per_peak_category_bp.tsv"

# Emit promoter/exon/intron/intergenic fractions for every peak, including zeros.
awk 'BEGIN{OFS="\t"}
     FNR==NR {
         chr[$4]=$1; start[$4]=$2; end[$4]=$3; original[$4]=$5;
         len[$4]=$3-$2; ids[++n]=$4; next
     }
     {bp[$1 SUBSEP $2]=$3}
     END {
         print "chr","start","end","peak_id","original_name","category","overlap_bp","peak_length_bp","fraction","percent";
         for (i=1; i<=n; i++) {
             id=ids[i]; assigned=0;
             assigned += bp[id SUBSEP "promoter"]+0;
             assigned += bp[id SUBSEP "exon"]+0;
             assigned += bp[id SUBSEP "intron"]+0;
             bp[id SUBSEP "intergenic"] = len[id]-assigned;

             categories[1]="promoter"; categories[2]="exon";
             categories[3]="intron"; categories[4]="intergenic";
             for (j=1; j<=4; j++) {
                 c=categories[j]; x=bp[id SUBSEP c]+0;
                 frac=(len[id]>0 ? x/len[id] : 0);
                 printf "%s\t%d\t%d\t%s\t%s\t%s\t%.0f\t%d\t%.8f\t%.6f\n", \
                        chr[id],start[id],end[id],id,original[id],c,x,len[id],frac,100*frac;
             }
         }
     }' "$OUTDIR/chip/ChIP.peaks.with_ids.bed" \
        "$OUTDIR/tmp/per_peak_category_bp.tsv" \
  > "$OUTDIR/chip/ChIP_per_peak_fractional_annotation.tsv"

# Compact version: only categories with nonzero overlap.
awk 'NR==1 || $7>0' "$OUTDIR/chip/ChIP_per_peak_fractional_annotation.tsv" \
  > "$OUTDIR/chip/ChIP_per_peak_fractional_annotation.nonzero.tsv"

# --------------------- COMBINED PERCENTAGE TABLE -----------------------------
echo "Combining genome, FAIRE, and ChIP percentages..."

awk 'BEGIN{OFS="\t"}
     NR==1 {next}
     tolower($1)!="total" {print "Genome",tolower($1),$2,$3}' \
    "$GENOME_PERCENTAGES" \
    > "$OUTDIR/tmp/genome.long.tsv"

awk 'BEGIN{OFS="\t"} NR>1 && $2!="total" {print $1,$2,$3,$4}' \
    "$OUTDIR/faire/FAIRE_coverage_by_feature.tsv" \
    > "$OUTDIR/tmp/faire.long.tsv"

awk 'BEGIN{OFS="\t"} NR>1 && $2!="total" {print $1,$2,$3,$4}' \
    "$OUTDIR/chip/ChIP_coverage_by_feature.tsv" \
    > "$OUTDIR/tmp/chip.long.tsv"

{
    echo -e "dataset\tcategory\tbp\tpercent"
    cat "$OUTDIR/tmp/genome.long.tsv" \
        "$OUTDIR/tmp/faire.long.tsv" \
        "$OUTDIR/tmp/chip.long.tsv"
} > "$OUTDIR/genome_FAIRE_ChIP_feature_percentages.tsv"

# Wide table useful for plotting and checking.
awk 'BEGIN{OFS="\t"}
     NR>1 {pct[$2 SUBSEP $1]=$4; bp[$2 SUBSEP $1]=$3}
     END {
         print "category","genome_bp","genome_percent","FAIRE_bp","FAIRE_percent","ChIP_bp","ChIP_percent";
         categories[1]="promoter"; categories[2]="exon";
         categories[3]="intron"; categories[4]="intergenic";
         for (i=1; i<=4; i++) {
             c=categories[i];
             print c,bp[c SUBSEP "Genome"]+0,pct[c SUBSEP "Genome"]+0, \
                     bp[c SUBSEP "FAIRE"]+0,pct[c SUBSEP "FAIRE"]+0, \
                     bp[c SUBSEP "ChIP"]+0,pct[c SUBSEP "ChIP"]+0;
         }
     }' "$OUTDIR/genome_FAIRE_ChIP_feature_percentages.tsv" \
  > "$OUTDIR/genome_FAIRE_ChIP_feature_percentages.wide.tsv"

# ------------------------- ENRICHMENT TABLE ----------------------------------
awk 'BEGIN{OFS="\t"}
     NR==1 {next}
     {
         c=$1; genome=$3; accessible=$5; chip=$7;
         fg=(genome>0 ? chip/genome : "NA");
         fa=(accessible>0 ? chip/accessible : "NA");
         lg=(fg=="NA" || fg<=0 ? "NA" : log(fg)/log(2));
         la=(fa=="NA" || fa<=0 ? "NA" : log(fa)/log(2));
         print c,genome,accessible,chip,fg,lg,fa,la;
     }' "$OUTDIR/genome_FAIRE_ChIP_feature_percentages.wide.tsv" \
  | {
      echo -e "category\tgenome_percent\tFAIRE_percent\tChIP_percent\tChIP_over_genome\tlog2_ChIP_over_genome\tChIP_over_FAIRE\tlog2_ChIP_over_FAIRE"
      cat
    } \
  > "$OUTDIR/ChIP_feature_enrichment.tsv"

# ----------------------------- FINAL CHECKS ----------------------------------
# Validation: the promoter + exon + intron + intergenic fractions for every
# ChIP peak must sum to 100% (within a small numerical tolerance). This catches
# annotation or arithmetic inconsistencies before results are accepted.
awk 'NR>1 {sum[$4]+=$10}
     END {
         bad=0;
         for (id in sum) {
             if (sum[id] < 99.999 || sum[id] > 100.001) {
                 print id,sum[id] > "/dev/stderr";
                 bad++
             }
         }
         exit bad>0
     }' "$OUTDIR/chip/ChIP_per_peak_fractional_annotation.tsv" || {
         echo "ERROR: one or more per-peak annotations do not sum to 100%." >&2
         exit 1
     }

rm -f "$OUTDIR/tmp/annotation_self_overlap_check.tsv" \
      "$OUTDIR/tmp/per_peak_category_bp.tsv" \
      "$OUTDIR/tmp/genome.long.tsv" \
      "$OUTDIR/tmp/faire.long.tsv" \
      "$OUTDIR/tmp/chip.long.tsv"

echo
echo "Done. Main outputs:"
echo "  $OUTDIR/genome_FAIRE_ChIP_feature_percentages.tsv"
echo "  $OUTDIR/genome_FAIRE_ChIP_feature_percentages.wide.tsv"
echo "  $OUTDIR/ChIP_feature_enrichment.tsv"
echo "  $OUTDIR/chip/ChIP_per_peak_fractional_annotation.tsv"
echo "  $OUTDIR/chip/ChIP_coverage_by_feature.tsv"
echo "  $OUTDIR/faire/FAIRE_coverage_by_feature.tsv"
