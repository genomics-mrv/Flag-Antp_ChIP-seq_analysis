# Genomic feature coverage analysis

This module quantifies the genomic distribution of Antennapedia ChIP-seq
regions using base-pair coverage.

The analysis uses four mutually exclusive genomic categories:

1. promoter
2. exon
3. intron
4. intergenic

Overlapping annotations are resolved using the hierarchy:

promoter > exon > intron > intergenic

Promoters are defined as ±1 kb around annotated transcription start sites.

## Scripts

### calculate_genomic_feature_percentages.sh

Calculates the fraction of the dm6 genome occupied by promoter, exon,
intron, and intergenic regions after converting the annotations into
mutually exclusive categories.

Usage:

```bash
bash calculate_genomic_feature_percentages.sh \
    ../../resources/annotation/Promoters.bed \
    ../../resources/annotation/Exons.bed \
    ../../resources/annotation/Introns.bed \
    ../../resources/reference/dm6.chrom.sizes


### calculate_feature_coverage.sh

Calculates the number and percentage of unique base pairs from a ChIP-seq
peak set and a FAIRE-seq accessible-chromatin dataset that overlap each
genomic feature.

Individual peaks are not assigned to a single category. A peak spanning
multiple annotations contributes proportionally to each category according
to the number of overlapping base pairs.

bash calculate_feature_coverage.sh \
    ../../resources/annotation/Promoters.bed \
    ../../resources/annotation/Exons.bed \
    ../../resources/annotation/Introns.bed \
    ChIP_peaks.bed \
    FAIRE_peaks.bed \
    genomic_feature_percentages/genomic_feature_percentages.tsv \
    genomic_coverage_annotation

The core overlap calculation uses: bedtools intersect -wo

which reports the exact number of overlapping base pairs between each
dataset interval and each mutually exclusive genomic annotation.

### Requirements:

BEDTools
awk
sort
bash

