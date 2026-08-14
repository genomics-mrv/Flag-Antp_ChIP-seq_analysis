# 07 Antp ChIP-seq dataset comparison

This module compares three Antennapedia ChIP-seq peak sets and annotates selected peak files for downstream gene-based analyses.

## Scripts

### `01_compare_antp_peak_sets.sh`

Compares the following BED files:

- `Flag-Antp_idr.bed`
- `GFP-Antp_Emb_idr.bed`
- `GFP-Antp_WD.bed`

The script first sorts all three BED files by chromosome and genomic coordinate.

It then uses `bedtools intersect` to divide the datasets into mutually exclusive overlap categories:

- Flag-Antp only
- GFP-Antp embryos only
- GFP-Antp wing discs only
- Flag-Antp + GFP-Antp embryos only
- Flag-Antp + GFP-Antp wing discs only
- GFP-Antp embryos + GFP-Antp wing discs only
- shared by all three datasets

`bedtools intersect -u` is used to retain intervals that overlap another dataset, while `bedtools intersect -v` is used to retain intervals that do not overlap another dataset.

No minimum overlap fraction is imposed; any genomic overlap detected by BEDTools is considered an overlap.

### `02_annotate_peaks_homer.sh`

Annotates a BED file using HOMER `annotatePeaks.pl`.

The resulting annotation table can be used for downstream assignment of peaks to genes and subsequent analyses.

## Requirements

- Bash
- BEDTools
- HOMER
- HOMER `dm6` genome annotation
- standard Unix utilities (`sort`, `wc`)
