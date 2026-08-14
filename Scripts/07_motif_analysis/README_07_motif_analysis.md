# Antp motif analysis

This module identifies previously described Antennapedia DNA-binding
motifs within high-confidence Flag-Antp ChIP-seq regions.

## Workflow

1. Generate 200-bp regions centered on IDR peak summits.
2. Extract genomic sequences from the dm6 reference genome.
3. Scan sequences using FIMO.
4. Classify peaks according to the presence or absence of the three
   Antp motif models.
5. Generate motif co-occurrence and sequence-logo visualizations.

## Motifs

Three motif models are used:

- Antp TTTAATKA motif (`Antp_SOLEXA_FBgn0000095`)
- Antp-Exd high-affinity motif (`TGATNNAY`)
- Antp-Exd low-affinity motif (`TGACNNAY`)

## Requirements

- BEDTools
- MEME Suite / FIMO
- R
- ComplexUpset
- ggseqlogo

## Genome assembly

All genomic coordinates and extracted sequences correspond to dm6.