# Flag-Antp ChIP-seq analysis

Reproducible bioinformatics workflow for the analysis of Antennapedia (Antp)
ChIP-seq data generated from an endogenously Flag-tagged *Drosophila
melanogaster* allele.

The experiment profiles Flag-Antp occupancy in embryos collected 3–7 hours
after egg laying (AEL), corresponding to embryonic stages 5–11.5. The dataset
contains two biological ChIP-seq replicates and one shared input control, all
sequenced as paired-end (PE150) libraries.

Raw sequencing data are available from the NCBI Gene Expression Omnibus:

**[GSE318263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318263)**

This repository accompanies the associated *Data in Brief* dataset article and
contains the scripts, small reference files, and documentation needed to
reproduce the computational analyses.

## Analysis overview

The workflow performs:

1. Download of Flag-Antp and published comparison datasets.
2. FASTQ quality assessment with FastQC.
3. Paired-end alignment to dm6 with Bowtie2.
4. Alignment filtering, mate processing, coordinate sorting, and duplicate
   removal with SAMtools.
5. Replicate-level peak calling with MACS2 in `BAMPE` mode.
6. Identification of reproducible binding regions using IDR.
7. Genomic-feature annotation and base-pair coverage analysis.
8. Generation of CPM-normalized signal tracks and replicate-correlation QC.
9. Antp motif scanning with FIMO.
10. Comparison with published embryonic and wing-disc GFP-Antp datasets.
11. GO Biological Process enrichment analysis.
12. Generation of motif-based figure panels.

The principal analysis parameters include:

- Reference genome: dm6
- Minimum alignment quality: MAPQ 30
- MACS2 input format: paired-end BAM (`BAMPE`)
- MACS2 candidate-peak threshold: `p ≤ 1 × 10⁻⁴`
- Reproducibility threshold: `IDR ≤ 0.05`
- FIMO motif-occurrence threshold: `p < 1 × 10⁻³`

The analysis produced 1,188 high-confidence reproducible Flag-Antp binding
regions.

## Repository structure

```text
Flag-Antp_ChIPseq_analysis/
├── README.md
├── LICENCE.md
├── Documents/
│   └── Calculation of genomic percentages.docx
├── Environment/
│   ├── README_enviroment.md
│   └── Software_versions.xlsx
├── Resources/
│   ├── annotation/
│   │   ├── Promoters.bed
│   │   ├── Exons.bed
│   │   └── Introns.bed
│   ├── motifs/
│   │   └── Antp_motif_matrices.txt
│   └── data/
└── Scripts/
    ├── 00_download/
    ├── 01_qc/
    ├── 02_alignment/
    ├── 03_bam_processing/
    ├── 04_peak_calling/
    ├── 05_genomic_features/
    ├── 06_signal_tracks_qc/
    ├── 07_motif_analysis/
    ├── 08_dataset_comparison/
    ├── 09_GO_analysis/
    └── 10_figures/
```

Large sequencing and alignment files are not stored in the repository. The
download scripts retrieve the public data required by the workflow.

## Workflow modules

Each numbered directory is a documented analysis module. Run the modules in
numerical order unless starting from an existing intermediate file.

| Step | Directory | Purpose |
|---:|---|---|
| 00 | [`00_download`](Scripts/00_download/README_00_download.md) | Download Flag-Antp data, published comparison datasets, dm6 resources, blacklist regions, and liftOver files |
| 01 | [`01_qc`](Scripts/01_qc/README_01_qc.md) | Assess raw FASTQ quality with FastQC |
| 02 | [`02_alignment`](Scripts/02_alignment/README_02_alignment.md) | Align paired-end reads to dm6 and retain MAPQ ≥30 alignments |
| 03 | [`03_bam_processing`](Scripts/03_bam_processing/README_03_bam_processing.md) | Add mate information, coordinate-sort alignments, remove duplicates, and index BAM files |
| 04 | [`04_peak_calling`](Scripts/04_peak_calling/README_04_peak_calling.md) | Call replicate peaks with MACS2 and identify reproducible regions using IDR |
| 05 | [`05_genomic_features`](Scripts/05_genomic_features/README_05_genomic_features.md) | Classify regions and calculate genomic-feature representation and base-pair coverage |
| 06 | [`06_signal_tracks_qc`](Scripts/06_signal_tracks_qc/README_06_signal_tracks_qc.md) | Generate CPM-normalized bigWigs and evaluate replicate correlation |
| 07 | [`07_motif_analysis`](Scripts/07_motif_analysis/README_07_motif_analysis.md) | Extract summit-centered sequences and identify Antp motif occurrences with FIMO |
| 08 | [`08_dataset_comparison`](Scripts/08_dataset_comparison/README_08_dataset_comparison.md) | Compare Flag-Antp regions with embryonic and wing-disc GFP-Antp datasets and annotate target genes |
| 09 | [`09_GO_analysis`](Scripts/09_GO_analysis/README_09_GO_analysis.md) | Test GO Biological Process enrichment among genes associated with Flag-Antp-specific regions |
| 10 | [`10_figures`](Scripts/10_figures/README_10_figures.md) | Generate motif probability logos and strand-aware motif representations |

The README inside each module describes its inputs, dependencies, command-line
usage, parameters, and expected outputs.

## Data sources

The download module retrieves or documents the following resources:

- Flag-Antp ChIP-seq replicates and input control from GEO accession
  [GSE318263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318263).
- Published embryonic GFP-Antp ChIP-seq data.
- Published wing-disc GFP-Antp ChIP-seq data.
- Published embryonic and wing-disc FAIRE-seq datasets used to represent
  accessible chromatin.
- The dm6 *Drosophila melanogaster* genome assembly.
- dm6 blacklist regions.
- The dm3-to-dm6 liftOver chain used where coordinate conversion was required.

Accession numbers, source URLs, and download commands are provided in
[`Scripts/00_download`](Scripts/00_download/README_00_download.md).

## Reference resources

### Genomic annotations

`Resources/annotation` contains the promoter, exon, and intron intervals used
for genomic-feature analyses. Features were treated as mutually exclusive using
the following hierarchy:

```text
promoter > exon > intron > intergenic
```

Promoters were defined as regions extending 1 kb upstream and downstream of
annotated transcription start sites.

### Antp motifs

`Resources/motifs/Antp_motif_matrices.txt` contains the three probability
matrices used for motif scanning and visualization:

- Antp AT-rich motif (`TTTAATKA`)
- Antp–Exd low-affinity motif (`TGACNNAY`)
- Antp–Exd high-affinity motif (`TGATNNAY`)

The low- and high-affinity Antp–Exd motifs were described by Kribelbauer and
colleagues.

## GO enrichment design

GO Biological Process enrichment was evaluated for genes associated with the
446 Flag-Antp-specific regions. Genes associated with the combined set of
Flag-Antp-specific and embryonic GFP-Antp-shared regions were used as the
background universe.

The analysis therefore tests for biological processes overrepresented among
Flag-Antp-specific targets relative to the broader Antp-associated gene set,
not relative to every annotated gene in the *Drosophila* genome.

## Figure generation

The repository includes scripts used to generate:

- Probability logos for the three Antp motifs in Figure 3C using `ggseqlogo`.
- Strand-aware motif representations for selected GO-associated peaks in
  Figure S5 using Python and Logomaker.

Other figures, including genome-browser panels and Venn diagrams, were
assembled manually from the documented analysis outputs.

## Software environment

The workflow uses command-line tools, R packages, and Python packages. Exact
versions used for the analysis are recorded in:

- [`Environment/Software_versions.xlsx`](Environment/Software_versions.xlsx)
- [`Environment/README_enviroment.md`](Environment/README_enviroment.md)

Major dependencies include FastQC, Bowtie2, SAMtools, BEDTools, MACS2, IDR,
deepTools, MEME Suite/FIMO, HOMER, R, `clusterProfiler`, `org.Dm.eg.db`,
`ggseqlogo`, Python, Matplotlib, pandas, NumPy, and Logomaker.

## Citation

If you use this repository or dataset, please cite the associated *Data in
Brief* article and the GEO dataset:

> Flag-Antp ChIP-seq in *Drosophila melanogaster* embryos. GEO accession
> GSE318263.

The complete article citation and DOI should be added here when available.

## License

This repository is distributed under the terms described in
[`LICENCE.md`](LICENCE.md).

## Contact

Questions, problems, and reproducibility issues can be submitted through the
GitHub repository's issue tracker.
