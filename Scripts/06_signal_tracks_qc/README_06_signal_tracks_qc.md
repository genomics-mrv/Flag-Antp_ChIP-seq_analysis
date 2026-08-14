# Signal-track generation and replicate correlation

CPM-normalized bigWig signal tracks were generated from the processed Flag-Antp BAM files using `bamCoverage` from deepTools. Genome-wide signal values from the two biological replicates were subsequently summarized using `multiBigwigSummary`, and their reproducibility was evaluated by Pearson correlation using `plotCorrelation`.

## Requirements

- deepTools
- Coordinate-sorted, duplicate-removed, and indexed BAM files for both Flag-Antp biological replicates

The following deepTools programs must be available:

```text
bamCoverage
multiBigwigSummary
plotCorrelation
```

## Directory structure

```text
scripts/
└── 06_signal_tracks_qc/
    ├── README.md
    └── run_signal_tracks_qc.sh
```

The script creates the following output directory:

```text
signal_tracks_qc_results/
```

## Usage

Make the script executable:

```bash
chmod +x scripts/06_signal_tracks_qc/run_signal_tracks_qc.sh
```

Run the script from the repository root:

```bash
scripts/06_signal_tracks_qc/run_signal_tracks_qc.sh \
  bam_results/Flag-Antp_rep1.bam \
  bam_results/Flag-Antp_rep2.bam
```

The script requires two arguments:

```text
run_signal_tracks_qc.sh REP1_BAM REP2_BAM
```

Where:

- `REP1_BAM` is the processed BAM file for Flag-Antp biological replicate 1.
- `REP2_BAM` is the processed BAM file for Flag-Antp biological replicate 2.

## CPM-normalized bigWig generation

The processed BAM files were converted into bigWig signal tracks using `bamCoverage`.

For each replicate, the script executes:

```bash
bamCoverage \
  --bam sample.bam \
  --outFileName sample_CPMnorm.bw \
  --outFileFormat bigwig \
  --binSize 25 \
  --normalizeUsing CPM \
  --minMappingQuality 30 \
  --smoothLength 75 \
  --skipNonCoveredRegions \
  --ignoreDuplicates \
  --numberOfProcessors 8
```

The options have the following functions:

- `--bam`: specifies the input BAM file.
- `--outFileName`: specifies the output bigWig file.
- `--outFileFormat bigwig`: writes the signal track in bigWig format.
- `--binSize 25`: calculates coverage in 25-bp genomic bins.
- `--normalizeUsing CPM`: normalizes the signal to counts per million mapped reads.
- `--minMappingQuality 30`: includes only alignments with a mapping quality of at least 30.
- `--smoothLength 75`: smooths the signal using a 75-bp window.
- `--skipNonCoveredRegions`: omits genomic regions without read coverage from the output.
- `--ignoreDuplicates`: excludes duplicate alignments from signal calculation.
- `--numberOfProcessors 8`: uses eight processor cores.

The mapping-quality and duplicate filters provide an additional safeguard, although these filters were already applied during upstream BAM processing.

The resulting CPM-normalized bigWig files can be visualized in genome browsers such as IGV or the UCSC Genome Browser.

## Genome-wide signal summarization

Signal values from the two CPM-normalized bigWig files were summarized across the genome using `multiBigwigSummary`.

```bash
multiBigwigSummary bins \
  --bwfiles \
    Flag-Antp_rep1_CPMnorm.bw \
    Flag-Antp_rep2_CPMnorm.bw \
  --labels \
    Flag-Antp_rep1 \
    Flag-Antp_rep2 \
  --binSize 10000 \
  --numberOfProcessors 8 \
  --outFileName correlation_results.npz \
  --outRawCounts correlation_bin_values.tsv
```

The options have the following functions:

- `bins`: divides the genome into consecutive bins.
- `--bwfiles`: specifies the two CPM-normalized bigWig files.
- `--labels`: assigns descriptive labels to the biological replicates.
- `--binSize 10000`: summarizes the signal in 10-kb genomic bins.
- `--outFileName`: generates a compressed matrix for use with `plotCorrelation`.
- `--outRawCounts`: writes the summarized signal values as a tab-delimited table.
- `--numberOfProcessors 8`: uses eight processor cores.

The 10-kb bin size corresponds to the default bin size used by `multiBigwigSummary` and is stated explicitly in the script to ensure reproducibility.

## Pearson correlation between biological replicates

Reproducibility between the two Flag-Antp biological replicates was evaluated using the Pearson correlation coefficient.

```bash
plotCorrelation \
  --corData correlation_results.npz \
  --corMethod pearson \
  --whatToPlot scatterplot \
  --labels \
    Flag-Antp_rep1 \
    Flag-Antp_rep2 \
  --plotTitle "Flag-Antp replicate correlation" \
  --plotFile replicate_correlation.svg \
  --outFileCorMatrix Pearson_correlation_matrix.tsv
```

The options have the following functions:

- `--corData`: specifies the matrix produced by `multiBigwigSummary`.
- `--corMethod pearson`: calculates the Pearson correlation coefficient.
- `--whatToPlot scatterplot`: generates a scatterplot comparing the two replicates.
- `--labels`: defines the sample names shown in the plot.
- `--plotTitle`: specifies the figure title.
- `--plotFile`: writes the correlation scatterplot in SVG format.
- `--outFileCorMatrix`: writes the calculated Pearson correlation coefficients to a tab-delimited matrix.

Genomic bins containing zero signal were retained because the original analysis did not use the `--skipZeros` option.

## Output files

The script produces:

```text
signal_tracks_qc_results/
├── Flag-Antp_rep1_CPMnorm.bw
├── Flag-Antp_rep2_CPMnorm.bw
├── correlation_results.npz
├── correlation_bin_values.tsv
├── Pearson_correlation_matrix.tsv
├── replicate_correlation.svg
├── Flag-Antp_rep1_bamCoverage.log
├── Flag-Antp_rep2_bamCoverage.log
├── multiBigwigSummary.log
└── plotCorrelation.log
```

The principal output files are:

- `Flag-Antp_rep1_CPMnorm.bw`: CPM-normalized signal track for replicate 1.
- `Flag-Antp_rep2_CPMnorm.bw`: CPM-normalized signal track for replicate 2.
- `replicate_correlation.svg`: Pearson correlation scatterplot.
- `Pearson_correlation_matrix.tsv`: numerical correlation matrix.
- `correlation_bin_values.tsv`: signal values summarized in 10-kb genomic bins.

## Reference

Ramírez F, Ryan DP, Grüning B, et al. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research*. 2016;44:W160–W165. https://doi.org/10.1093/nar/gkw257