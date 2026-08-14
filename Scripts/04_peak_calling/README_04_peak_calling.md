# Peak calling with MACS2

Antp-enriched genomic regions were identified independently in each biological replicate using MACS2 v2.1.1. The corresponding input sample was used as the control.

Because the sequencing libraries were paired-end, peak calling was performed in `BAMPE` mode. In this mode, MACS2 determines genomic fragments directly from the coordinates of properly paired alignments; therefore, no fixed fragment extension size was specified.

Candidate peaks were called using a permissive p-value threshold of `1 × 10⁻⁴` before assessing reproducibility between biological replicates using the Irreproducible Discovery Rate (IDR) framework.

## Requirements

- MACS2 v2.1.1
- Coordinate-sorted and duplicate-removed BAM files
- One BAM file for each Flag-Antp ChIP-seq replicate
- The corresponding input-control BAM file

## Directory structure

```text
scripts/
└── 04_peak_calling/
    ├── README.md
    └── run_macs2.sh
```

The script creates the following results structure:

```text
peak_calling_results/
├── Flag-Antp_rep1/
└── Flag-Antp_rep2/
```

## Usage

The script requires three arguments:

```bash
run_macs2.sh SAMPLE CHIP_BAM INPUT_BAM
```

Where:

- `SAMPLE` is the name assigned to the ChIP-seq replicate.
- `CHIP_BAM` is the duplicate-removed BAM file for that replicate.
- `INPUT_BAM` is the duplicate-removed BAM file for the input control.

Make the script executable:

```bash
chmod +x scripts/04_peak_calling/run_macs2.sh
```

Run MACS2 independently for each biological replicate.

### Replicate 1

```bash
scripts/04_peak_calling/run_macs2.sh \
  Flag-Antp_rep1 \
  bam_results/Flag-Antp_rep1.bam \
  bam_results/Input.bam
```

### Replicate 2

```bash
scripts/04_peak_calling/run_macs2.sh \
  Flag-Antp_rep2 \
  bam_results/Flag-Antp_rep2.bam \
  bam_results/Input.bam
```

Modify the paths if the BAM files are stored in a different directory.

## MACS2 command

For each biological replicate, the script executes:

```bash
macs2 callpeak \
  --treatment "${chip_bam}" \
  --control "${input_bam}" \
  --format BAMPE \
  --gsize dm \
  --name "${sample}" \
  --pvalue 1e-4 \
  --outdir "${output_dir}"
```

The options have the following functions:

- `--treatment`: specifies the Flag-Antp ChIP-seq BAM file.
- `--control`: specifies the input-control BAM file.
- `--format BAMPE`: instructs MACS2 to use paired-end fragment coordinates.
- `--gsize dm`: uses the effective genome-size setting for *Drosophila melanogaster*.
- `--name`: defines the prefix used for the output files.
- `--pvalue 1e-4`: retains candidate peaks with a p-value of `1 × 10⁻⁴` or lower.
- `--outdir`: specifies the output directory.

## Output files

MACS2 generates several files for each replicate, including:

```text
Flag-Antp_rep1_peaks.narrowPeak
Flag-Antp_rep1_peaks.xls
Flag-Antp_rep1_summits.bed
Flag-Antp_rep1_treat_pileup.bdg
Flag-Antp_rep1_control_lambda.bdg
```

Corresponding files are generated for replicate 2.

The principal files used in the subsequent reproducibility analysis are:

```text
Flag-Antp_rep1_peaks.narrowPeak
Flag-Antp_rep2_peaks.narrowPeak
```

These files contain the candidate enriched regions and their statistical scores and are used as input for the subsequent IDR analysis.

## Notes

- Peak calling must be performed independently for each biological replicate.
- The replicate BAM files should not be merged before peak calling.
- The same input-control BAM file can be used for both replicates because this experiment contained one shared input sample.
- The `--pvalue 1e-4` option should not be replaced with `--qvalue 0.1`, because p-value and q-value thresholds are different statistical criteria.
- A fixed fragment extension size should not be specified in `BAMPE` mode because fragment boundaries are obtained directly from the paired-end alignments.

## Reproducibility assessment using IDR

Peaks identified independently in the two biological replicates were compared using the Irreproducible Discovery Rate (IDR) framework. Peaks were ranked using the MACS2 p-value score, corresponding to the `pValue` column of the `narrowPeak` files.

Regions with a global IDR of 0.05 or lower were retained as high-confidence reproducible Flag-Antp binding regions.

### Requirements

- IDR
- MACS2 `narrowPeak` files from both biological replicates
- Matplotlib for generation of the IDR diagnostic plot

### Usage

Make the script executable:

```bash
chmod +x scripts/04_peak_calling/run_idr.sh
```

Run IDR using the peak files from both replicates:

```bash
scripts/04_peak_calling/run_idr.sh \
  peak_calling_results/Flag-Antp_rep1/Flag-Antp_rep1_peaks.narrowPeak \
  peak_calling_results/Flag-Antp_rep2/Flag-Antp_rep2_peaks.narrowPeak \
  Flag-Antp
```

### IDR command

The script executes:

```bash
idr \
  --samples "${rep1_peaks}" "${rep2_peaks}" \
  --input-file-type narrowPeak \
  --output-file-type narrowPeak \
  --rank p.value \
  --idr-threshold 0.05 \
  --soft-idr-threshold 0.05 \
  --output-file "${output_dir}/${output_prefix}_IDR0.05.narrowPeak" \
  --plot
```

The options have the following functions:

- `--samples`: specifies the MACS2 peak files from the two biological replicates.
- `--input-file-type narrowPeak`: identifies the input files as `narrowPeak` files.
- `--output-file-type narrowPeak`: writes the reproducible regions in `narrowPeak` format.
- `--rank p.value`: ranks the peaks using the MACS2 p-value score.
- `--idr-threshold 0.05`: retains peaks with a global IDR of 0.05 or lower.
- `--soft-idr-threshold 0.05`: uses an IDR threshold of 0.05 for reporting summary statistics and generating the diagnostic plot.
- `--output-file`: specifies the name and location of the IDR-filtered peak file.
- `--plot`: generates an IDR diagnostic plot.

### Output files

The IDR analysis produces:

```text
idr_results/
├── Flag-Antp_IDR0.05.narrowPeak
├── Flag-Antp_IDR0.05.narrowPeak.png
└── Flag-Antp_IDR.log
```

The file:

```text
Flag-Antp_IDR0.05.narrowPeak
```

contains the high-confidence reproducible Flag-Antp binding regions that pass the global IDR threshold of 0.05. This analysis produced the final set of 1,188 Flag-Antp binding regions used in the downstream analyses.

The diagnostic plot summarizes the correspondence between peak ranks in the two biological replicates and the relationship between peak rank and IDR.

## Notes

- IDR must be run using the peak files generated independently for the two biological replicates.
- Peaks are ranked using `p.value` because MACS2 peak calling was performed using a p-value threshold.
- In a MACS2 `narrowPeak` file, the `pValue` column contains the `−log10(p-value)` score; higher values therefore represent stronger statistical evidence.
- `--idr-threshold 0.05` filters the output and retains only reproducible peaks.
- `--soft-idr-threshold 0.05` alone would report statistics at this threshold but would not filter the output.

## References

Zhang Y, Liu T, Meyer CA, et al. Model-based analysis of ChIP-Seq (MACS). *Genome Biology*. 2008;9:R137. https://doi.org/10.1186/gb-2008-9-9-r137

Li Q, Brown JB, Huang H, Bickel PJ. Measuring reproducibility of high-throughput experiments. *The Annals of Applied Statistics*. 2011;5:1752–1779. https://doi.org/10.1214/11-AOAS466