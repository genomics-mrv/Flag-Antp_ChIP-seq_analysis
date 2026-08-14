# BAM processing and duplicate removal

Name-sorted paired-end alignments were processed using SAMtools. Mate
information and mate-related alignment tags were added with `samtools
fixmate -m`. The alignments were then sorted by genomic coordinate, and
PCR duplicates were identified and removed using `samtools markdup -r`.
The resulting BAM files were indexed, and alignment summary statistics
were generated using `samtools flagstat`.

## Usage

```bash
run_bam_processing.sh SAMPLE INPUT_NAMESORT_BAM