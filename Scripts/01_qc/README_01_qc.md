# Quality control with FastQC

Raw paired-end FASTQ files were assessed using FastQC to evaluate sequencing
quality before read alignment.

## Requirements

- FastQC
- Paired-end FASTQ files compressed with gzip (`*.fastq.gz`)

## Usage

Run the script from the directory containing the FASTQ files:

```bash
/path/to/Antp_ChIPseq_analysis/scripts/01_qc/run_fastqc.sh