# Read alignment with Bowtie2

Paired-end reads were aligned to the *Drosophila melanogaster* reference
genome (dm6) using Bowtie2 in `--very-sensitive` mode. Discordant and
mixed-pair alignments were excluded. Alignments with a mapping quality
(MAPQ) below 30 were removed using SAMtools. The retained alignments were
converted to BAM format and sorted by read name for subsequent mate
processing.

## Requirements

- Bowtie2
- SAMtools
- Bowtie2 index for the dm6 reference genome

## Usage

```bash
run_alignment.sh SAMPLE R1.fastq.gz R2.fastq.gz /path/to/dm6_index