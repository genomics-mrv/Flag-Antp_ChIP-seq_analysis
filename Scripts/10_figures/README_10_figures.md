# Figure-generation scripts

This directory contains the scripts used to generate motif visualizations that were subsequently incorporated into the final figure layouts. Assembly of the complete multi-panel figures, including screenshots and Venn diagrams, was performed manually in Adobe Photoshop and is therefore not represented as a fully scripted workflow.

## Contents

### `plot_antp_motif_probability_logos.R`

Generates the probability sequence logos used in Figure 3C for the three motif classes scanned in the Flag-Antp ChIP-seq peaks:

- AT-rich Antp motif (`TTTAATKA`)
- low-affinity Antp-Exd motif (`TGACNNAY`)
- high-affinity Antp-Exd motif (`TGATNNAY`)

The position-frequency matrices are embedded directly in the script. Small deviations of rounded column totals from 1 are normalized before plotting.

Run from this directory with:

```bash
Rscript plot_antp_motif_probability_logos.R
```

An alternative output directory can be supplied as the first argument:

```bash
Rscript plot_antp_motif_probability_logos.R path/to/output
```

The script creates:

```text
output/Antp_motif_probability_logos.svg
output/Antp_motif_probability_logos.eps
```

Requirements:

- R
- ggseqlogo
- ggplot2
- svglite

Install the required R packages, if necessary, with:

```r
install.packages(c("ggseqlogo", "ggplot2", "svglite"))
```

### `motif_logo_strands.py`

Generates the strand-aware motif-scanning representation used in Figure S5. Unlike the probability logos in Figure 3C, this visualization places motif matches in their genomic context and distinguishes hits on the positive and negative strands. The script represents occurrences of the same three motif classes listed above within selected Flag-Antp peaks associated with the Gene Ontology analysis.

Run with:

```bash
python motif_logo_strands.py
```

Requirements:

- Python 3
- logomaker
- matplotlib
- numpy
- pandas

Install the Python packages, if necessary, with:

```bash
python -m pip install logomaker matplotlib numpy pandas
```

## Notes

- The R script generates probability logos and does not display genomic strand orientation.
- The Python script is strand-aware and was used for the genomic motif-hit representation.
- Deprecation warnings referring to `guides()` or `aes_string()` originate from the current `ggseqlogo` implementation and do not prevent figure generation.
- Final figure assembly and labeling were performed manually outside these scripts.
