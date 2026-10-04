## Quantitative analysis of Flag-Antp ChIP/input enrichment

This analysis compares the ChIP/input enrichment of the 1,188 reproducible Flag-Antp IDR peaks according to two classifications:

- **Motif status:** motif-containing (n = 561) versus motif-free (n = 627).
- **Dataset overlap:** peaks shared with the published GFP-Antp embryonic ChIP-seq dataset (n = 739) versus non-shared peaks (n = 449).

### Signal normalization and enrichment calculation

The two Flag-Antp ChIP-seq biological replicates were merged using SAMtools. ChIP/input ratios were calculated using deepTools `bamCompare`, with CPM normalization, a 10-bp bin size and a pseudocount of 1.

The pseudocount is added to both the ChIP and input signals before calculating their ratio:

\[
\text{ChIP/input ratio} =
\frac{\text{CPM-normalized ChIP signal}+1}
{\text{CPM-normalized input signal}+1}
\]

This prevents division by zero and reduces disproportionately large ratios in regions with very low input coverage. Because the pseudocount has a greater influence on low-coverage regions, the resulting values should be interpreted as normalized, pseudocount-adjusted ChIP/input ratios rather than MACS2 fold-enrichment scores.

The mean ChIP/input ratio across each IDR peak was extracted from the resulting bigWig file using deepTools `multiBigwigSummary` in BED-file mode.

### Statistical analysis and visualization

Descriptive statistics were calculated in R using `dplyr`. Enrichment distributions were compared using two-sided Mann–Whitney U tests without continuity correction (`wilcox.test`, `exact = FALSE`, `correct = FALSE`).

Violin plots with embedded box plots and individual peak values were generated using `ggplot2`. The two panels were combined using `patchwork` and exported as a single SVG figure using `svglite`.

### Results

No significant difference in ChIP/input enrichment was detected between motif-containing and motif-free peaks (median ratios of 1.844 and 1.859, respectively; p = 0.875).

In contrast, peaks shared with the published GFP-Antp embryonic dataset exhibited significantly higher enrichment than non-shared peaks (median ratios of 2.014 and 1.672, respectively; p = 1.25 × 10⁻⁸³).

These comparisons characterize enrichment within the Flag-Antp dataset. They do not constitute a direct quantitative comparison of enrichment between the Flag-Antp and GFP-Antp experiments.