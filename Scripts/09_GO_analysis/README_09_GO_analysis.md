# Gene Ontology enrichment analysis

Gene Ontology (GO) enrichment analysis was performed to identify biological processes associated with genes linked to the Flag-Antp-specific binding regions.

Genes associated with the 446 Flag-Antp-specific regions were used as the foreground gene set. Genes associated with the combined set of Flag-Antp-specific and embryonic GFP-Antp-shared regions were used as the background universe.

This analysis therefore identifies GO terms that are overrepresented among genes associated with the Flag-Antp-specific regions relative to the broader set of Antp-associated genes. It is not a comparison against all annotated genes in the *Drosophila melanogaster* genome.

## Directory structure

```text
scripts/
└── 09_go_analysis/
    ├── README.md
    └── run_go_enrichment.R

data/
└── go_analysis/
    ├── Flag-Antp_specific_genes.txt
    └── Flag-GFP_Antp_universe_genes.txt
```

## Input files

### Foreground gene set

```text
Flag-Antp_specific_genes.txt
```

This file contains the nonredundant gene symbols associated with the 446 Flag-Antp-specific binding regions.

### Background universe

```text
Flag-GFP_Antp_universe_genes.txt
```

This file contains the nonredundant gene symbols associated with the combined set of Flag-Antp-specific and embryonic GFP-Antp-shared binding regions.

The foreground genes must be contained within the background universe.

Both input files must contain one official *Drosophila melanogaster* gene symbol per line:

```text
Antp
engrailed
mid
ptc
```

Blank lines and duplicated gene symbols are removed automatically by the R script.

## Requirements

The analysis requires R and the following packages:

- `clusterProfiler`
- `org.Dm.eg.db`
- `AnnotationDbi`
- `dplyr`
- `stringr`

The Bioconductor packages can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler",
  "org.Dm.eg.db",
  "AnnotationDbi"
))
```

The CRAN packages can be installed with:

```r
install.packages(c(
  "dplyr",
  "stringr"
))
```

## Usage

Run the analysis from the repository root:

```bash
Rscript scripts/09_go_analysis/run_go_enrichment.R \
  data/go_analysis/Flag-Antp_specific_genes.txt \
  data/go_analysis/Flag-GFP_Antp_universe_genes.txt \
  go_analysis_results
```

The script requires three arguments:

```text
run_go_enrichment.R FOREGROUND_GENES UNIVERSE_GENES OUTPUT_DIRECTORY
```

Where:

- `FOREGROUND_GENES` is the gene list associated with the 446 Flag-Antp-specific regions.
- `UNIVERSE_GENES` is the background set containing genes associated with the combined Flag-Antp-specific and embryonic GFP-Antp-shared regions.
- `OUTPUT_DIRECTORY` specifies where the results will be written.

## Gene identifier processing

The input lists are interpreted as official gene symbols using:

```r
keyType = "SYMBOL"
```

Gene symbols are validated against the `org.Dm.eg.db` annotation database.

The script reports:

- Number of genes in the original foreground.
- Number of genes in the original universe.
- Number of valid foreground gene symbols.
- Number of valid universe gene symbols.
- Foreground genes not found in the annotation database.
- Universe genes not found in the annotation database.
- Foreground genes absent from the background universe.

The analysis stops if any valid foreground genes are absent from the background universe, because the foreground must be a subset of the universe.

## GO enrichment

GO Biological Process enrichment is performed using `enrichGO` from the `clusterProfiler` package:

```r
ego <- enrichGO(
  gene          = foreground,
  universe      = universe,
  OrgDb         = org.Dm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
```

The analysis uses:

- GO ontology: Biological Process (`BP`)
- Identifier type: official gene symbol (`SYMBOL`)
- Multiple-testing correction: Benjamini–Hochberg (`BH`)
- P-value cutoff: `0.05`
- Q-value cutoff: `0.20`

The reported enrichment is conditional on the specified Antp-associated background universe.

## Developmental and neural terms

In addition to the complete GO enrichment results, the script generates a descriptive subset of enriched terms related to embryonic development and nervous-system biology.

Terms are selected when their descriptions contain one or more of the following keywords:

```text
embryo
embryonic
brain
neurogenesis
neuron
nervous system
central nervous system
CNS
head development
```

This keyword-based selection is applied after GO enrichment. It is intended to facilitate interpretation and does not represent an independent statistical enrichment analysis.

All statistically enriched terms remain available in the complete results table.

## Output files

The script creates:

```text
go_analysis_results/
├── Flag-Antp_GO_BP_all.csv
├── Flag-Antp_GO_BP_development_neural_subset.csv
├── Flag-Antp_candidate_genes_development_neural.csv
├── foreground_valid_genes.txt
├── universe_valid_genes.txt
├── foreground_invalid_genes.txt
├── universe_invalid_genes.txt
├── GO_analysis_summary.txt
└── R_sessionInfo.txt
```

### Complete enrichment results

```text
Flag-Antp_GO_BP_all.csv
```

Contains all GO Biological Process terms passing the enrichment thresholds.

### Developmental and neural subset

```text
Flag-Antp_GO_BP_development_neural_subset.csv
```

Contains enriched terms selected by the developmental and neural keyword filter.

### Candidate genes

```text
Flag-Antp_candidate_genes_development_neural.csv
```

Contains the nonredundant foreground genes contributing to the selected developmental and neural GO terms.

### Identifier-validation files

```text
foreground_valid_genes.txt
universe_valid_genes.txt
foreground_invalid_genes.txt
universe_invalid_genes.txt
```

These files document which input gene symbols were recognized by `org.Dm.eg.db`.

### Analysis summary

```text
GO_analysis_summary.txt
```

Reports gene-list sizes, identifier-validation results, and the number of enriched GO terms.

### Software information

```text
R_sessionInfo.txt
```

Records the R version, operating system, and package versions used for the analysis.

## Interpretation

A significantly enriched GO term indicates that the corresponding biological process is represented more frequently among genes associated with the 446 Flag-Antp-specific regions than expected relative to the specified Antp-associated background universe.

Because the background is restricted to Antp-associated genes, the results should not be interpreted as enrichment relative to the complete *Drosophila melanogaster* genome.

## References

Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS*. 2012;16:284–287. https://doi.org/10.1089/omi.2011.0118

Carlson M. `org.Dm.eg.db`: Genome-wide annotation for *Drosophila melanogaster*. Bioconductor annotation package.