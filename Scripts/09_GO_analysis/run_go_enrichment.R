#!/usr/bin/env Rscript

# GO Biological Process enrichment for genes associated with the
# 446 Flag-Antp-specific regions.
#
# Usage:
# Rscript run_go_enrichment.R \
#   FOREGROUND_GENES \
#   UNIVERSE_GENES \
#   OUTPUT_DIRECTORY

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage:",
      "Rscript run_go_enrichment.R",
      "FOREGROUND_GENES UNIVERSE_GENES OUTPUT_DIRECTORY"
    ),
    call. = FALSE
  )
}

foreground_file <- normalizePath(args[[1]], mustWork = TRUE)
universe_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

required_packages <- c(
  "clusterProfiler",
  "org.Dm.eg.db",
  "AnnotationDbi",
  "dplyr",
  "stringr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(stringr)
})

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE)
  genes <- trimws(gsub("\r", "", genes, fixed = TRUE))
  sort(unique(genes[nzchar(genes)]))
}

write_gene_list <- function(genes, filename) {
  writeLines(genes, file.path(output_dir, filename))
}

foreground <- read_gene_list(foreground_file)
universe <- read_gene_list(universe_file)

keytype_to_use <- "SYMBOL"
valid_keys <- AnnotationDbi::keys(
  org.Dm.eg.db,
  keytype = keytype_to_use
)

foreground_valid <- intersect(foreground, valid_keys)
universe_valid <- intersect(universe, valid_keys)
foreground_invalid <- setdiff(foreground, valid_keys)
universe_invalid <- setdiff(universe, valid_keys)
foreground_missing_from_universe <- setdiff(
  foreground_valid,
  universe_valid
)

write_gene_list(foreground_valid, "foreground_valid_genes.txt")
write_gene_list(universe_valid, "universe_valid_genes.txt")
write_gene_list(foreground_invalid, "foreground_invalid_genes.txt")
write_gene_list(universe_invalid, "universe_invalid_genes.txt")
write_gene_list(
  foreground_missing_from_universe,
  "foreground_missing_from_universe.txt"
)

summary_lines <- c(
  paste("Foreground input file:", foreground_file),
  paste("Universe input file:", universe_file),
  paste("Foreground genes after cleaning:", length(foreground)),
  paste("Universe genes after cleaning:", length(universe)),
  paste("Valid foreground genes:", length(foreground_valid)),
  paste("Valid universe genes:", length(universe_valid)),
  paste("Invalid foreground genes:", length(foreground_invalid)),
  paste("Invalid universe genes:", length(universe_invalid)),
  paste(
    "Valid foreground genes absent from universe:",
    length(foreground_missing_from_universe)
  )
)

writeLines(
  summary_lines,
  file.path(output_dir, "GO_analysis_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"), "\n")

if (length(foreground_valid) < 10) {
  stop(
    "Fewer than 10 valid foreground genes remain after validation.",
    call. = FALSE
  )
}

if (length(foreground_missing_from_universe) > 0) {
  stop(
    paste0(
      "The foreground is not fully contained in the universe. See ",
      file.path(output_dir, "foreground_missing_from_universe.txt"),
      "."
    ),
    call. = FALSE
  )
}

ego <- clusterProfiler::enrichGO(
  gene = foreground_valid,
  universe = universe_valid,
  OrgDb = org.Dm.eg.db,
  keyType = keytype_to_use,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ego_df <- if (is.null(ego)) {
  data.frame()
} else {
  as.data.frame(ego)
}

development_neural_pattern <- paste(
  c(
    "embryo",
    "embryonic",
    "brain",
    "neurogenesis",
    "neuron",
    "nervous system",
    "central nervous system",
    "\\bCNS\\b",
    "head development"
  ),
  collapse = "|"
)

if (nrow(ego_df) > 0) {
  development_neural_terms <- ego_df %>%
    filter(
      str_detect(
        Description,
        regex(development_neural_pattern, ignore_case = TRUE)
      )
    )
} else {
  development_neural_terms <- ego_df
}

candidate_genes <- character(0)
if (
  nrow(development_neural_terms) > 0 &&
  "geneID" %in% colnames(development_neural_terms)
) {
  candidate_genes <- sort(
    unique(unlist(strsplit(development_neural_terms$geneID, "/")))
  )
}
candidate_df <- data.frame(Gene = candidate_genes)

write.csv(
  ego_df,
  file.path(output_dir, "Flag-Antp_GO_BP_all.csv"),
  row.names = FALSE
)
write.csv(
  development_neural_terms,
  file.path(
    output_dir,
    "Flag-Antp_GO_BP_development_neural_subset.csv"
  ),
  row.names = FALSE
)
write.csv(
  candidate_df,
  file.path(
    output_dir,
    "Flag-Antp_candidate_genes_development_neural.csv"
  ),
  row.names = FALSE
)

result_summary <- c(
  paste("Enriched GO Biological Process terms:", nrow(ego_df)),
  paste(
    "Developmental/neural terms selected by keyword:",
    nrow(development_neural_terms)
  ),
  paste("Genes contributing to selected terms:", nrow(candidate_df))
)

write(
  result_summary,
  file = file.path(output_dir, "GO_analysis_summary.txt"),
  append = TRUE
)

cat(paste(result_summary, collapse = "\n"), "\n")
cat("Results written to:", output_dir, "\n")
