#!/usr/bin/env Rscript

# Reproduce the probability sequence logos used in Figure 3C.
#
# Usage:
#   Rscript plot_antp_motif_probability_logos.R [output_directory]
#
# If no output directory is supplied, files are written to "output".

required_packages <- c("ggseqlogo", "ggplot2", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pfm_text <- "
>Kribelbauer2020_AntpExd_high_TGATNNAY
0.041 0.037 0.900 0.041 0.291 0.291 0.900 0.058
0.029 0.026 0.029 0.029 0.209 0.209 0.029 0.450
0.029 0.900 0.029 0.029 0.209 0.209 0.029 0.042
0.900 0.037 0.041 0.900 0.291 0.291 0.041 0.450

>Kribelbauer2020_AntpExd_low_TGACNNAY
0.041 0.037 0.900 0.037 0.291 0.291 0.900 0.058
0.029 0.026 0.029 0.900 0.209 0.209 0.029 0.450
0.029 0.900 0.029 0.026 0.209 0.209 0.029 0.042
0.900 0.037 0.041 0.037 0.291 0.291 0.041 0.450

>Antp_SOLEXA_FBgn0000095
0.209 0.141 0.068 0.794 0.962 0.001 0.065 0.581
0.257 0.195 0.000 0.048 0.038 0.075 0.033 0.013
0.206 0.085 0.010 0.010 0.000 0.037 0.462 0.362
0.328 0.579 0.922 0.149 0.000 0.887 0.440 0.044
"

read_pfm_block <- function(lines4) {
  if (length(lines4) != 4L) {
    stop("Each motif must contain exactly four rows (A, C, G and T).")
  }

  matrix_rows <- lapply(
    lines4,
    function(line) as.numeric(strsplit(trimws(line), "\\s+")[[1L]])
  )

  if (any(vapply(matrix_rows, function(x) anyNA(x), logical(1)))) {
    stop("A motif matrix contains a non-numeric value.")
  }

  matrix_widths <- vapply(matrix_rows, length, integer(1))
  if (length(unique(matrix_widths)) != 1L) {
    stop("All rows within a motif matrix must have the same length.")
  }

  matrix <- do.call(rbind, matrix_rows)
  rownames(matrix) <- c("A", "C", "G", "T")

  column_totals <- colSums(matrix)
  if (any(column_totals <= 0)) {
    stop("Every motif position must have a positive probability total.")
  }

  # Correct small deviations from 1 that result from rounded probabilities.
  sweep(matrix, 2L, column_totals, "/")
}

parse_pfms <- function(text) {
  lines <- trimws(unlist(strsplit(text, "\n", fixed = TRUE)))
  lines <- lines[nzchar(lines)]

  motifs <- list()
  i <- 1L

  while (i <= length(lines)) {
    if (!startsWith(lines[[i]], ">")) {
      stop("Expected a motif name beginning with '>' at line ", i, ".")
    }

    if (i + 4L > length(lines)) {
      stop("Incomplete motif matrix for ", lines[[i]], ".")
    }

    motif_name <- sub("^>", "", lines[[i]])
    motifs[[motif_name]] <- read_pfm_block(lines[(i + 1L):(i + 4L)])
    i <- i + 5L
  }

  motifs
}

motifs <- parse_pfms(pfm_text)

# Keep the motif order and terminology consistent with the manuscript.
motifs <- motifs[c(
  "Antp_SOLEXA_FBgn0000095",
  "Kribelbauer2020_AntpExd_low_TGACNNAY",
  "Kribelbauer2020_AntpExd_high_TGATNNAY"
)]
names(motifs) <- c(
  "AT-rich Antp motif (TTTAATKA)",
  "Low-affinity Antp-Exd motif (TGACNNAY)",
  "High-affinity Antp-Exd motif (TGATNNAY)"
)

logo_plot <- ggseqlogo::ggseqlogo(
  motifs,
  method = "prob",
  ncol = 3
) +
  ggplot2::labs(
    title = "Antp motifs (probability logos)",
    x = "Position",
    y = "Probability"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank()
  )

svg_file <- file.path(output_dir, "Antp_motif_probability_logos.svg")
eps_file <- file.path(output_dir, "Antp_motif_probability_logos.eps")

ggplot2::ggsave(
  filename = svg_file,
  plot = logo_plot,
  device = svglite::svglite,
  width = 9,
  height = 3,
  units = "in",
  bg = "white"
)

save_eps <- function(filename, plot, width, height) {
  grDevices::postscript(
    file = filename,
    width = width,
    height = height,
    onefile = FALSE,
    horizontal = FALSE,
    paper = "special"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

save_eps(eps_file, logo_plot, width = 9, height = 3)

message("Created: ", normalizePath(svg_file, mustWork = FALSE))
message("Created: ", normalizePath(eps_file, mustWork = FALSE))

