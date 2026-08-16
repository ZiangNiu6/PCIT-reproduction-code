#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
project_root <- normalizePath(file.path(script_dir, ".."))
figure_dir <- file.path(
  project_root, "project-progress", "combine", "figures-and-tables"
)

expected_figures <- c(
  "ccle-typeI-synthY-EIF5A-calibration.pdf",
  "illustrative-example-ci-misspecified.pdf",
  "nb-permutation-signflip-ci-focused-3x3.pdf",
  "pcm_additive.pdf",
  "ccle-typeI-synthY-30targets-typeI-combined.pdf",
  "frangieh-score5-fdr-power.pdf",
  "nb-permutation-signflip-mi-focused-3x3.pdf",
  "cv-locked-tuning-simulation.pdf",
  "ccle-typeI-synthY-30targets-pvalqq.pdf",
  "frangieh-score5-pvalue-histograms.pdf",
  "regularization-robustness-grid.pdf"
)

v4_composites <- c(
  "ccle-typeI-synthY-intro-calibration-power.pdf",
  "ccle-typeI-synthY-30targets-typeI-power-combined.pdf"
)
v4_composite_paths <- file.path(figure_dir, v4_composites)
if (any(file.exists(v4_composite_paths))) {
  missing_v4 <- v4_composites[!file.exists(v4_composite_paths)]
  if (length(missing_v4) > 0L) {
    stop(
      "Incomplete matched-v4 manuscript composite set:\n  ",
      paste(missing_v4, collapse = "\n  "),
      call. = FALSE
    )
  }
  expected_figures <- c(
    expected_figures,
    v4_composites,
    "ccle-alternative-synthY-30targets-rejection.pdf"
  )
}

paths <- file.path(figure_dir, expected_figures)
missing <- expected_figures[!file.exists(paths)]
if (length(missing) > 0L) {
  stop(
    "Missing manuscript figure(s):\n  ",
    paste(missing, collapse = "\n  "),
    call. = FALSE
  )
}

sizes <- file.info(paths)$size
empty <- expected_figures[is.na(sizes) | sizes <= 0]
if (length(empty) > 0L) {
  stop(
    "Empty or unreadable manuscript figure(s):\n  ",
    paste(empty, collapse = "\n  "),
    call. = FALSE
  )
}

summary <- data.frame(
  figure = expected_figures,
  bytes = sizes,
  md5 = unname(tools::md5sum(paths)),
  stringsAsFactors = FALSE
)
cat(sprintf("Validated %d manuscript figures in %s\n", nrow(summary), figure_dir))
print(summary, row.names = FALSE)
