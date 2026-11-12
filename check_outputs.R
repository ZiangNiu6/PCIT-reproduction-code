#!/usr/bin/env Rscript
# Validate publication figure files and, optionally, the bundled result archive.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_file), "reproduction-helpers.R"))
package_root <- reproduction_root(script_file)

usage <- function() {
  cat(
    "Usage: Rscript --vanilla check_outputs.R [options]\n",
    "  --figure-dir DIR   Figure directory (default: <results-root>/figures-and-tables)\n",
    "  --results-root DIR Result directory (default: LOCAL_PDRT_DATA_DIR or package/results)\n",
    "  --verify-results   Verify the 24 canonical bundled RDS files against artifact-manifest.json\n",
    "  --allow-legacy     Permit the eight legacy PDFs when neither v4 composite exists\n",
    "  --help             Show this help\n",
    "Explicit paths may be absolute or relative to the current working directory.\n",
    "--verify-results checks the released archive's identity; fresh numerical reruns may differ.\n",
    sep = ""
  )
}
arguments <- commandArgs(TRUE)
if (identical(arguments, "--help")) {
  usage()
  quit(save = "no", status = 0L)
}
figure_argument <- NULL
results_argument <- NULL
verify_results <- FALSE
allow_legacy <- FALSE
seen <- character()
while (length(arguments)) {
  option <- arguments[[1L]]
  arguments <- arguments[-1L]
  if (!option %in% c("--figure-dir", "--results-root", "--verify-results", "--allow-legacy") ||
      option %in% seen) {
    stop("Unknown or repeated option: ", option, ". Use --help for usage.", call. = FALSE)
  }
  seen <- c(seen, option)
  if (option %in% c("--figure-dir", "--results-root")) {
    if (!length(arguments) || startsWith(arguments[[1L]], "--") ||
        !nzchar(arguments[[1L]])) {
      stop("Missing directory after ", option, call. = FALSE)
    }
    value <- normalizePath(path.expand(arguments[[1L]]), mustWork = FALSE)
    arguments <- arguments[-1L]
    if (option == "--figure-dir") figure_argument <- value else results_argument <- value
  } else if (option == "--verify-results") {
    verify_results <- TRUE
  } else {
    allow_legacy <- TRUE
  }
}
results_root <- if (is.null(results_argument)) {
  reproduction_results_root(package_root)
} else results_argument
figure_dir <- if (!is.null(figure_argument)) {
  figure_argument
} else if (!is.null(results_argument)) {
  file.path(results_root, "figures-and-tables")
} else reproduction_figure_dir(package_root)

baseline_figures <- c(
  "ccle-typeI-synthY-EIF5A-calibration.pdf",
  "illustrative-example-ci-misspecified.pdf",
  "nb-permutation-signflip-ci-focused-3x3.pdf",
  "pcm_additive.pdf",
  "ccle-typeI-synthY-30targets-typeI-combined.pdf",
  "frangieh-score5-fdr-power.pdf",
  "nb-permutation-signflip-mi-focused-3x3.pdf",
  "frangieh-score5-pvalue-histograms.pdf"
)
v4_composites <- c(
  "ccle-typeI-synthY-intro-calibration-power.pdf",
  "ccle-typeI-synthY-30targets-typeI-power-combined.pdf"
)
expected_figures <- c(baseline_figures, v4_composites)
if (allow_legacy && !any(file.exists(file.path(figure_dir, v4_composites)))) {
  expected_figures <- baseline_figures
}
paths <- file.path(figure_dir, expected_figures)
missing <- expected_figures[!file.exists(paths)]
if (length(missing)) {
  stop("Missing publication figure(s):\n  ", paste(missing, collapse = "\n  "),
       call. = FALSE)
}
present_figures <- list.files(
  figure_dir, pattern = "[.](pdf|png|jpe?g|eps|svg)$", ignore.case = TRUE
)
unexpected <- setdiff(present_figures, expected_figures)
if (length(unexpected)) {
  stop("Unexpected publication figure(s):\n  ",
       paste(unexpected, collapse = "\n  "),
       "\nKeep diagnostics with their analysis results instead.", call. = FALSE)
}
sizes <- file.info(paths)$size
empty <- expected_figures[is.na(sizes) | sizes < 5L | file.info(paths)$isdir]
if (length(empty)) {
  stop("Empty or unreadable publication figure(s):\n  ",
       paste(empty, collapse = "\n  "), call. = FALSE)
}
is_pdf <- vapply(paths, function(path) {
  tryCatch({
    connection <- file(path, "rb")
    on.exit(close(connection))
    identical(readBin(connection, what = "raw", n = 5L), charToRaw("%PDF-"))
  }, error = function(e) FALSE)
}, logical(1L))
if (any(!is_pdf)) {
  stop("Invalid PDF header in publication figure(s):\n  ",
       paste(expected_figures[!is_pdf], collapse = "\n  "), call. = FALSE)
}

if (verify_results) {
  .libPaths(c(reproduction_library_paths(package_root), .libPaths()))
  required_packages <- c("jsonlite", "digest")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Restore the pinned environment to verify results; missing package(s): ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  manifest_file <- file.path(results_root, "artifact-manifest.json")
  if (!file.exists(manifest_file)) {
    stop("Bundled-result manifest not found: ", manifest_file, call. = FALSE)
  }
  manifest <- jsonlite::fromJSON(manifest_file, simplifyVector = FALSE)
  if (!is.list(manifest) || !identical(manifest$schema, "pcit_bundled_results_v1") ||
      !is.list(manifest$artifacts) || !length(manifest$artifacts)) {
    stop("Malformed or unsupported bundled-result manifest.", call. = FALSE)
  }
  required_artifacts <- c(
    "results/simulation-results/glm-variance-robustness/illustrative-example-results.rds",
    "results/simulation-results/nb-permutation-signflip-continuous-dispersion-estimated-theta/nb-permutation-signflip_results.rds",
    "results/simulation-results/pcm-additive-matched-dcrt/pcm-additive_results.rds",
    "results/ccle/ccle-typeI-synthY-30targets-matched-v4-results.rds",
    sprintf("results/frangieh-2021/score-calibration-nonpositive-nulls/score_fdr_rep%03d.rds", 1:20)
  )
  entries_valid <- vapply(manifest$artifacts, function(entry) {
    is.list(entry) && is.character(entry$path) && length(entry$path) == 1L &&
      !is.na(entry$path) && startsWith(entry$path, "results/") &&
      !grepl("\\\\|(^|/)[.][.]?(/|$)", entry$path) &&
      is.numeric(entry$bytes) && length(entry$bytes) == 1L &&
      is.finite(entry$bytes) && entry$bytes > 0 && entry$bytes == floor(entry$bytes) &&
      is.character(entry$sha256) && length(entry$sha256) == 1L &&
      !is.na(entry$sha256) && grepl("^[a-fA-F0-9]{64}$", entry$sha256)
  }, logical(1L))
  if (!all(entries_valid)) stop("Malformed bundled-result manifest entry.", call. = FALSE)
  artifact_names <- vapply(manifest$artifacts, `[[`, character(1L), "path")
  if (anyDuplicated(artifact_names) || !setequal(artifact_names, required_artifacts)) {
    stop("Manifest must list exactly the 24 canonical bundled result files.", call. = FALSE)
  }
  normalized_root <- normalizePath(results_root, mustWork = TRUE)
  for (entry in manifest$artifacts) {
    path <- file.path(results_root, sub("^results/", "", entry$path))
    if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
      stop("Missing bundled result: ", entry$path, call. = FALSE)
    }
    normalized_path <- normalizePath(path, mustWork = TRUE)
    if (!startsWith(normalized_path, paste0(normalized_root, .Platform$file.sep))) {
      stop("Bundled result resolves outside its result directory: ", entry$path,
           call. = FALSE)
    }
    actual_sha <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
    if (file.info(path)$size != entry$bytes ||
        !identical(tolower(actual_sha), tolower(entry$sha256))) {
      stop("Bundled result checksum/size mismatch: ", entry$path, call. = FALSE)
    }
  }
  cat(sprintf("Verified SHA-256 and sizes of %d canonical bundled result files.\n",
              length(manifest$artifacts)))
}
summary <- data.frame(
  figure = expected_figures, bytes = sizes,
  md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE
)
cat(sprintf("Validated %d publication PDF files in %s\n", nrow(summary), figure_dir))
print(summary, row.names = FALSE)
