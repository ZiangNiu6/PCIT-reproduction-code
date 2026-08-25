#!/usr/bin/env Rscript

# Rebuild the two manuscript GLM-score figures from the saved null simulations.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", script_argument))
script_directory <- dirname(script_file)
project_root <- normalizePath(file.path(script_directory, "..", "..", ".."))
simulation_script <- file.path(
  dirname(script_directory),
  "glm-variance-robustness.R"
)

status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(simulation_script), "--plot-only")
)
if (!identical(status, 0L)) {
  stop("GLM variance-robustness plotting failed with status ", status, ".")
}

results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
output_directory <- file.path(
  results_root,
  "simulation-results",
  "glm-variance-robustness"
)
manuscript_directory <- file.path(
  project_root,
  "manuscript",
  "figures-and-tables"
)
dir.create(manuscript_directory, recursive = TRUE, showWarnings = FALSE)

figure_names <- c(
  "illustrative-example-ci-misspecified.pdf"
)
for (figure_name in figure_names) {
  source_file <- file.path(output_directory, figure_name)
  target_file <- file.path(manuscript_directory, figure_name)
  if (!file.exists(source_file)) {
    stop("Expected figure was not created: ", source_file)
  }
  if (!file.copy(source_file, target_file, overwrite = TRUE)) {
    stop("Could not copy manuscript figure to ", target_file)
  }
  message("Copied: ", normalizePath(target_file))
}
