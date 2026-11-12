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
source(file.path(script_directory, "..", "..", "reproduction-helpers.R"))
package_root <- reproduction_root(script_file)
arguments <- commandArgs(trailingOnly = TRUE)
if (any(!arguments %in% "--validate-only") || anyDuplicated(arguments)) {
  stop("Supported argument: --validate-only")
}
validate_only <- "--validate-only" %in% arguments
simulation_script <- file.path(
  dirname(script_directory),
  "glm-variance-robustness.R"
)

status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(simulation_script),
    if (validate_only) "--validate-only" else "--plot-only")
)
if (!identical(status, 0L)) {
  stop("GLM variance-robustness plotting failed with status ", status, ".")
}

if (validate_only) {
  quit(save = "no", status = 0L)
}
results_root <- reproduction_results_root(package_root)
output_directory <- file.path(
  results_root,
  "simulation-results",
  "glm-variance-robustness"
)
manuscript_directory <- reproduction_figure_dir(package_root)
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
