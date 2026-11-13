#!/usr/bin/env Rscript
## Read-only intermediate or saved-result validation; no resume planning.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "reproduction-helpers.R"))
root <- reproduction_root(script_dir)
source(file.path(script_dir, "method_config.R"))
source(file.path(script_dir, "workflow-io.R"))
source(file.path(script_dir, "result-validation.R"))
args <- commandArgs(TRUE)
if (identical(args, "--intermediates")) {
  if (!requireNamespace("ondisc", quietly = TRUE)) stop("ondisc is required to validate intermediate objects")
  directory <- Sys.getenv("FRANGIEH_INTERMEDIATE_DIR", unset = file.path(reproduction_results_root(root), "frangieh-2021", PERTURB_INTERMEDIATE_SUBDIR))
  frangieh_validate_intermediates(directory)
} else if (!length(args)) {
  directory <- Sys.getenv("FRANGIEH_SCORE_RESULTS_DIR", unset = file.path(reproduction_results_root(root), "frangieh-2021", PERTURB_RESULTS_SUBDIR))
  invisible(frangieh_validate_collection(frangieh_result_files(directory)))
  cat("All 20 saved Frangieh replicate outputs are valid.\n")
} else stop("Usage: validate-results.R [--intermediates]")
