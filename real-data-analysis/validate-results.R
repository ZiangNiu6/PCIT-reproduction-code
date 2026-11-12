#!/usr/bin/env Rscript
## Read-only resume planning and intermediate checks; stdout contains only IDs.
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
  directory <- Sys.getenv("FRANGIEH_INTERMEDIATE_DIR", unset = file.path(reproduction_results_root(root), "frangieh-2021", "intermediate-data"))
  frangieh_validate_intermediates(directory, script_dir = script_dir)
  quit(save = "no", status = 0L)
}
force <- "--force" %in% args
args <- args[args != "--force"]
if (length(args) != 2L || args[[1L]] != "--missing") stop("Usage: validate-results.R --missing N [--force] | --intermediates")
n <- suppressWarnings(as.numeric(args[[2L]]))
if (length(n) != 1L || !is.finite(n) || n != as.integer(n) || !n %in% 1:20) stop("N must be an integer in 1:20")
directory <- Sys.getenv("FRANGIEH_SCORE_RESULTS_DIR", unset = file.path(reproduction_results_root(root), "frangieh-2021", PERTURB_RESULTS_SUBDIR))
files <- frangieh_result_files(directory)
if (force) {
  selected <- sprintf("score_fdr_rep%03d.rds", seq_len(as.integer(n)))
  invisible(frangieh_validate_collection(files[!basename(files) %in% selected], require_complete = FALSE))
}
if (!force) {
  intermediate <- Sys.getenv("FRANGIEH_INTERMEDIATE_DIR")
  if (!nzchar(intermediate)) intermediate <- file.path(reproduction_results_root(root), "frangieh-2021", "intermediate-data")
  intermediate <- frangieh_intermediate_directory(intermediate)
  if (!all(file.exists(file.path(intermediate, frangieh_intermediate_files)))) intermediate <- NULL
  invisible(frangieh_validate_collection(files, require_complete = FALSE, intermediate_dir = intermediate))
}
missing <- seq_len(as.integer(n))
if (!force) missing <- missing[!sprintf("score_fdr_rep%03d.rds", missing) %in% basename(files)]
if (length(missing)) cat(paste(missing, collapse = " "), "\n", sep = "")
