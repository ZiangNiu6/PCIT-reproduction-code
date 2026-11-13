#!/usr/bin/env Rscript
## Filesystem/collection checks use only temporary synthetic q10 fixtures.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "test-result-contracts.R"))
run_io_checks <- function() {
  sandbox <- tempfile("frangieh-q10-io-");dir.create(sandbox);on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  inputs <- file.path(sandbox, "inputs");object <- make_frangieh_fixture(inputs)
  files <- file.path(sandbox, sprintf("score_fdr_rep%03d.rds", 1:20))
  for (i in 1:20) { copy <- object;copy$config$rep_idx <- i;copy$config$seed <- 1000L+i;copy$result$rep <- i;copy$result$seed <- 1000L+i;saveRDS(copy, files[[i]]) }
  hashes <- tools::md5sum(files);set.seed(123);rng <- .Random.seed
  stopifnot(length(frangieh_validate_collection(frangieh_result_files(sandbox))) == 20L, identical(.Random.seed, rng))
  expect_error(frangieh_validate_collection(files[-20]));expect_error(validate_rep(object, files[[2]]))
  unexpected <- file.path(sandbox, "score_fdr_rep021.rds");file.copy(files[[1]], unexpected);expect_error(frangieh_result_files(sandbox));unlink(unexpected)
  atomic <- file.path(sandbox, "atomic.rds");frangieh_atomic_rds("old", atomic);expect_error(frangieh_atomic_rds("replacement", atomic))
  stopifnot(identical(readRDS(atomic), "old"));frangieh_atomic_rds("new", atomic, replace = TRUE);stopifnot(identical(readRDS(atomic), "new"))
  lock <- paste0(atomic, ".publish-lock");dir.create(lock);expect_error(frangieh_atomic_rds("concurrent", atomic, replace = TRUE));unlink(lock, recursive = TRUE)
  # A first intermediate generation is readable; no generation replacement
  # or resume path is part of the numerical workflow.
  bundle <- file.path(sandbox, "bundle")
  folder <- file.path(bundle, ".generations", "generation-first")
  dir.create(folder, recursive = TRUE)
  saveRDS("unchanged", file.path(folder, "sentinel.rds"))
  pointer <- file.path(bundle, "current-generation.rds")
  frangieh_atomic_rds("generation-first", pointer)
  stopifnot(readRDS(file.path(frangieh_intermediate_directory(bundle), "sentinel.rds")) == "unchanged")
  saveRDS("../unsafe", pointer)
  expect_error(frangieh_intermediate_directory(bundle))
  provenance <- object$provenance
  provenance$inputs$md5[[1]] <- paste(rep("0", 32), collapse = "")
  frangieh_validate_provenance(provenance, script_dir, inputs)

  # Existing intermediates/replicates are refused before loading analysis
  # packages or fitting models, and the dispatcher never skips an old result.
  previous <- Sys.getenv(c("FRANGIEH_INTERMEDIATE_DIR", "FRANGIEH_SCORE_RESULTS_DIR"), unset = NA_character_)
  on.exit(for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(previous[[name]]), name))
  }, add = TRUE)
  Sys.setenv(FRANGIEH_INTERMEDIATE_DIR = inputs, FRANGIEH_SCORE_RESULTS_DIR = sandbox)
  expect_cli_refusal <- function(executable, args, message) {
    output <- suppressWarnings(system2(executable, args, stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    stopifnot(!is.null(status), status != 0L, any(grepl(message, output, fixed = TRUE)))
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  expect_cli_refusal(rscript, c("--vanilla", shQuote(file.path(script_dir, "run_combined_fdr.R")), "1"), "Refusing to overwrite")
  expect_cli_refusal(rscript, c("--vanilla", shQuote(file.path(script_dir, "data_preprocessing.R")), "1000"), "fresh empty intermediate directory")
  expect_cli_refusal("bash", c(shQuote(file.path(script_dir, "submit_combined_fdr.sh")), "1", "--local"), "fresh empty directory")
  dispatch <- file.path(sandbox, "previous-dispatch")
  Sys.setenv(FRANGIEH_SCORE_RESULTS_DIR = dispatch)
  # Invalid synthetic intermediates stop the first dispatch without fitting;
  # its permanent claim must still prevent a second attempt with no results.
  attempt <- suppressWarnings(system2("bash", c(shQuote(file.path(script_dir, "submit_combined_fdr.sh")), "1", "--local"), stdout = TRUE, stderr = TRUE))
  stopifnot(!is.null(attr(attempt, "status")), attr(attempt, "status") != 0L,
            dir.exists(file.path(dispatch, ".replicates-dispatched")))
  expect_cli_refusal("bash", c(shQuote(file.path(script_dir, "submit_combined_fdr.sh")), "1", "--local"), "fresh empty directory")
  stopifnot(identical(tools::md5sum(files), hashes));cat("Single-q10 fresh-run refusals, collections and atomic publication passed; synthetic fixtures unchanged.\n")
}
run_io_checks()
