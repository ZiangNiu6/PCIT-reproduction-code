#!/usr/bin/env Rscript
## Focused read-only validation/atomic-publication regression checks. All test
## mutations are isolated in tempdir(); no numerical experiment is executed.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "reproduction-helpers.R"))
root <- reproduction_root(script_dir)
source(file.path(script_dir, "method_config.R"))
source(file.path(script_dir, "workflow-io.R"))
source(file.path(script_dir, "result-validation.R"))

expect_error <- function(expression) {
  failed <- inherits(tryCatch(force(expression), error = identity), "error")
  if (!failed) stop("Expected a fail-closed validation error")
}

run_checks <- function() {
  fixture <- file.path(root, "results", "frangieh-2021", PERTURB_RESULTS_SUBDIR)
  files <- frangieh_result_files(fixture)
  if (length(files) != 20L) stop("Tests require the 20 bundled v2 replicate fixtures")
  original_hashes <- tools::md5sum(files)
  sandbox <- tempfile("frangieh-io-test-")
  dir.create(sandbox)
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  stopifnot(all(file.copy(files, sandbox)))
  copied <- frangieh_result_files(sandbox)
  set.seed(123)
  rng <- .Random.seed
  stopifnot(length(frangieh_validate_collection(copied)) == 20L)
  stopifnot(identical(.Random.seed, rng))
  expect_error(frangieh_validate_collection(copied[-20L]))
  obj <- readRDS(copied[[1L]])
  invalid <- obj
  invalid$config$B_max <- 1L
  expect_error(validate_rep(invalid, copied[[1L]]))
  expect_error(validate_rep(obj, file.path(sandbox, "score_fdr_rep002.rds")))
  invalid <- obj
  invalid$config$schema_version <- "legacy"
  expect_error(validate_rep(invalid, copied[[1L]]))
  invalid <- obj
  invalid$result$methods[[2L]]$by_q[[1L]]$p_value[[1L]] <- NA_real_
  expect_error(validate_rep(invalid, copied[[1L]]))
  unexpected <- file.path(sandbox, "score_fdr_rep021.rds")
  file.copy(copied[[1L]], unexpected)
  expect_error(frangieh_result_files(sandbox))
  unlink(unexpected)

  atomic <- file.path(sandbox, "atomic.rds")
  frangieh_atomic_rds("old", atomic)
  expect_error(frangieh_atomic_rds("unauthorized replacement", atomic))
  stopifnot(identical(readRDS(atomic), "old"))
  frangieh_atomic_rds("new", atomic, replace = TRUE)
  stopifnot(identical(readRDS(atomic), "new"))
  publication_lock <- paste0(atomic, ".publish-lock")
  stopifnot(!dir.exists(publication_lock))
  dir.create(publication_lock)
  expect_error(frangieh_atomic_rds("concurrent replacement", atomic, replace = TRUE))
  stopifnot(identical(readRDS(atomic), "new"), dir.exists(publication_lock))
  unlink(publication_lock, recursive = TRUE)
  stopifnot(!length(list.files(sandbox, pattern = "^[.]atomic", all.files = TRUE)))

  bundle <- file.path(sandbox, "bundle")
  for (generation in c("generation-old", "generation-new")) {
    directory <- file.path(bundle, ".generations", generation)
    dir.create(directory, recursive = TRUE)
    saveRDS(generation, file.path(directory, "sentinel.rds"))
  }
  pointer <- file.path(bundle, "current-generation.rds")
  frangieh_atomic_rds("generation-old", pointer)
  old_reader <- frangieh_intermediate_directory(bundle)
  frangieh_atomic_rds("generation-new", pointer, replace = TRUE)
  stopifnot(identical(readRDS(file.path(old_reader, "sentinel.rds")), "generation-old"))
  stopifnot(identical(readRDS(file.path(frangieh_intermediate_directory(bundle), "sentinel.rds")), "generation-new"))
  frangieh_atomic_rds("../unsafe", pointer, replace = TRUE)
  expect_error(frangieh_intermediate_directory(bundle))

  inputs <- file.path(sandbox, "inputs")
  dir.create(inputs)
  for (filename in frangieh_intermediate_files) saveRDS("fixture", file.path(inputs, filename))
  provenance <- list(schema = "frangieh_io_provenance_v1",
                     inputs = frangieh_file_manifest(file.path(inputs, frangieh_intermediate_files)),
                     code = frangieh_file_manifest(frangieh_code_files(script_dir)),
                     packages = setNames(rep("fixture", 5L), c("pcmt", "ondisc", "MASS", "dplyr", "Matrix")),
                     R = as.character(getRversion()),
                     cohort_policy = "original_v2_cohort_retained_without_scientific_changes")
  frangieh_validate_provenance(provenance, script_dir)
  provenance$code$md5[[1L]] <- paste(rep("0", 32L), collapse = "")
  expect_error(frangieh_validate_provenance(provenance, script_dir))
  stopifnot(identical(tools::md5sum(files), original_hashes))
  cat("Frangieh workflow I/O tests passed; bundled fixtures unchanged.\n")
}
run_checks()
