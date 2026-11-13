#!/usr/bin/env Rscript
## Data-independent q10 schema and self-contained execution-adoption tests.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
for (name in c("method_config.R", "workflow-io.R", "result-validation.R")) source(file.path(script_dir, name))
expect_error <- function(expression) {
  if (!inherits(tryCatch(force(expression), error = identity), "error")) stop("Expected a fail-closed contract error")
}
make_frangieh_fixture <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  contract <- list(schema = "frangieh_gb2024_preprocessing_v1", policy = PERTURB_COHORT_POLICY,
    feature_stat_policy = "historical_stored_counts", cell_rule = "exactly_one_original_assigned_guide_before_feature_qc",
    reference_commit = "005c28b9ba30e7d39c3ee18d51111722f49d32ec", ondisc_commit = "5adcc537664aada08c3e79f08f2a944a0d31598d",
    gene_fraction_min = .005, guide_cells_min = 10L, covariates = c("log_n_umis", "log_n_nonzero"), id_cleanup = "underscore_to_dash",
    counts = list(input_cells = 100L, retained_cells = 80L, pure_nt_cells = 20L, input_genes = 2000L,
      retained_genes = 1500L, input_guides = 210L, eligible_guides = 200L, eligible_target_groups = 50L))
  for (file in frangieh_intermediate_files) saveRDS("synthetic fixture", file.path(directory, file))
  saveRDS(contract, file.path(directory, "preprocessing-contract.rds"))
  config <- perturb_run_config(1L, contract, 20L, 3L, paste0("synthetic-null-", 1:950))
  provenance <- list(schema = PERTURB_PROVENANCE_SCHEMA, kind = "native_run", role = "native_execution",
    inputs = frangieh_file_manifest(file.path(directory, frangieh_intermediate_files)),
    code = frangieh_file_manifest(frangieh_code_files(script_dir)),
    packages = setNames(rep("synthetic-fixture", 5L), c("pcmt", "ondisc", "MASS", "dplyr", "Matrix")),
    R = as.character(getRversion()), cohort_policy = PERTURB_COHORT_POLICY, preprocessing_contract = contract)
  pairs <- data.frame(pair_idx = 1:1000, gene = c(paste0("synthetic-alt-", 1:50), config$null_gene_pool),
    pair_type = c(rep("alt", 50), rep("null", 950)), is_alt = c(rep(TRUE, 50), rep(FALSE, 950)),
    nt_subset = NA_character_, n_cells = c(rep(21L, 50), rep(20L, 950)), n_exposed = 1L,
    theta_hat = 1, score_sum = 0, information_denom_sq = 1, empirical_denom_sq = 1,
    information_stat = 0, empirical_stat = 0, score_eff_n = 1, max_score_share = 1)
  key <- perturb_q_key(PERTURB_Q_LEVELS)
  methods <- setNames(lapply(PERTURB_METHOD_LEVELS, function(method) {
    if (method == "parametric_score") return(list(observed = rep(0, 1000), p_value = rep(1, 1000), rejected = setNames(list(integer()), key)))
    adaptive <- list(p_value = rep(1, 1000), rejected = integer(), futility = 1:1000, active = integer(),
      K = rep(15L, 1000), stop_round = rep(15L, 1000), t_total = 15L, R_threshold = 0)
    list(observed = rep(0, 1000), by_q = setNames(list(adaptive), key))
  }), PERTURB_METHOD_LEVELS)
  timings <- c(list(cache = 0), setNames(rep(list(setNames(0, key)), 4), PERTURB_RESAMPLING_METHODS))
  list(config = config, result = list(rep = 1L, seed = 1001L, pairs = pairs, methods = methods, timings = timings), provenance = provenance)
}
make_adopted_fixture <- function(object, directory) {
  cfg <- list(schema_version = "frangieh_score5_q10_extension_v1", pilot = FALSE, q = .1,
    h = 15L, B_max = 25000L, method_levels = PERTURB_METHOD_LEVELS, RNGkind = PERTURB_RNG_KIND,
    seed_reference_B_max = 25000L, seed_reference_q = .3, seed_reference_q_levels = c(.2, .3))
  nms <- c("run_q10_rep.R", "score_helpers.R", "method_config.R")
  snapshots <- setNames(lapply(nms, function(name) charToRaw(paste("Synthetic source fixture:", name))), nms)
  sources <- lapply(nms, function(name) list(path = paste0("/nonexistent-execution/", name), bytes = length(snapshots[[name]]),
    sha256 = digest::digest(snapshots[[name]], algo = "sha256", serialize = FALSE)))
  inputs <- list(list(path = "/nonexistent-execution/input", bytes = 1L, sha256 = paste(rep("0", 64), collapse = "")))
  fp <- digest::digest(list(config = cfg, inputs = inputs, sources = sources), algo = "sha256")
  result <- object$result;result$timings <- setNames(as.list(rep(0, 4)), PERTURB_RESAMPLING_METHODS)
  original <- list(config = cfg, result = result, provenance = list(inputs = inputs, sources = sources, fingerprint = fp))
  file <- file.path(directory, "synthetic-selected.rds");saveRDS(original, file)
  bytes <- readBin(file, "raw", n = file.info(file)$size);unlink(file)
  record <- list(path = "/nonexistent-execution/selected.rds", bytes = length(bytes), sha256 = digest::digest(bytes, algo = "sha256", serialize = FALSE))
  receipt <- list(status = "complete", fingerprint = fp, config = cfg, inputs = inputs, sources = sources, outputs = list(record))
  frangieh_adopt_verified_q10(original, record, receipt, list(original_rds = bytes, source_snapshots = snapshots), object$provenance, object$config)
}
run_contract_checks <- function() {
  sandbox <- tempfile("frangieh-q10-contract-");dir.create(sandbox);on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  object <- make_frangieh_fixture(file.path(sandbox, "inputs"));path <- file.path(sandbox, "score_fdr_rep001.rds")
  invisible(validate_rep(object, path))
  stopifnot(identical(PERTURB_SCHEMA_VERSION, "perturb_score5_v4_gb2024_q10"), identical(PERTURB_Q_LEVELS, .1))
  for (field in c("schema_version", "rng_layout", "preprocessing_policy")) { bad <- object;bad$config[[field]] <- "incompatible";expect_error(validate_rep(bad, path)) }
  bad <- object;bad$config$q_levels <- .05;expect_error(validate_rep(bad, path))
  bad <- object;bad$config$B_max <- 19L;expect_error(validate_rep(bad, path))
  bad <- object;bad$config$covariate_names <- c("log_n_umis", "batch");expect_error(validate_rep(bad, path))
  bad <- object;bad$result$pairs$n_cells[[51]] <- 21L;expect_error(validate_rep(bad, path))
  bad <- object;bad$result$methods[[2]]$by_q[[1]]$p_value[[1]] <- .5;expect_error(validate_rep(bad, path))
  bad <- object;bad$provenance <- NULL;expect_error(validate_rep(bad, path))
  # Recorded source/input fingerprints may differ from today's files. They
  # describe the run and do not turn saved-result validation into a rerun gate.
  historical <- object
  historical$provenance$code$md5[[1]] <- paste(rep("0", 32), collapse = "")
  historical$provenance$inputs$md5[[1]] <- paste(rep("0", 32), collapse = "")
  invisible(validate_rep(historical, path, intermediate_dir = file.path(sandbox, "inputs")))
  adopted <- make_adopted_fixture(object, sandbox);invisible(validate_rep(adopted, path))
  stopifnot(identical(adopted$provenance$kind, "imported_verified_q10"), is.null(adopted$result$timings$cache))
  bad <- adopted;bad$result$pairs$score_sum[[1]] <- 1;expect_error(validate_rep(bad, path))
  bad <- adopted;bad$provenance$adoption$archive$original_rds[[1]] <- as.raw(0);expect_error(validate_rep(bad, path))
  bad <- adopted;bad$provenance$adoption$archive$source_snapshots[[1]][[1]] <- as.raw(0);expect_error(validate_rep(bad, path))
  bad <- adopted;bad$provenance$kind <- "native_run";bad$provenance$role <- "native_execution";expect_error(validate_rep(bad, path))
  cat("Single-q10 scientific/result and origin-archive contracts passed; execution fingerprints are metadata.\n")
}
if (sys.nframe() == 0L) run_contract_checks()
