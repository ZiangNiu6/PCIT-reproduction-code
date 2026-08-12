#!/usr/bin/env Rscript
## Aggregate the five-method Frangieh score-calibration experiment.
##
## Input:
##   results/score-calibration-nonpositive-nulls/score_fdr_repNNN.rds
##
## Outputs in results/score-calibration-nonpositive-nulls/:
##   score_pair_long.rds       pair x method x q p-values and decisions
##   score_fdr_summary.rds     replicate-level FDR and power
##   score_fdr_aggregated.rds  across-replicate means and standard errors
##   score_diagnostics.rds     adaptive-run timing and stopping diagnostics
##
## Rejection sets are always the sets stored by the replicate runner.  The
## adaptive methods' p-values are q-specific, so the long output includes q
## rather than choosing one p-value vector after the fact.

suppressPackageStartupMessages({
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
results_root <- Sys.getenv(
  "LOCAL_PDRT_DATA_DIR",
  unset = file.path(project_root, "results")
)
source(file.path(script_dir, "method_config.R"))
results_dir <- Sys.getenv(
  "FRANGIEH_SCORE_RESULTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_RESULTS_SUBDIR)
)

rep_files <- list.files(
  results_dir,
  pattern = "^score_fdr_rep[0-9]{3}\\.rds$",
  full.names = TRUE
)
if (length(rep_files) == 0L) {
  stop("No score_fdr_repNNN.rds files found in ", results_dir)
}
rep_files <- sort(rep_files)

schema_error <- function(path, ...) {
  stop(basename(path), ": incompatible ", PERTURB_SCHEMA_VERSION,
       " result: ", paste0(..., collapse = ""), call. = FALSE)
}

check_rejection_indices <- function(x, n_pairs, path, label) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != as.integer(x)) || any(x < 1L | x > n_pairs) ||
      anyDuplicated(x)) {
    schema_error(path, label, " must contain unique pair indices in 1:",
                 n_pairs)
  }
  invisible(TRUE)
}

check_p_values <- function(x, n_pairs, path, label) {
  if (!is.numeric(x) || length(x) != n_pairs || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0 | x > 1)) {
    schema_error(path, label, " must be a finite numeric vector of length ",
                 n_pairs, " with values in [0, 1]")
  }
  invisible(TRUE)
}

validate_rep <- function(obj, path) {
  if (!is.list(obj) || !all(c("config", "result") %in% names(obj))) {
    schema_error(path, "top level must contain config and result")
  }
  cfg <- obj$config
  res <- obj$result
  if (!identical(cfg$schema_version, PERTURB_SCHEMA_VERSION)) {
    schema_error(
      path, "schema_version is ", deparse(cfg$schema_version), "; expected ",
      deparse(PERTURB_SCHEMA_VERSION),
      ". Do not mix legacy and corrected result files"
    )
  }
  if (!identical(cfg$method_levels, PERTURB_METHOD_LEVELS)) {
    schema_error(path, "config$method_levels does not match method_config.R")
  }
  required_config <- c(
    "n_reps_max", "rep_idx", "seed", "n_pos", "n_null", "n_total",
    "median_size", "n_cells_full", "n_pure_nt", "q_levels", "h",
    "B_max", "alternative", "null_gene_policy", "null_gene_pool_size",
    "null_gene_pool"
  )
  if (!all(required_config %in% names(cfg))) {
    schema_error(path, "config is missing: ",
                 paste(setdiff(required_config, names(cfg)), collapse = ", "))
  }
  if (length(cfg$h) != 1L || !is.numeric(cfg$h) || cfg$h < 1L ||
      length(cfg$B_max) != 1L || !is.numeric(cfg$B_max) || cfg$B_max < 1L ||
      !identical(cfg$alternative, "two.sided")) {
    schema_error(path, "config has invalid h, B_max, or alternative")
  }
  if (!identical(cfg$null_gene_policy, PERTURB_NULL_GENE_POLICY) ||
      !is.character(cfg$null_gene_pool) || anyNA(cfg$null_gene_pool) ||
      any(cfg$null_gene_pool == "") || anyDuplicated(cfg$null_gene_pool) ||
      !identical(
        as.integer(cfg$null_gene_pool_size),
        as.integer(length(cfg$null_gene_pool))
      )) {
    schema_error(path, "config has an invalid non-positive null-gene pool")
  }

  q_levels <- as.numeric(cfg$q_levels)
  if (length(q_levels) == 0L || anyNA(q_levels) ||
      any(q_levels <= 0 | q_levels >= 1) || anyDuplicated(q_levels)) {
    schema_error(path, "config$q_levels must be unique probabilities in (0, 1)")
  }
  q_keys <- perturb_q_key(q_levels)

  required_result <- c("rep", "seed", "pairs", "methods", "timings")
  if (!all(required_result %in% names(res))) {
    schema_error(path, "result is missing: ",
                 paste(setdiff(required_result, names(res)), collapse = ", "))
  }
  if (length(res$rep) != 1L || !is.numeric(res$rep) || is.na(res$rep)) {
    schema_error(path, "result$rep must be one non-missing number")
  }
  if (!identical(as.integer(res$rep), as.integer(cfg$rep_idx))) {
    schema_error(path, "result$rep and config$rep_idx disagree")
  }
  if (length(res$seed) != 1L || !is.numeric(res$seed) || is.na(res$seed) ||
      !identical(as.integer(res$seed), as.integer(cfg$seed))) {
    schema_error(path, "result$seed and config$seed must agree")
  }

  required_pair_columns <- c(
    "pair_idx", "gene", "pair_type", "is_alt", "nt_subset", "n_cells",
    "n_exposed", "theta_hat", "score_sum", "information_denom_sq",
    "empirical_denom_sq", "information_stat", "empirical_stat",
    "score_eff_n", "max_score_share"
  )
  pairs <- res$pairs
  if (!is.data.frame(pairs) ||
      !all(required_pair_columns %in% names(pairs))) {
    schema_error(path, "result$pairs is missing: ",
                 paste(setdiff(required_pair_columns, names(pairs)),
                       collapse = ", "))
  }
  n_pairs <- nrow(pairs)
  if (n_pairs == 0L || !identical(as.integer(pairs$pair_idx), seq_len(n_pairs))) {
    schema_error(path, "pairs$pair_idx must equal 1:nrow(pairs)")
  }
  if (!is.logical(pairs$is_alt) || anyNA(pairs$is_alt)) {
    schema_error(path, "pairs$is_alt must be a non-missing logical vector")
  }
  if (anyNA(pairs$pair_type) ||
      !all(as.character(pairs$pair_type) == ifelse(pairs$is_alt, "alt", "null"))) {
    schema_error(path, "pairs$pair_type and pairs$is_alt disagree")
  }
  if (!is.null(cfg$n_total) && as.integer(cfg$n_total) != n_pairs) {
    schema_error(path, "config$n_total and nrow(result$pairs) disagree")
  }
  if (sum(pairs$is_alt) != as.integer(cfg$n_pos) ||
      sum(!pairs$is_alt) != as.integer(cfg$n_null)) {
    schema_error(path, "configured alternative/null counts disagree with pairs")
  }
  alt_genes <- unique(as.character(pairs$gene[pairs$is_alt]))
  null_genes <- as.character(pairs$gene[!pairs$is_alt])
  if (anyNA(null_genes) || any(null_genes == "") ||
      any(null_genes %in% alt_genes) ||
      !all(null_genes %in% cfg$null_gene_pool)) {
    schema_error(
      path,
      "null pairs must use only configured non-positive response genes"
    )
  }

  methods <- res$methods
  if (!is.list(methods) || !identical(names(methods), PERTURB_METHOD_LEVELS)) {
    schema_error(path, "result$methods must have the configured methods in order")
  }
  for (method in PERTURB_METHOD_LEVELS) {
    method_result <- methods[[method]]
    if (!is.list(method_result) || !is.numeric(method_result$observed) ||
        length(method_result$observed) != n_pairs) {
      schema_error(path, "methods$", method,
                   "$observed must be numeric with one value per pair")
    }
    if (method == "parametric_score") {
      check_p_values(method_result$p_value, n_pairs, path,
                     "methods$parametric_score$p_value")
      if (!is.list(method_result$rejected) ||
          !setequal(names(method_result$rejected), q_keys)) {
        schema_error(path, "methods$parametric_score$rejected must be keyed by q")
      }
      for (key in q_keys) {
        check_rejection_indices(
          method_result$rejected[[key]], n_pairs, path,
          paste0("methods$parametric_score$rejected[[", key, "]]")
        )
      }
    } else {
      if (!is.list(method_result$by_q) ||
          !setequal(names(method_result$by_q), q_keys)) {
        schema_error(path, "methods$", method, "$by_q must be keyed by q")
      }
      for (key in q_keys) {
        by_q <- method_result$by_q[[key]]
        required_by_q <- c(
          "p_value", "rejected", "futility", "active", "K",
          "stop_round", "t_total", "R_threshold"
        )
        if (!is.list(by_q) || !all(required_by_q %in% names(by_q))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]] is missing adaptive-result fields")
        }
        check_p_values(by_q$p_value, n_pairs, path,
                       paste0("methods$", method, "$by_q[[", key,
                              "]]$p_value"))
        check_rejection_indices(
          by_q$rejected, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$rejected")
        )
        check_rejection_indices(
          by_q$futility, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$futility")
        )
        check_rejection_indices(
          by_q$active, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$active")
        )
        partition <- c(by_q$rejected, by_q$futility, by_q$active)
        if (!identical(sort(as.integer(partition)), seq_len(n_pairs))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]] rejection/futility/active sets do not partition pairs")
        }
        if (!is.numeric(by_q$K) || length(by_q$K) != n_pairs ||
            anyNA(by_q$K) || any(by_q$K != as.integer(by_q$K)) ||
            any(by_q$K < 0L | by_q$K > as.integer(cfg$h))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$K must contain one loss count in 0:h per pair")
        }
        if (!is.numeric(by_q$stop_round) || length(by_q$stop_round) != n_pairs) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$stop_round must have one value per pair")
        }
        if (!all(is.na(by_q$stop_round[by_q$active])) ||
            anyNA(by_q$stop_round[c(by_q$rejected, by_q$futility)])) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$stop_round must be NA exactly for active pairs")
        }
      }
    }
  }

  timings <- res$timings
  if (!is.list(timings) || !is.numeric(timings$cache) ||
      length(timings$cache) != 1L || !is.finite(timings$cache)) {
    schema_error(path, "result$timings$cache must be one finite number")
  }
  for (method in PERTURB_RESAMPLING_METHODS) {
    method_timings <- timings[[method]]
    if (!is.numeric(method_timings) ||
        !setequal(names(method_timings), q_keys) ||
        anyNA(method_timings) || any(!is.finite(method_timings))) {
      schema_error(path, "result$timings$", method,
                   " must be a finite numeric vector keyed by q")
    }
  }

  obj
}

cat(sprintf("Reading and validating %d replicate files...\n", length(rep_files)))
rep_objects <- Map(
  function(path) validate_rep(readRDS(path), path),
  rep_files
)
# unname(): Map() names rep_objects by rep_files, and vapply() propagates those
# names, so without this sort(rep_ids) carries a names attribute and
# identical(..., seq_len(20L)) fails even when the values are exactly 1:20.
rep_ids <- unname(vapply(rep_objects, function(x) as.integer(x$result$rep), integer(1)))
if (anyDuplicated(rep_ids)) stop("Duplicate replicate IDs in score result files")
expected_rep_ids <- seq_len(20L)
if (!identical(sort(rep_ids), expected_rep_ids)) {
  stop(
    "Canonical manuscript outputs require exactly replicate IDs 1:20; found: ",
    paste(sort(rep_ids), collapse = ", "),
    call. = FALSE
  )
}

analysis_config_fields <- c(
  "schema_version", "method_levels", "n_reps_max", "n_pos", "n_null",
  "n_total", "median_size", "n_cells_full", "n_pure_nt", "q_levels",
  "h", "B_max", "alternative", "null_gene_policy",
  "null_gene_pool_size", "null_gene_pool"
)
reference_config <- rep_objects[[1L]]$config[analysis_config_fields]
for (i in seq_along(rep_objects)[-1L]) {
  candidate_config <- rep_objects[[i]]$config[analysis_config_fields]
  if (!identical(candidate_config, reference_config)) {
    differing <- analysis_config_fields[vapply(
      analysis_config_fields,
      function(field) {
        !identical(candidate_config[[field]], reference_config[[field]])
      },
      logical(1)
    )]
    stop(
      basename(rep_files[[i]]),
      ": analysis configuration differs from the first replicate: ",
      paste(differing, collapse = ", "),
      call. = FALSE
    )
  }
}
q_levels <- as.numeric(reference_config$q_levels)
q_keys <- perturb_q_key(q_levels)

get_q_result <- function(res, method, key) {
  if (method == "parametric_score") {
    list(
      p_value = res$methods[[method]]$p_value,
      rejected = res$methods[[method]]$rejected[[key]]
    )
  } else {
    res$methods[[method]]$by_q[[key]]
  }
}

## Pair-level long table -----------------------------------------------------
build_pair_long <- function(obj) {
  res <- obj$result
  n_pairs <- nrow(res$pairs)
  rows <- vector(
    "list", length(PERTURB_METHOD_LEVELS) * length(q_levels)
  )
  row_idx <- 0L
  for (q_idx in seq_along(q_levels)) {
    q <- q_levels[[q_idx]]
    key <- q_keys[[q_idx]]
    for (method in PERTURB_METHOD_LEVELS) {
      row_idx <- row_idx + 1L
      q_result <- get_q_result(res, method, key)
      rows[[row_idx]] <- data.frame(
        rep = as.integer(res$rep),
        seed = as.integer(res$seed),
        res$pairs,
        method = method,
        q = q,
        observed = res$methods[[method]]$observed,
        p_value = q_result$p_value,
        rejected = seq_len(n_pairs) %in% q_result$rejected,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  bind_rows(rows)
}

cat("Building q-specific pair-level table...\n")
pair_long <- bind_rows(lapply(rep_objects, build_pair_long))
pair_long$method <- factor(pair_long$method, levels = PERTURB_METHOD_LEVELS)
pair_long <- pair_long |>
  mutate(schema_version = PERTURB_SCHEMA_VERSION, .before = 1L)
saveRDS(pair_long, file.path(results_dir, "score_pair_long.rds"))

## Replicate-level FDR and power --------------------------------------------
evaluate_rejections <- function(rejected, is_alt) {
  n_rejected <- length(rejected)
  n_true_pos <- sum(is_alt[rejected])
  n_false_pos <- n_rejected - n_true_pos
  data.frame(
    n_rejected = n_rejected,
    n_true_pos = n_true_pos,
    n_false_pos = n_false_pos,
    power = if (sum(is_alt) == 0L) NA_real_ else n_true_pos / sum(is_alt),
    fdr = if (n_rejected == 0L) 0 else n_false_pos / n_rejected
  )
}

build_summary <- function(obj) {
  res <- obj$result
  rows <- vector(
    "list", length(PERTURB_METHOD_LEVELS) * length(q_levels)
  )
  row_idx <- 0L
  for (q_idx in seq_along(q_levels)) {
    q <- q_levels[[q_idx]]
    key <- q_keys[[q_idx]]
    for (method in PERTURB_METHOD_LEVELS) {
      row_idx <- row_idx + 1L
      q_result <- get_q_result(res, method, key)
      rows[[row_idx]] <- bind_cols(
        data.frame(rep = as.integer(res$rep), method = method, q = q),
        evaluate_rejections(q_result$rejected, res$pairs$is_alt)
      )
    }
  }
  bind_rows(rows)
}

cat("Computing replicate-level FDR and power from stored rejection sets...\n")
summary_df <- bind_rows(lapply(rep_objects, build_summary))
summary_df$method <- factor(summary_df$method, levels = PERTURB_METHOD_LEVELS)
summary_df <- summary_df |>
  mutate(schema_version = PERTURB_SCHEMA_VERSION, .before = 1L)
saveRDS(summary_df, file.path(results_dir, "score_fdr_summary.rds"))

se_mean <- function(x) {
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

aggregated_df <- summary_df |>
  group_by(schema_version, method, q) |>
  summarise(
    mean_fdr = mean(fdr),
    se_fdr = se_mean(fdr),
    mean_power = mean(power),
    se_power = se_mean(power),
    mean_rejected = mean(n_rejected),
    se_rejected = se_mean(n_rejected),
    n_reps = n(),
    .groups = "drop"
  )
saveRDS(aggregated_df, file.path(results_dir, "score_fdr_aggregated.rds"))

## Adaptive-method diagnostics ---------------------------------------------
build_diagnostics <- function(obj) {
  res <- obj$result
  rows <- vector(
    "list", length(PERTURB_RESAMPLING_METHODS) * length(q_levels)
  )
  row_idx <- 0L
  for (q_idx in seq_along(q_levels)) {
    q <- q_levels[[q_idx]]
    key <- q_keys[[q_idx]]
    for (method in PERTURB_RESAMPLING_METHODS) {
      row_idx <- row_idx + 1L
      adaptive <- res$methods[[method]]$by_q[[key]]
      stop_round <- as.numeric(adaptive$stop_round)
      completed_stop_round <- stop_round[!is.na(stop_round)]
      rows[[row_idx]] <- data.frame(
        rep = as.integer(res$rep),
        q = q,
        method = method,
        cache_secs = as.numeric(res$timings$cache),
        method_secs = as.numeric(res$timings[[method]][[key]]),
        t_total = as.numeric(adaptive$t_total),
        n_rejected = length(adaptive$rejected),
        n_futility = length(adaptive$futility),
        n_active = length(adaptive$active),
        mean_stop_round = if (length(completed_stop_round) == 0L) {
          NA_real_
        } else {
          mean(completed_stop_round)
        },
        max_stop_round = if (length(completed_stop_round) == 0L) {
          NA_real_
        } else {
          max(completed_stop_round)
        },
        stringsAsFactors = FALSE
      ) |>
        mutate(R_threshold = list(adaptive$R_threshold))
    }
  }
  bind_rows(rows)
}

cat("Collecting adaptive-run diagnostics...\n")
diagnostics_df <- bind_rows(lapply(rep_objects, build_diagnostics))
diagnostics_df$method <- factor(
  diagnostics_df$method, levels = PERTURB_RESAMPLING_METHODS
)
diagnostics_df <- diagnostics_df |>
  mutate(schema_version = PERTURB_SCHEMA_VERSION, .before = 1L)
saveRDS(diagnostics_df, file.path(results_dir, "score_diagnostics.rds"))

cat("\nFDR and power summary (mean +/- SE):\n")
print(
  aggregated_df |>
    mutate(
      method = PERTURB_METHOD_DISPLAY[as.character(method)],
      fdr = sprintf("%.3f +/- %.3f", mean_fdr, se_fdr),
      power = sprintf("%.3f +/- %.3f", mean_power, se_power)
    ) |>
    select(method, q, mean_rejected, fdr, power)
)

cat("\nWrote:\n")
for (filename in c(
  "score_pair_long.rds", "score_fdr_summary.rds",
  "score_fdr_aggregated.rds", "score_diagnostics.rds"
)) {
  cat("  ", file.path(results_dir, filename), "\n", sep = "")
}
