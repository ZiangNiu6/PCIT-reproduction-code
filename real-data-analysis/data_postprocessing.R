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
source(file.path(script_dir, "..", "reproduction-helpers.R"))
project_root <- reproduction_root(script_dir)
results_root <- reproduction_results_root(project_root)
source(file.path(script_dir, "method_config.R"))
results_dir <- Sys.getenv(
  "FRANGIEH_SCORE_RESULTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_RESULTS_SUBDIR)
)

source(file.path(script_dir, "workflow-io.R"))
source(file.path(script_dir, "result-validation.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) && !identical(args, "--validate-only"))) {
  stop("Usage: Rscript data_postprocessing.R [--validate-only]")
}
rep_files <- frangieh_result_files(results_dir)
cat(sprintf("Reading and validating %d replicate files...\n", length(rep_files)))
rep_objects <- frangieh_validate_collection(rep_files)
reference_config <- rep_objects[[1L]]$config
if (identical(args, "--validate-only")) {
  cat("All 20 Frangieh v2 replicate outputs are valid; no files written.\n")
  quit(save = "no", status = 0L)
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
frangieh_atomic_rds(pair_long, file.path(results_dir, "score_pair_long.rds"), replace = TRUE)

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
frangieh_atomic_rds(summary_df, file.path(results_dir, "score_fdr_summary.rds"), replace = TRUE)

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
frangieh_atomic_rds(aggregated_df, file.path(results_dir, "score_fdr_aggregated.rds"), replace = TRUE)

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
frangieh_atomic_rds(diagnostics_df, file.path(results_dir, "score_diagnostics.rds"), replace = TRUE)

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
