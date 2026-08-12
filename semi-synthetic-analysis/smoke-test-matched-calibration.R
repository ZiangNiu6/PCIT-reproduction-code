#!/usr/bin/env Rscript

# Fast, data-independent checks for the CCLE matched-calibration machinery.
# This intentionally uses tiny random forests and a tiny capped adaptive
# permutation budget; it validates code invariants rather than performance.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not determine this script's path", call. = FALSE)
}
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "matched-calibration-utils.R"))

required_packages <- c("pcmt", "comets", "ranger")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Smoke test requires package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

seeds <- ccle_seed_schedule(setup_idx = 2L, rep_idx = 3L)
stopifnot(
  identical(names(seeds), names(CCLE_SEED_COMPONENT_OFFSETS)),
  length(unique(unname(seeds))) == length(seeds)
)

n <- 48L
synthetic <- ccle_with_seed(seeds[["data_seed"]], {
  Z <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
  X <- 0.7 * Z[, 1L] - 0.3 * Z[, 2L] + stats::rnorm(n)
  Y <- Z[, 1L]^2 + Z[, 3L] + stats::rnorm(n)
  list(Y = Y, X = X, Z = Z)
})
split_D2 <- ccle_make_pcm_split(
  n = n,
  split_seed = seeds[["pcm_split_seed"]],
  frac = 0.5
)
stopifnot(
  identical(
    split_D2,
    ccle_make_pcm_split(n, seeds[["pcm_split_seed"]], frac = 0.5)
  ),
  length(split_D2) == ceiling(n / 2)
)

run_tests <- function(order) {
  ccle_run_matched_calibration(
    Y = synthetic$Y,
    X = synthetic$X,
    Z = synthetic$Z,
    seeds = seeds,
    pcm_D2_indices = split_D2,
    B_max = 7L,
    adaptive_h = 3L,
    alpha = 0.05,
    n_cores = 1L,
    alternative = "two.sided",
    pcm_frac = 0.5,
    rf_num_trees = 20L,
    test_order = order
  )
}

set.seed(991L)
forward_1 <- run_tests(c("gcm", "pcm"))
invisible(stats::runif(25L))
forward_2 <- run_tests(c("gcm", "pcm"))
reverse_order <- run_tests(c("pcm", "gcm"))

for (result in list(forward_1, forward_2, reverse_order)) {
  stopifnot(
    isTRUE(result$matched_calibration),
    identical(result$alternative, "two.sided"),
    result$diagnostics$gcm$actual_B >= 1L,
    result$diagnostics$gcm$actual_B <= 7L,
    result$diagnostics$pcm$actual_B >= 1L,
    result$diagnostics$pcm$actual_B <= 7L,
    result$diagnostics$gcm$stop_reason %in% c("rejection", "futility", "cap"),
    result$diagnostics$pcm$stop_reason %in% c("rejection", "futility", "cap"),
    result$diagnostics$gcm$loss_count >= 0L,
    result$diagnostics$gcm$loss_count <= 3L,
    result$diagnostics$pcm$loss_count >= 0L,
    result$diagnostics$pcm$loss_count <= 3L,
    result$diagnostics$gcm$adaptive_backend %in%
      c("pcmt_native", "ordered_perm_stats_fallback"),
    result$diagnostics$pcm$adaptive_backend %in%
      c("pcmt_native", "ordered_perm_stats_fallback"),
    is.na(result$diagnostics$gcm$error),
    is.na(result$diagnostics$pcm$error),
    identical(
      result$diagnostics$gcm$normal_statistic,
      result$diagnostics$gcm$permutation_statistic
    ),
    identical(
      result$diagnostics$pcm$normal_statistic,
      result$diagnostics$pcm$permutation_statistic
    )
  )

  expected_asym <- c(
    asym_gcm = 2 * stats::pnorm(
      abs(result$statistics[["gcm"]]),
      lower.tail = FALSE
    ),
    asym_pcm = 2 * stats::pnorm(
      abs(result$statistics[["pcm"]]),
      lower.tail = FALSE
    )
  )
  stopifnot(isTRUE(all.equal(
    unname(result$pvalues[names(expected_asym)]),
    unname(expected_asym),
    tolerance = 0
  )))
}

# Check the anytime-valid formula and stopping metadata independently of the
# fitted learners.
adaptive_check <- ccle_adaptive_calibration(
  statistic = 2,
  perm_stats = c(0, 0, 3, 0, 4, 0),
  h = 2L,
  alpha = 0.05,
  alternative = "two.sided"
)
stopifnot(
  identical(adaptive_check$adaptive$stop_reason, "futility"),
  identical(nrow(adaptive_check$adaptive$trajectory), 5L),
  identical(adaptive_check$adaptive$K, 2L),
  identical(adaptive_check$adaptive$t_star, 5L),
  identical(adaptive_check$p.value, 2 / 5)
)

# Component-local RNG streams make the result repeatable despite unrelated RNG
# use and invariant to whether GCM or PCM is executed first.
stopifnot(
  identical(forward_1$pvalues, forward_2$pvalues),
  identical(forward_1$statistics, forward_2$statistics),
  identical(forward_1$pvalues, reverse_order$pvalues),
  identical(forward_1$statistics, reverse_order$statistics)
)

cat("Matched-calibration smoke test passed.\n")
