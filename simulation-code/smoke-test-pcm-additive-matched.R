#!/usr/bin/env Rscript

# Lightweight contract test for the Section 5.2 simulation methods.  This does
# not reproduce the simulation; it checks that one GCM and one PCM statistic
# each supply both their normal and permutation p-values, that dCRT is a
# separate simulatr method with a valid Gaussian model-X reference, and that
# unrelated use of R's ambient random-number stream cannot change either
# result.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_argument))
suppressPackageStartupMessages(library(simulatr))

required_pcmt_formals <- c(
  "args_YonZ", "args_XonZ", "B", "n_cores", "h", "alpha", "alternative"
)
missing_pcmt_formals <- setdiff(
  required_pcmt_formals,
  names(formals(pcmt::perm_gcm))
)
if (length(missing_pcmt_formals) > 0L) {
  stop(
    "The active pcmt build is not the repository-pinned adaptive version. ",
    "Run this check through run_all_simulation.sh after restoring renv.",
    call. = FALSE
  )
}

specifier_file <- file.path(
  dirname(script_file),
  "sim-spec-formal",
  "pcm-additive.R"
)
source(specifier_file, local = FALSE)

expected_run_methods <- c("matched_gcm_pcm_rf", "dcrt_rf")
if (!identical(names(run_method_functions), expected_run_methods) ||
    !identical(names(simulatr_spec@run_method_functions), expected_run_methods)) {
  stop("Section 5.2 simulatr methods are missing or out of order.",
       call. = FALSE)
}

set.seed(7701L)
smoke_data <- generate_data(
  n = 40L,
  dgp = "setting1_xindep",
  gamma_signal = 0,
  B_perm = 3L
)
set.seed(7702L)
second_seed_data <- generate_data(
  n = 40L,
  dgp = "setting1_xindep",
  gamma_signal = 0,
  B_perm = 3L
)
if (!is.numeric(smoke_data$dcrt_resampling_seed) ||
    length(smoke_data$dcrt_resampling_seed) != 1L ||
    smoke_data$dcrt_resampling_seed < 1L ||
    identical(
      smoke_data$dcrt_resampling_seed,
      second_seed_data$dcrt_resampling_seed
    )) {
  stop("dCRT resampling seeds are not valid dataset-specific seeds.",
       call. = FALSE)
}

first_result <- matched_gcm_pcm_rf_f(smoke_data)
first_dcrt <- dcrt_rf_f(smoke_data)
invisible(stats::runif(25L))
second_result <- matched_gcm_pcm_rf_f(smoke_data)
second_dcrt <- dcrt_rf_f(smoke_data)

if (!identical(first_result, second_result)) {
  stop(
    "Matched calibration changed after unrelated ambient RNG use.",
    call. = FALSE
  )
}
if (!identical(first_dcrt, second_dcrt)) {
  stop(
    "dCRT changed after unrelated ambient RNG use.",
    call. = FALSE
  )
}
if (!identical(first_result$schema_version, "pcm_additive_matched") ||
    !isTRUE(first_result$matched[["all"]]) ||
    !all(first_result$flags)) {
  stop("Matched-calibration flags or schema are invalid.", call. = FALSE)
}

expected_normal <- c(
  asym_gcm_rf = 2 * stats::pnorm(
    abs(first_result$statistics[["gcm"]]),
    lower.tail = FALSE
  ),
  asym_pcm_rf = stats::pnorm(
    first_result$statistics[["pcm"]],
    lower.tail = FALSE
  )
)
if (!isTRUE(all.equal(
  unname(first_result$p_values[names(expected_normal)]),
  unname(expected_normal),
  tolerance = 1e-12,
  check.attributes = FALSE
))) {
  stop(
    "Normal p-values were not derived from the stored observed statistics.",
    call. = FALSE
  )
}

if (any(!is.finite(first_result$p_values)) ||
    any(first_result$p_values < 0 | first_result$p_values > 1) ||
    first_result$pcm_split$D1_n + first_result$pcm_split$D2_n != 40L) {
  stop("Smoke-test output failed its numerical contract.", call. = FALSE)
}

if (!identical(first_dcrt$schema_version, "pcm_additive_dcrt_rf_v1") ||
    !identical(first_dcrt$alternative, "two.sided") ||
    !identical(names(first_dcrt$p_values), "dcrt_rf") ||
    !isTRUE(all.equal(
      unname(first_dcrt$p.value),
      unname(first_dcrt$p_values[["dcrt_rf"]]),
      tolerance = 0
    )) ||
    any(!is.finite(first_dcrt$p_values)) ||
    any(first_dcrt$p_values < 0 | first_dcrt$p_values > 1) ||
    !all(first_dcrt$flags)) {
  stop("dCRT output failed its schema or probability contract.", call. = FALSE)
}

expected_dcrt_p <- (
  1 + first_dcrt$monte_carlo$n_exceed
) / (
  1 + first_dcrt$monte_carlo$M
)
if (first_dcrt$monte_carlo$M != smoke_data$B_perm ||
    !isTRUE(all.equal(first_dcrt$p.value, expected_dcrt_p, tolerance = 0)) ||
    !isTRUE(all.equal(
      first_dcrt$monte_carlo$minimum_p_value,
      1 / (smoke_data$B_perm + 1),
      tolerance = 0
    )) ||
    !identical(first_dcrt$model_x$family, "gaussian") ||
    !identical(first_dcrt$model_x$mean, "random_forest_prediction") ||
    !identical(first_dcrt$model_x$variance, "homoskedastic") ||
    !identical(
      first_dcrt$model_x$variance_estimator,
      "mean_squared_in_sample_residual"
    ) ||
    !is.finite(first_dcrt$model_x$sigma2_hat) ||
    first_dcrt$model_x$sigma2_hat <= 0 ||
    !is.finite(first_dcrt$nuisance$analytic_reference_sd) ||
    first_dcrt$nuisance$analytic_reference_sd <= 0 ||
    !isTRUE(all.equal(
      first_dcrt$nuisance$reference_standardized_statistic,
      first_dcrt$statistic / first_dcrt$nuisance$analytic_reference_sd,
      tolerance = 1e-12
    )) ||
    !is.finite(first_dcrt$monte_carlo$reference_mean) ||
    !is.finite(first_dcrt$monte_carlo$reference_sd) ||
    first_dcrt$monte_carlo$reference_sd < 0 ||
    !identical(first_dcrt$configuration$mean_learner, "random_forest") ||
    !identical(first_dcrt$configuration$fit_scope, "in_sample") ||
    !identical(first_dcrt$configuration$model_x_family, "gaussian") ||
    !identical(first_dcrt$configuration$model_x_variance, "homoskedastic") ||
    !identical(
      first_dcrt$configuration$variance_estimator,
      "mean_squared_in_sample_residual_divisor_n"
    ) ||
    first_dcrt$configuration$M != smoke_data$B_perm ||
    !identical(first_dcrt$configuration$refit_model_x_mean, FALSE) ||
    !identical(
      first_dcrt$configuration$two_sided_rule,
      "absolute_statistic_plus_one"
    )) {
  stop("dCRT model-X or Monte Carlo diagnostics are inconsistent.",
       call. = FALSE)
}

# dCRT deliberately refits the two nuisances as its own simulatr method, but
# its purpose-specific RF seeds match GCM.  The normalized observed statistic
# must therefore recover the GCM statistic exactly on the shared dataset.
if (!isTRUE(all.equal(
  first_dcrt$nuisance$studentized_statistic,
  unname(first_result$statistics[["gcm"]]),
  tolerance = 1e-12
))) {
  stop("dCRT and GCM did not share the intended nuisance fits.", call. = FALSE)
}

cat("Section 5.2 matched-calibration and dCRT smoke test passed.\n")
