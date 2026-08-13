#!/usr/bin/env Rscript

# Fast contract test for the Section 5.1 continuous-dispersion NB simulation.
# This test does not estimate rejection rates and never writes manuscript files.

suppressPackageStartupMessages({
  library(simulatr)
})

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", script_argument))
specifier_file <- file.path(
  dirname(script_file), "sim-spec-formal", "nb-permutation-signflip.R"
)
source(specifier_file, local = FALSE)

stopifnot(
  nrow(parameter_grid) == 30L,
  identical(
    unique(parameter_grid$dgp_schema),
    "continuous_nb_dispersion_estimated_theta"
  ),
  identical(
    parameter_grid$gamma_signal[
      parameter_grid$null_design == "MI" &
        parameter_grid$scenario_id == "severe_variance"
    ],
    -seq(0, 0.22, length.out = 5L)
  ),
  identical(
    parameter_grid$gamma_signal[
      parameter_grid$null_design == "CI" &
        parameter_grid$scenario_id == "severe_variance"
    ],
    -seq(0, 0.27, length.out = 5L)
  )
)

make_dataset <- function(null_design, scenario_id, seed) {
  row <- parameter_grid[
    parameter_grid$null_design == null_design &
      parameter_grid$scenario_id == scenario_id &
      !parameter_grid$ground_truth,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Expected exactly one null grid row for ", null_design, "/", scenario_id)
  }
  arguments <- as.list(row[1L, intersect(
    formalArgs(generate_data), names(row)
  ), drop = FALSE])
  arguments$B_resample <- 19L
  set.seed(seed)
  do.call(generate_data, arguments)
}

ci_data <- make_dataset("CI", "severe_variance", 20240813L)
mi_data <- make_dataset("MI", "severe_variance", 20240814L)
ordinary_data <- make_dataset("CI", "common_nb", 20240815L)

stopifnot(
  isTRUE(ci_data$success),
  isTRUE(mi_data$success),
  isTRUE(ordinary_data$success),
  identical(
    ci_data$dgp_schema,
    "continuous_nb_dispersion_estimated_theta"
  ),
  identical(ci_data$working_dispersion_mode, "estimated"),
  is.na(ci_data$working_nb_size),
  is.finite(ci_data$theta_hat),
  ci_data$theta_hat > 0,
  is.finite(mi_data$theta_hat),
  mi_data$theta_hat > 0,
  is.finite(ordinary_data$theta_hat),
  ordinary_data$theta_hat > 0,
  sd(ci_data$propensity) > 0,
  max(abs(mi_data$propensity - 0.14)) < 1e-12,
  all(is.finite(46 * exp(-2.1 * ci_data$Z))),
  all(46 * exp(-2.1 * ci_data$Z) > 0)
)

method_functions <- list(
  parametric_score = parametric_score_f,
  information_studentized_permutation =
    information_studentized_permutation_f,
  studentized_permutation = studentized_permutation_f,
  hgf_signflip_mc = hgf_signflip_mc_f,
  standardized_signflip = standardized_signflip_f
)
outputs <- lapply(method_functions, function(method) method(ci_data))

stopifnot(
  all(vapply(outputs, function(x) isTRUE(x$success), logical(1L))),
  all(vapply(
    outputs,
    function(x) identical(
      x$implementation_schema,
      "continuous_nb_dispersion_estimated_theta"
    ),
    logical(1L)
  )),
  all(vapply(
    outputs,
    function(x) is.finite(x$p.value) && x$p.value >= 0 && x$p.value <= 1,
    logical(1L)
  )),
  all(vapply(outputs, function(x) is.finite(x$statistic), logical(1L)))
)

# Recreating a dataset from the same seed must be bitwise stable.
ci_repeat <- make_dataset("CI", "severe_variance", 20240813L)
stopifnot(
  identical(ci_data$Y, ci_repeat$Y),
  identical(ci_data$X, ci_repeat$X),
  identical(ci_data$Z, ci_repeat$Z),
  identical(ci_data$perm_seed, ci_repeat$perm_seed),
  identical(ci_data$flip_seed, ci_repeat$flip_seed)
)

cat("Section 5.1 continuous-dispersion NB smoke test passed.\n")
