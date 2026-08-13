#!/usr/bin/env Rscript

# Lightweight contract test for the matched-calibration implementation used in
# Section 5.2.  This does not reproduce the simulation; it checks that one GCM
# and one PCM statistic each supply both their normal and permutation p-values,
# and that unrelated use of R's ambient random-number stream cannot change the
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

set.seed(7701L)
smoke_data <- generate_data(
  n = 40L,
  dgp = "setting1_xindep",
  gamma_signal = 0,
  B_perm = 3L
)

first_result <- matched_gcm_pcm_rf_f(smoke_data)
invisible(stats::runif(25L))
second_result <- matched_gcm_pcm_rf_f(smoke_data)

if (!identical(first_result, second_result)) {
  stop(
    "Matched calibration changed after unrelated ambient RNG use.",
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

cat("Section 5.2 matched-calibration smoke test passed.\n")
