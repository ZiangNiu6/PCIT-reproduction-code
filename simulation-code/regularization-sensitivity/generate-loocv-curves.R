#!/usr/bin/env Rscript

# Generate the LOOCV curves used in the locked-tuning simulation figure.
# The calibration rates themselves come from the independent R=1000,
# B=499 confirmation runs recorded in calibration-summary.csv. This script
# supplies a focused tuning diagnostic: mean exact-LOOCV loss for the observed
# X and the exact conditional mean of that loss over its raw-X permutation
# orbit in each of 1000 null datasets.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", script_argument))
script_directory <- dirname(script_file)
project_root <- normalizePath(file.path(script_directory, "..", "..", ".."))
results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
output_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)

n <- 300L
q_penalized <- as.integer(Sys.getenv(
  "PCIT_LOOCV_Q_PENALIZED",
  unset = "178"
))
unpenalized_rank <- 2L
signal_variance <- as.numeric(Sys.getenv(
  "PCIT_LOOCV_SIGNAL_VARIANCE",
  unset = "7.5"
))
signal_pattern <- Sys.getenv(
  "PCIT_LOOCV_SIGNAL_PATTERN",
  unset = "leading"
)
signal_rank <- as.integer(Sys.getenv(
  "PCIT_LOOCV_SIGNAL_RANK",
  unset = "16"
))
n_rep <- as.integer(Sys.getenv("PCIT_LOOCV_R", unset = "1000"))
curve_seed <- as.integer(Sys.getenv(
  "PCIT_LOOCV_SEED",
  unset = "20260804"
))
if (length(signal_variance) != 1L ||
    !is.finite(signal_variance) ||
    signal_variance < 0) {
  stop("PCIT_LOOCV_SIGNAL_VARIANCE must be one nonnegative number.")
}
if (length(q_penalized) != 1L ||
    is.na(q_penalized) ||
    q_penalized < 2L ||
    q_penalized >= n - unpenalized_rank ||
    q_penalized %% 2L != 0L) {
  stop("PCIT_LOOCV_Q_PENALIZED must be even and between 2 and n - 3.")
}
if (!signal_pattern %in% c("spread", "leading")) {
  stop("PCIT_LOOCV_SIGNAL_PATTERN must be 'spread' or 'leading'.")
}
if (length(signal_rank) != 1L ||
    is.na(signal_rank) ||
    signal_rank < 2L ||
    signal_rank > q_penalized ||
    signal_rank %% 2L != 0L) {
  stop("PCIT_LOOCV_SIGNAL_RANK must be even and between 2 and q_penalized.")
}
if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L) {
  stop("PCIT_LOOCV_R must be a positive integer.")
}
if (length(curve_seed) != 1L || is.na(curve_seed)) {
  stop("PCIT_LOOCV_SEED must be an integer.")
}

effective_dimension <- function(weights) {
  sum(weights^2)^2 / sum(weights^4)
}

find_penalty <- function(weight_function, target, label) {
  objective <- function(log_ratio) {
    effective_dimension(weight_function(exp(log_ratio))) - target
  }
  bracket <- c(-30, 30)
  endpoint_values <- vapply(bracket, objective, numeric(1))
  if (prod(sign(endpoint_values)) > 0) {
    stop(label, " cannot attain effective dimension ", target, ".")
  }
  exp(stats::uniroot(objective, bracket, tol = 1e-12)$root)
}

make_fourier_geometry <- function(n) {
  theta <- 2 * pi * (seq_len(n) - 1L) / n
  unpenalized <- cbind(
    intercept = rep(1 / sqrt(n), n),
    confounder = sqrt(2 / n) * cos(theta)
  )
  perpendicular <- list(sqrt(2 / n) * sin(theta))
  for (frequency in 2L:(n / 2L - 1L)) {
    perpendicular[[length(perpendicular) + 1L]] <-
      sqrt(2 / n) * cos(frequency * theta)
    perpendicular[[length(perpendicular) + 1L]] <-
      sqrt(2 / n) * sin(frequency * theta)
  }
  perpendicular[[length(perpendicular) + 1L]] <-
    cos((n / 2L) * theta) / sqrt(n)
  perpendicular <- do.call(cbind, perpendicular)
  complete_basis <- cbind(unpenalized, perpendicular)
  orthogonality_error <- max(abs(
    crossprod(complete_basis) - diag(n)
  ))
  if (orthogonality_error > 1e-10) {
    stop("The Fourier basis is not numerically orthogonal.")
  }
  list(
    unpenalized = unpenalized,
    perpendicular = perpendicular,
    confounder = sqrt(2) * cos(theta),
    orthogonality_error = orthogonality_error
  )
}

make_linear_geometry <- function(n, design_seed = 20260804L) {
  # Fixed Gaussian linear design.  The first orthogonal direction is an
  # unpenalized linear confounder; the next q_penalized directions are the
  # left singular vectors of the penalized design.  The singular values are
  # supplied by make_grid(), so the generated conditional mean is exactly
  # linear in the same covariates used by LRR.
  set.seed(design_seed)
  intercept <- rep(1 / sqrt(n), n)
  raw_design <- matrix(stats::rnorm(n * (n - 1L)), nrow = n)
  raw_design <- raw_design - intercept %*%
    crossprod(intercept, raw_design)
  orthogonal <- qr.Q(qr(raw_design, LAPACK = TRUE), complete = FALSE)
  unpenalized <- cbind(
    intercept = intercept,
    confounder = orthogonal[, 1L]
  )
  perpendicular <- orthogonal[, -1L, drop = FALSE]
  complete_basis <- cbind(unpenalized, perpendicular)
  orthogonality_error <- max(abs(
    crossprod(complete_basis) - diag(n)
  ))
  if (orthogonality_error > 1e-10) {
    stop("The fixed linear design basis is not numerically orthogonal.")
  }
  list(
    unpenalized = unpenalized,
    perpendicular = perpendicular,
    confounder = sqrt(n) * unpenalized[, "confounder"],
    orthogonality_error = orthogonality_error,
    design = "fixed Gaussian linear"
  )
}

make_grid <- function(geometry, learner) {
  residual_rank <- n - unpenalized_rank
  hard_residual_rank <- residual_rank - q_penalized
  eigenvalues <- if (learner == "LRR") {
    if (signal_pattern == "leading") {
      exp(seq(
        log(4 * n),
        log(0.25 * n),
        length.out = q_penalized
      ))
    } else {
      exp(seq(
        log(0.25 * n),
        log(4 * n),
        length.out = q_penalized
      ))
    }
  } else if (learner == "KRR") {
    exp(seq(log(1), log(0.01), length.out = q_penalized))
  } else {
    stop("Unknown learner: ", learner)
  }
  penalty_scale <- stats::median(eigenvalues)
  weight_function <- function(ratio) {
    lambda <- ratio * penalty_scale
    c(
      lambda / (eigenvalues + lambda),
      rep(1, hard_residual_rank)
    )
  }
  standard_dimensions <- c(seq(90, 290, by = 10), 295)
  interior_dimensions <- standard_dimensions[
    standard_dimensions > hard_residual_rank &
      standard_dimensions < residual_rank
  ]
  target_dimensions <- c(
    hard_residual_rank,
    interior_dimensions,
    residual_rank
  )
  interior_ratios <- vapply(interior_dimensions, function(target) {
    find_penalty(
      weight_function,
      target,
      paste0(learner, " tuning path")
    )
  }, numeric(1))
  penalty_ratios <- c(0, interior_ratios, Inf)
  weights <- do.call(cbind, lapply(penalty_ratios, function(ratio) {
    if (ratio == 0) {
      c(rep(0, q_penalized), rep(1, hard_residual_rank))
    } else if (is.infinite(ratio)) {
      rep(1, residual_rank)
    } else {
      weight_function(ratio)
    }
  }))
  attained_dimensions <- apply(weights, 2L, effective_dimension)
  if (max(abs(attained_dimensions - target_dimensions)) > 1e-6) {
    stop(learner, " effective-dimension path is incorrect.")
  }
  list(
    learner = learner,
    basis = geometry$perpendicular,
    active_basis = geometry$perpendicular[
      ,
      seq_len(q_penalized),
      drop = FALSE
    ],
    active_weights = weights[seq_len(q_penalized), , drop = FALSE],
    residual_diagonal = (geometry$perpendicular^2) %*% weights,
    residual_leverage = (geometry$perpendicular^2) %*% (weights^2),
    effective_dimension = attained_dimensions,
    penalty_ratio = penalty_ratios
  )
}

make_signal <- function(grid) {
  active_x <- if (signal_pattern == "leading") {
    seq.int(1L, signal_rank, by = 2L)
  } else {
    seq.int(1L, q_penalized, by = 2L)
  }
  coefficients <- rep(
    1 / sqrt(length(active_x)),
    length(active_x)
  )
  as.vector(
    sqrt(n * signal_variance) *
      grid$basis[, active_x, drop = FALSE] %*% coefficients
  )
}

loocv_score_matrix <- function(grid, geometry, values) {
  residual_infinite <- values - geometry$unpenalized %*%
    crossprod(geometry$unpenalized, values)
  active_coordinates <- crossprod(grid$active_basis, values)
  scores <- matrix(
    NA_real_,
    nrow = ncol(values),
    ncol = length(grid$effective_dimension)
  )
  standard_errors <- scores
  for (candidate in seq_along(grid$effective_dimension)) {
    adjustment <- grid$active_weights[, candidate] - 1
    residual <- residual_infinite + grid$active_basis %*% (
      adjustment * active_coordinates
    )
    loo_residual <- sweep(
      residual,
      1L,
      grid$residual_diagonal[, candidate],
      "/"
    )
    losses <- loo_residual^2
    scores[, candidate] <- colMeans(losses)
    standard_errors[, candidate] <- apply(
      losses,
      2L,
      stats::sd
    ) / sqrt(nrow(losses))
  }
  list(mean = scores, se = standard_errors)
}

permutation_orbit_score_matrix <- function(grid, values) {
  # For a fixed response vector x and a uniform permutation P, the covariance
  # of Px on the intercept-orthogonal subspace is var(x) I, where var uses the
  # n-1 denominator. Since every candidate annihilates the intercept, this is
  # the exact conditional permutation-orbit mean of squared LOOCV residuals.
  sample_variance <- apply(values, 2L, stats::var)
  orbit_multiplier <- colMeans(
    grid$residual_leverage / grid$residual_diagonal^2
  )
  outer(sample_variance, orbit_multiplier)
}

summarize_scores <- function(
    scores,
    grid,
    sample_type,
    within_sample_se = NULL) {
  selected_index <- max.col(-scores, ties.method = "last")
  one_se_index <- rep(NA_integer_, nrow(scores))
  if (!is.null(within_sample_se)) {
    row_index <- seq_len(nrow(scores))
    threshold <- scores[cbind(row_index, selected_index)] +
      within_sample_se[cbind(row_index, selected_index)]
    one_se_index <- max.col(
      scores <= threshold,
      ties.method = "last"
    )
  }
  score_mean <- colMeans(scores)
  score_se <- apply(scores, 2L, stats::sd) / sqrt(nrow(scores))
  relative_scores <- scores / apply(scores, 1L, min)
  relative_mean <- colMeans(relative_scores)
  relative_se <- apply(relative_scores, 2L, stats::sd) /
    sqrt(nrow(relative_scores))
  data.frame(
    learner = grid$learner,
    sample_type = sample_type,
    q_penalized = q_penalized,
    signal_variance = signal_variance,
    signal_pattern = signal_pattern,
    signal_rank = signal_rank,
    effective_dimension = grid$effective_dimension,
    penalty_ratio = grid$penalty_ratio,
    mean_loocv_mse = score_mean,
    se_loocv_mse = score_se,
    lower_95 = score_mean - 1.96 * score_se,
    upper_95 = score_mean + 1.96 * score_se,
    mean_within_sample_se = if (is.null(within_sample_se)) {
      NA_real_
    } else {
      colMeans(within_sample_se)
    },
    mean_relative_loocv_mse = relative_mean,
    se_relative_loocv_mse = relative_se,
    relative_lower_95 = relative_mean - 1.96 * relative_se,
    relative_upper_95 = relative_mean + 1.96 * relative_se,
    selection_rate = tabulate(
      selected_index,
      nbins = length(grid$effective_dimension)
    ) / length(selected_index),
    mean_selected_effective_dimension = mean(
      grid$effective_dimension[selected_index]
    ),
    one_se_selection_rate = if (all(is.na(one_se_index))) {
      NA_real_
    } else {
      tabulate(
        one_se_index,
        nbins = length(grid$effective_dimension)
      ) / length(one_se_index)
    },
    mean_one_se_effective_dimension = if (all(is.na(one_se_index))) {
      NA_real_
    } else {
      mean(grid$effective_dimension[one_se_index])
    },
    stringsAsFactors = FALSE
  )
}

geometries <- list(
  LRR = make_linear_geometry(n),
  KRR = make_fourier_geometry(n)
)
grids <- lapply(names(geometries), function(learner) {
  make_grid(geometries[[learner]], learner)
})
names(grids) <- names(geometries)

# The learner-specific conditional means are matched to their fitted design:
# linear in the fixed LRR covariates and smooth in the KRR eigenbasis. Common
# innovations reduce Monte Carlo noise without forcing the two DGPs to share a
# covariate geometry.
set.seed(curve_seed)
innovation_x <- matrix(rnorm(n * n_rep), nrow = n)

curve_rows <- lapply(names(grids), function(learner) {
  grid <- grids[[learner]]
  geometry <- geometries[[learner]]
  mean_x <- 0.5 + 2 * geometry$confounder + make_signal(grid)
  observed_x <- mean_x + innovation_x
  observed_loocv <- loocv_score_matrix(grid, geometry, observed_x)
  permutation_orbit_scores <- permutation_orbit_score_matrix(
    grid,
    observed_x
  )
  rbind(
    summarize_scores(
      observed_loocv$mean,
      grid,
      "Observed X",
      within_sample_se = observed_loocv$se
    ),
    summarize_scores(
      permutation_orbit_scores,
      grid,
      "Permutation orbit"
    )
  )
})
curves <- do.call(rbind, curve_rows)

output_file <- Sys.getenv(
  "PCIT_LOOCV_OUTPUT",
  unset = file.path(output_directory, "loocv-curves.csv")
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(curves, output_file, row.names = FALSE)

selection_summary <- unique(curves[, c(
  "learner",
  "sample_type",
  "mean_selected_effective_dimension",
  "mean_one_se_effective_dimension"
)])
print(selection_summary, row.names = FALSE, digits = 5)
message("Saved: ", normalizePath(output_file))
