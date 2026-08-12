#!/usr/bin/env Rscript

# Locked-only fixed-lambda KRR path confirmation.
#
# For each of 20 fixed residual-effective-dimension targets, this script uses
# the same fixed KRR residual operator for observed X, observed Y, and all
# raw-X permutations.  There is no cross-validation or tuning re-selection.

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
regularization_output_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)

parse_positive_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

n <- 300L
unpenalized_rank <- 2L
residual_rank <- n - unpenalized_rank
q_penalized <- 178L
q_total <- q_penalized + unpenalized_rank
hard_residual_rank <- residual_rank - q_penalized
signal_rank <- 16L
signal_variance <- 7.5
confounding_strength <- 2
nominal_level <- 0.05
critical_value <- stats::qnorm(1 - nominal_level / 2)

two_sided_normal_probability <- function(location, scale, cutoff) {
  if (any(!is.finite(location)) ||
      any(!is.finite(scale)) ||
      any(scale <= 0) ||
      any(!is.finite(cutoff)) ||
      any(cutoff <= 0)) {
    stop("Invalid Gaussian approximation parameter.")
  }
  stats::pnorm((-cutoff - location) / scale) +
    stats::pnorm((location - cutoff) / scale)
}

R <- parse_positive_integer("PCIT_FIXED_PATH_R", "1000")
B <- parse_positive_integer("PCIT_FIXED_PATH_B", "499")
if (B < 19L) {
  stop("PCIT_FIXED_PATH_B must be at least 19.")
}
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(detected_cores) != 1L ||
    !is.finite(detected_cores) ||
    detected_cores < 1L) {
  detected_cores <- 1L
}
cores <- parse_positive_integer(
  "PCIT_FIXED_PATH_CORES",
  as.character(min(4L, detected_cores))
)
output_file <- Sys.getenv(
  "PCIT_FIXED_PATH_OUTPUT",
  unset = file.path(regularization_output_directory, "krr-fixed-path-summary.csv")
)
raw_file <- sub("\\.csv$", "-raw.rds", output_file)
if (identical(raw_file, output_file)) {
  raw_file <- paste0(output_file, "-raw.rds")
}
simulation_seed <- 20260803L + 2400000L
cell_seed <- simulation_seed + q_penalized * 1009L

effective_dimension <- function(weights) {
  sum(weights^2)^2 / sum(weights^4)
}

find_penalty <- function(weight_function, target) {
  objective <- function(log_ratio) {
    effective_dimension(weight_function(exp(log_ratio))) - target
  }
  endpoints <- c(-30, 30)
  endpoint_values <- vapply(endpoints, objective, numeric(1))
  if (prod(sign(endpoint_values)) > 0) {
    stop("Cannot attain residual effective dimension ", target, ".")
  }
  exp(stats::uniroot(objective, endpoints, tol = 1e-12)$root)
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
  orthogonality_error <- max(abs(
    crossprod(cbind(unpenalized, perpendicular)) - diag(n)
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

make_krr_grid <- function(geometry) {
  eigenvalues <- exp(seq(
    log(1),
    log(0.01),
    length.out = q_penalized
  ))
  penalty_scale <- stats::median(eigenvalues)
  weight_function <- function(ratio) {
    lambda <- ratio * penalty_scale
    c(
      lambda / (eigenvalues + lambda),
      rep(1, hard_residual_rank)
    )
  }
  target_dimensions <- c(
    seq.int(hard_residual_rank, 290L, by = 10L),
    295,
    residual_rank
  )
  interior_targets <- target_dimensions[-c(1L, length(target_dimensions))]
  interior_ratios <- vapply(interior_targets, function(target) {
    find_penalty(weight_function, target)
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
  dimension_error <- max(abs(attained_dimensions - target_dimensions))
  if (dimension_error > 1e-6) {
    stop("The KRR path failed its effective-dimension check.")
  }
  residual_diagonal <- (geometry$perpendicular^2) %*% weights
  residual_leverage <- (geometry$perpendicular^2) %*% (weights^2)
  if (any(!is.finite(residual_diagonal)) ||
      any(residual_diagonal <= sqrt(.Machine$double.eps))) {
    stop("The KRR path has an unstable residual diagonal.")
  }
  list(
    basis = geometry$perpendicular,
    active_basis = geometry$perpendicular[
      , seq_len(q_penalized), drop = FALSE
    ],
    hard_basis = geometry$perpendicular[
      , q_penalized + seq_len(hard_residual_rank), drop = FALSE
    ],
    weights = weights,
    active_weights = weights[seq_len(q_penalized), , drop = FALSE],
    residual_diagonal = residual_diagonal,
    residual_leverage = residual_leverage,
    target_dimensions = target_dimensions,
    effective_dimensions = attained_dimensions,
    penalty_ratios = penalty_ratios,
    penalty_absolute = penalty_ratios * penalty_scale,
    penalty_scale = penalty_scale,
    eigenvalues = eigenvalues,
    dimension_error = dimension_error
  )
}

make_orthogonal_signals <- function(grid) {
  active_x <- seq.int(1L, signal_rank, by = 2L)
  active_y <- seq.int(2L, signal_rank, by = 2L)
  scale <- sqrt(n * signal_variance)
  list(
    x = as.vector(
      scale * grid$basis[, active_x, drop = FALSE] %*%
        rep(1 / sqrt(length(active_x)), length(active_x))
    ),
    y = as.vector(
      scale * grid$basis[, active_y, drop = FALSE] %*%
        rep(1 / sqrt(length(active_y)), length(active_y))
    )
  )
}

make_permutations <- function(values, permutation_count) {
  draws <- matrix(NA_real_, nrow = length(values), ncol = permutation_count)
  for (draw_index in seq_len(permutation_count)) {
    draws[, draw_index] <- values[sample.int(length(values))]
  }
  draws
}

studentized_statistic <- function(residual_x, residual_y) {
  scores <- residual_x * residual_y
  numerator <- sum(scores)
  numerator / sqrt(max(
    sum(scores^2) - numerator^2 / length(scores),
    .Machine$double.xmin
  ))
}

exact_interval <- function(rejection_count, trials, confidence = 0.95) {
  unname(stats::binom.test(
    rejection_count,
    trials,
    conf.level = confidence
  )$conf.int)
}

geometry <- make_fourier_geometry(n)
grid <- make_krr_grid(geometry)
signals <- make_orthogonal_signals(grid)
mean_x <- 0.5 + confounding_strength * geometry$confounder + signals$x
mean_y <- -0.5 + confounding_strength * geometry$confounder + signals$y
candidate_count <- length(grid$effective_dimensions)

# Population annihilation/reference benchmark, integrating over the Y noise.
mean_y_coordinates <- as.vector(crossprod(grid$basis, mean_y))
mean_x_coordinates <- as.vector(crossprod(grid$basis, mean_x))
population_residual_x <- grid$basis %*% (
  grid$weights * mean_x_coordinates
)
population_residual_y <- grid$basis %*% (
  grid$weights * mean_y_coordinates
)
population_numerator_variance <- colSums(
  (grid$weights^2 * mean_y_coordinates)^2
) + colSums(grid$weights^4)
population_studentizer <- colSums(
  grid$residual_leverage *
    (population_residual_y^2 + grid$residual_leverage)
) - population_numerator_variance / n
population_gamma <- population_numerator_variance / population_studentizer
population_alpha <- 2 * stats::pnorm(
  -critical_value / sqrt(population_gamma)
)

one_replicate <- function(replicate_index) {
  set.seed(cell_seed + replicate_index * 7919L)
  X <- mean_x + stats::rnorm(n)
  Y <- mean_y + stats::rnorm(n)
  permuted_x <- make_permutations(X, B)

  coordinates_x <- as.vector(crossprod(grid$basis, X))
  coordinates_y <- as.vector(crossprod(grid$basis, Y))
  residual_x <- grid$basis %*% (grid$weights * coordinates_x)
  residual_y <- grid$basis %*% (grid$weights * coordinates_y)

  observed_statistics <- vapply(seq_len(candidate_count), function(candidate) {
    studentized_statistic(
      residual_x[, candidate],
      residual_y[, candidate]
    )
  }, numeric(1))

  # Transform every raw-X permutation once.  The active contributions for all
  # 20 fixed filters are then evaluated in one large matrix multiplication.
  permutation_coordinates <- crossprod(grid$basis, permuted_x)
  active_coordinates <- permutation_coordinates[
    seq_len(q_penalized), , drop = FALSE
  ]
  hard_coordinates <- permutation_coordinates[
    q_penalized + seq_len(hard_residual_rank), , drop = FALSE
  ]
  hard_residuals <- grid$hard_basis %*% hard_coordinates
  weighted_active_coordinates <- do.call(cbind, lapply(
    seq_len(candidate_count),
    function(candidate) {
      grid$active_weights[, candidate] * active_coordinates
    }
  ))
  active_residuals <- grid$active_basis %*% weighted_active_coordinates

  locked_pvalues <- numeric(candidate_count)
  locked_q95 <- numeric(candidate_count)
  max_direct_error <- 0
  for (candidate in seq_len(candidate_count)) {
    columns <- (candidate - 1L) * B + seq_len(B)
    permutation_residuals <- hard_residuals +
      active_residuals[, columns, drop = FALSE]
    scores <- permutation_residuals * residual_y[, candidate]
    numerators <- colSums(scores)
    reference_statistics <- numerators / sqrt(pmax(
      colSums(scores^2) - numerators^2 / n,
      .Machine$double.xmin
    ))
    locked_pvalues[candidate] <- (
      1 + sum(abs(reference_statistics) >= abs(observed_statistics[candidate]))
    ) / (B + 1)
    locked_q95[candidate] <- unname(stats::quantile(
      abs(reference_statistics),
      0.95,
      type = 8
    ))

    if (replicate_index == 1L && candidate %in% c(1L, 10L, candidate_count)) {
      direct_residual <- grid$basis %*% (
        grid$weights[, candidate] *
          crossprod(grid$basis, permuted_x[, 1L])
      )
      max_direct_error <- max(
        max_direct_error,
        max(abs(direct_residual - permutation_residuals[, 1L]))
      )
    }
  }

  projected_y_coordinates <- grid$weights^2 * coordinates_y
  numerator_variance <- colSums(projected_y_coordinates^2)
  expected_studentizer <- colSums(
    grid$residual_leverage * residual_y^2
  ) - numerator_variance / n
  gamma_plugin <- numerator_variance / expected_studentizer
  alpha_plugin <- 2 * stats::pnorm(
    -critical_value / sqrt(gamma_plugin)
  )
  conditional_mean_scores <- population_residual_x * residual_y
  conditional_numerator_mean <- colSums(conditional_mean_scores)
  conditional_mean_studentizer <- colSums(
    conditional_mean_scores^2
  ) - conditional_numerator_mean^2 / n
  observed_expected_studentizer <-
    expected_studentizer + conditional_mean_studentizer
  observed_location <- conditional_numerator_mean /
    sqrt(observed_expected_studentizer)
  observed_scale <- sqrt(
    numerator_variance / observed_expected_studentizer
  )
  normal_theory_dgp <- two_sided_normal_probability(
    observed_location,
    observed_scale,
    critical_value
  )
  permutation_theory_dgp <- two_sided_normal_probability(
    observed_location,
    observed_scale,
    critical_value * sqrt(gamma_plugin)
  )

  list(
    normal_reject = abs(observed_statistics) > critical_value,
    locked_reject = locked_pvalues <= nominal_level,
    observed_statistic = observed_statistics,
    locked_pvalue = locked_pvalues,
    locked_q95 = locked_q95,
    gamma_plugin = gamma_plugin,
    alpha_plugin = alpha_plugin,
    normal_theory_dgp = normal_theory_dgp,
    permutation_theory_dgp = permutation_theory_dgp,
    max_direct_error = max_direct_error
  )
}

message(
  "Running fixed KRR path: R=", R,
  ", B=", B,
  ", candidates=", candidate_count,
  ", q_penalized=", q_penalized,
  ", rank=", signal_rank,
  ", V=", signal_variance,
  ", cores=", cores, "."
)

start_time <- proc.time()[["elapsed"]]
replicate_indices <- seq_len(R)
batch_size <- max(cores * 5L, min(50L, R))
batches <- split(
  replicate_indices,
  ceiling(replicate_indices / batch_size)
)
replicate_results <- vector("list", R)
for (batch_index in seq_along(batches)) {
  batch <- batches[[batch_index]]
  batch_results <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(
      batch,
      one_replicate,
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(batch, one_replicate)
  }
  replicate_results[batch] <- batch_results
  message(
    "Completed ", max(batch), "/", R,
    " replicates (batch ", batch_index, "/", length(batches), ")."
  )
}
simulation_elapsed_seconds <- proc.time()[["elapsed"]] - start_time

bind_result <- function(name) {
  do.call(rbind, lapply(replicate_results, `[[`, name))
}
normal_reject <- bind_result("normal_reject")
locked_reject <- bind_result("locked_reject")
observed_statistic <- bind_result("observed_statistic")
locked_pvalue <- bind_result("locked_pvalue")
locked_q95 <- bind_result("locked_q95")
gamma_plugin <- bind_result("gamma_plugin")
alpha_plugin <- bind_result("alpha_plugin")
normal_theory_dgp <- bind_result("normal_theory_dgp")
permutation_theory_dgp <- bind_result("permutation_theory_dgp")
max_direct_error <- max(vapply(
  replicate_results,
  `[[`,
  numeric(1),
  "max_direct_error"
))

summary_rows <- lapply(seq_len(candidate_count), function(candidate) {
  normal_count <- sum(normal_reject[, candidate])
  locked_count <- sum(locked_reject[, candidate])
  normal_interval <- exact_interval(normal_count, R)
  locked_interval <- exact_interval(locked_count, R)
  data.frame(
    learner = "KRR",
    tuning_regime = "fixed lambda",
    candidate = candidate,
    n = n,
    q_penalized = q_penalized,
    q_total = q_total,
    hard_residual_rank = hard_residual_rank,
    signal_pattern = "leading",
    signal_rank = signal_rank,
    signal_variance = signal_variance,
    target_effective_dimension = grid$target_dimensions[candidate],
    effective_dimension = grid$effective_dimensions[candidate],
    penalty_ratio = grid$penalty_ratios[candidate],
    penalty_absolute = grid$penalty_absolute[candidate],
    R = R,
    B = B,
    normal_rejections = normal_count,
    normal_empirical = normal_count / R,
    normal_lower_95 = normal_interval[1L],
    normal_upper_95 = normal_interval[2L],
    normal_annihilation_plugin = mean(alpha_plugin[, candidate]),
    normal_theory_dgp = mean(normal_theory_dgp[, candidate]),
    normal_annihilation_plugin_mc_se = stats::sd(
      alpha_plugin[, candidate]
    ) / sqrt(R),
    normal_annihilation_population = population_alpha[candidate],
    mean_gamma_reference_plugin = mean(gamma_plugin[, candidate]),
    median_gamma_reference_plugin = stats::median(
      gamma_plugin[, candidate]
    ),
    locked_permutation_rejections = locked_count,
    locked_permutation_empirical = locked_count / R,
    locked_permutation_lower_95 = locked_interval[1L],
    locked_permutation_upper_95 = locked_interval[2L],
    nominal_level = nominal_level,
    permutation_theory_dgp = mean(
      permutation_theory_dgp[, candidate]
    ),
    observed_statistic_sd = stats::sd(observed_statistic[, candidate]),
    mean_locked_pvalue = mean(locked_pvalue[, candidate]),
    mean_locked_abs_q95 = mean(locked_q95[, candidate]),
    simulation_elapsed_seconds = simulation_elapsed_seconds,
    orthogonality_error = geometry$orthogonality_error,
    effective_dimension_max_error = grid$dimension_error,
    batched_residual_max_error = max_direct_error,
    seed = simulation_seed,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)

validation <- list(
  row_count = nrow(summary),
  candidate_count = candidate_count,
  effective_dimension_strictly_increasing = all(
    diff(summary$effective_dimension) > 0
  ),
  effective_dimension_max_error = grid$dimension_error,
  orthogonality_error = geometry$orthogonality_error,
  batched_residual_max_error = max_direct_error,
  finite_statistics = all(is.finite(observed_statistic)),
  finite_plugin = all(is.finite(alpha_plugin)),
  finite_dgp_theory = all(is.finite(normal_theory_dgp)) &&
    all(is.finite(permutation_theory_dgp)),
  pvalues_on_mc_grid = all(abs(
    locked_pvalue * (B + 1) - round(locked_pvalue * (B + 1))
  ) < 1e-10),
  nominal_level_is_constant = all(summary$nominal_level == 0.05)
)
if (!identical(validation$row_count, 20L) ||
    !validation$effective_dimension_strictly_increasing ||
    validation$effective_dimension_max_error > 1e-6 ||
    validation$orthogonality_error > 1e-10 ||
    validation$batched_residual_max_error > 1e-10 ||
    !validation$finite_statistics ||
    !validation$finite_plugin ||
    !validation$finite_dgp_theory ||
    !validation$pvalues_on_mc_grid ||
    !validation$nominal_level_is_constant) {
  stop("One or more fixed-path validation checks failed.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, output_file, row.names = FALSE)
saveRDS(
  list(
    summary = summary,
    validation = validation,
    raw = list(
      normal_reject = normal_reject,
      locked_reject = locked_reject,
      observed_statistic = observed_statistic,
      locked_pvalue = locked_pvalue,
      locked_q95 = locked_q95,
      gamma_plugin = gamma_plugin,
      alpha_plugin = alpha_plugin,
      normal_theory_dgp = normal_theory_dgp,
      permutation_theory_dgp = permutation_theory_dgp
    ),
    grid = list(
      target_effective_dimension = grid$target_dimensions,
      effective_dimension = grid$effective_dimensions,
      penalty_ratio = grid$penalty_ratios,
      penalty_absolute = grid$penalty_absolute
    ),
    design = list(
      n = n,
      q_penalized = q_penalized,
      q_total = q_total,
      hard_residual_rank = hard_residual_rank,
      signal_pattern = "leading",
      signal_rank = signal_rank,
      signal_variance = signal_variance,
      R = R,
      B = B,
      seed = simulation_seed
    ),
    simulation_elapsed_seconds = simulation_elapsed_seconds
  ),
  raw_file
)

print(summary[, c(
  "candidate",
  "effective_dimension",
  "penalty_ratio",
  "normal_empirical",
  "normal_theory_dgp",
  "locked_permutation_empirical",
  "permutation_theory_dgp"
)], row.names = FALSE, digits = 5)
print(validation)
message("Simulation elapsed seconds: ", round(simulation_elapsed_seconds, 2))
message("Saved: ", output_file)
message("Saved: ", raw_file)
