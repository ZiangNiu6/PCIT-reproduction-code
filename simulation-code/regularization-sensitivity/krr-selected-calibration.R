#!/usr/bin/env Rscript

# Paired confirmation of minimum-MSE and +1SE exact-LOOCV rules for the
# rank-q KRR simulation. For each rule, this script compares normal
# calibration with locked-tuning raw-X permutation: select the X- and Y-side
# candidates once from the observed data, then hold the observed X-side
# candidate fixed for every raw-X permutation.
#
# Default experiment:
#   n = 300, q_penalized = 178, q_total = 180, rank = 16, V = 7.5,
#   R = 1000 Monte Carlo datasets, B = 499 raw-X permutations.
#
# The +1SE rule selects the largest-penalty (equivalently, largest residual
# effective-dimension) candidate whose exact-LOOCV mean squared error is no
# greater than the minimum-MSE value plus the standard error of the 300
# per-observation squared LOOCV residuals at the minimizing candidate.
#
# Optional environment variables:
#   PCIT_KRR_CONFIRM_Q_PENALIZED=178
#   PCIT_KRR_CONFIRM_SIGNAL_VARIANCE=7.5
#   PCIT_KRR_CONFIRM_SIGNAL_PATTERN=leading
#   PCIT_KRR_CONFIRM_SIGNAL_RANK=16
#   PCIT_KRR_CONFIRM_R=1000
#   PCIT_KRR_CONFIRM_B=499
#   PCIT_KRR_CONFIRM_CORES=4
#   PCIT_KRR_CONFIRM_OUTPUT_PREFIX=/path/to/krr-min-1se

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
  if (length(value) != 1L || !is.finite(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

n <- 300L
q_penalized <- parse_positive_integer(
  "PCIT_KRR_CONFIRM_Q_PENALIZED",
  "178"
)
unpenalized_rank <- 2L
q_total <- q_penalized + unpenalized_rank
signal_variance <- as.numeric(Sys.getenv(
  "PCIT_KRR_CONFIRM_SIGNAL_VARIANCE",
  unset = "7.5"
))
signal_pattern <- Sys.getenv(
  "PCIT_KRR_CONFIRM_SIGNAL_PATTERN",
  unset = "leading"
)
signal_rank <- as.integer(Sys.getenv(
  "PCIT_KRR_CONFIRM_SIGNAL_RANK",
  unset = "16"
))
confounding_strength <- 2
R <- parse_positive_integer("PCIT_KRR_CONFIRM_R", "1000")
B <- parse_positive_integer("PCIT_KRR_CONFIRM_B", "499")
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(detected_cores) != 1L ||
    !is.finite(detected_cores) ||
    detected_cores < 1L) {
  detected_cores <- 1L
}
cores <- parse_positive_integer(
  "PCIT_KRR_CONFIRM_CORES",
  as.character(min(4L, detected_cores))
)
output_prefix <- Sys.getenv(
  "PCIT_KRR_CONFIRM_OUTPUT_PREFIX",
  unset = file.path(regularization_output_directory, "krr-selected")
)
nominal_level <- 0.05
critical_value <- stats::qnorm(1 - nominal_level / 2)
if (length(signal_variance) != 1L ||
    !is.finite(signal_variance) ||
    signal_variance < 0) {
  stop("PCIT_KRR_CONFIRM_SIGNAL_VARIANCE must be one nonnegative number.")
}
if (q_penalized >= n - unpenalized_rank || q_penalized %% 2L != 0L) {
  stop(
    "PCIT_KRR_CONFIRM_Q_PENALIZED must be even and smaller than n - ",
    "unpenalized_rank."
  )
}
if (!signal_pattern %in% c("spread", "leading")) {
  stop("PCIT_KRR_CONFIRM_SIGNAL_PATTERN must be 'spread' or 'leading'.")
}
if (length(signal_rank) != 1L ||
    is.na(signal_rank) ||
    signal_rank < 2L ||
    signal_rank > q_penalized ||
    signal_rank %% 2L != 0L) {
  stop(
    "PCIT_KRR_CONFIRM_SIGNAL_RANK must be even and between 2 and ",
    "q_penalized."
  )
}

# Pair the minimum-CV and +1SE rules on identical datasets and permutations.
simulation_seed <- 20260803L + 2400000L
cell_seed <- simulation_seed + q_penalized * 1009L

effective_dimension <- function(weights) {
  sum(weights^2)^2 / sum(weights^4)
}

find_penalty <- function(weight_function, target, label) {
  objective <- function(log_ratio) {
    effective_dimension(weight_function(exp(log_ratio))) - target
  }
  endpoints <- c(-30, 30)
  endpoint_values <- vapply(endpoints, objective, numeric(1))
  if (prod(sign(endpoint_values)) > 0) {
    stop(label, " cannot attain effective dimension ", target, ".")
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
  complete_basis <- cbind(unpenalized, perpendicular)
  orthogonality_error <- max(abs(
    crossprod(complete_basis) - diag(n)
  ))
  if (orthogonality_error > 1e-10) {
    stop("The deterministic Fourier basis is not numerically orthogonal.")
  }
  list(
    unpenalized = unpenalized,
    perpendicular = perpendicular,
    confounder = sqrt(2) * cos(theta),
    orthogonality_error = orthogonality_error
  )
}

make_krr_grid <- function(geometry) {
  residual_rank <- n - unpenalized_rank
  hard_residual_rank <- residual_rank - q_penalized
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
    find_penalty(weight_function, target, "KRR exact-LOOCV path")
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
    stop("The KRR candidate path failed its effective-dimension check.")
  }
  residual_diagonal <- (geometry$perpendicular^2) %*% weights
  if (any(!is.finite(residual_diagonal)) ||
      any(residual_diagonal <= sqrt(.Machine$double.eps))) {
    stop("The KRR path has an unstable exact-LOOCV denominator.")
  }
  list(
    basis = geometry$perpendicular,
    weights = weights,
    residual_diagonal = residual_diagonal,
    residual_leverage = (geometry$perpendicular^2) %*% (weights^2),
    candidate = seq_along(penalty_ratios),
    effective_dimension = attained_dimensions,
    penalty_ratio = penalty_ratios,
    eigenvalue = eigenvalues,
    hard_residual_rank = hard_residual_rank
  )
}

make_orthogonal_signals <- function(grid) {
  active_limit <- if (signal_pattern == "leading") {
    signal_rank
  } else {
    q_penalized
  }
  active_x <- seq.int(1L, active_limit, by = 2L)
  active_y <- seq.int(2L, active_limit, by = 2L)
  coefficient_x <- rep(1 / sqrt(length(active_x)), length(active_x))
  coefficient_y <- rep(1 / sqrt(length(active_y)), length(active_y))
  scale <- sqrt(n * signal_variance)
  list(
    x = as.vector(
      scale * grid$basis[, active_x, drop = FALSE] %*% coefficient_x
    ),
    y = as.vector(
      scale * grid$basis[, active_y, drop = FALSE] %*% coefficient_y
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

apply_candidate <- function(grid, candidate, values) {
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  coordinates <- crossprod(grid$basis, values)
  grid$basis %*% (grid$weights[, candidate] * coordinates)
}

column_statistics <- function(residual_x, residual_y) {
  scores <- residual_x * residual_y
  numerator <- colSums(scores)
  denominator_squared <- colSums(scores^2) -
    numerator^2 / nrow(scores)
  numerator / sqrt(pmax(denominator_squared, .Machine$double.xmin))
}

mc_pvalue <- function(reference, observed) {
  (1 + sum(abs(reference) >= abs(observed))) /
    (length(reference) + 1)
}

selected_residuals <- function(candidate_residuals, selected_indices) {
  sample_count <- dim(candidate_residuals)[2L]
  output <- matrix(NA_real_, nrow = n, ncol = sample_count)
  for (candidate in seq_len(dim(candidate_residuals)[3L])) {
    columns <- which(selected_indices == candidate)
    if (length(columns) > 0L) {
      output[, columns] <- candidate_residuals[
        , columns, candidate, drop = FALSE
      ][, , 1L]
    }
  }
  output
}

select_loocv_rules <- function(grid, values) {
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  sample_count <- ncol(values)
  candidate_count <- length(grid$candidate)
  coordinates <- crossprod(grid$basis, values)
  scores <- matrix(NA_real_, nrow = sample_count, ncol = candidate_count)
  score_se <- scores
  candidate_residuals <- array(
    NA_real_,
    dim = c(n, sample_count, candidate_count)
  )

  for (candidate in grid$candidate) {
    residual <- grid$basis %*% (
      grid$weights[, candidate] * coordinates
    )
    candidate_residuals[, , candidate] <- residual
    leave_one_out_residual <- sweep(
      residual,
      1L,
      grid$residual_diagonal[, candidate],
      "/"
    )
    losses <- leave_one_out_residual^2
    score <- colMeans(losses)
    loss_sum_of_squares <- colSums(losses^2)
    loss_variance <- pmax(
      (loss_sum_of_squares - n * score^2) / (n - 1),
      0
    )
    scores[, candidate] <- score
    score_se[, candidate] <- sqrt(loss_variance / n)
  }

  row_index <- seq_len(sample_count)
  minimum_score <- apply(scores, 1L, min)
  min_tolerance <- 1e-12 * pmax(1, abs(minimum_score))
  minimum_eligible <- sweep(
    scores,
    1L,
    minimum_score + min_tolerance,
    "<="
  )
  min_index <- max.col(minimum_eligible, ties.method = "last")

  one_se_threshold <- scores[cbind(row_index, min_index)] +
    score_se[cbind(row_index, min_index)]
  one_se_tolerance <- 1e-12 * pmax(1, abs(one_se_threshold))
  one_se_eligible <- sweep(
    scores,
    1L,
    one_se_threshold + one_se_tolerance,
    "<="
  )
  one_se_index <- max.col(one_se_eligible, ties.method = "last")

  list(
    score = scores,
    score_se = score_se,
    candidate_residuals = candidate_residuals,
    rules = list(
      `Minimum CV` = list(
        index = min_index,
        residual = selected_residuals(candidate_residuals, min_index),
        threshold = minimum_score
      ),
      `+1SE` = list(
        index = one_se_index,
        residual = selected_residuals(candidate_residuals, one_se_index),
        threshold = one_se_threshold
      )
    )
  )
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

one_replicate <- function(replicate_index) {
  set.seed(cell_seed + replicate_index * 7919L)
  X <- mean_x + stats::rnorm(n)
  Y <- mean_y + stats::rnorm(n)
  selection <- select_loocv_rules(grid, cbind(Y, X))
  permuted_x <- make_permutations(X, B)

  rule_rows <- lapply(names(selection$rules), function(rule_name) {
    rule <- selection$rules[[rule_name]]
    y_index <- rule$index[1L]
    x_index <- rule$index[2L]
    residual_y <- rule$residual[, 1L]
    residual_x <- rule$residual[, 2L]
    locked_permutation_residuals <- apply_candidate(
      grid,
      x_index,
      permuted_x
    )

    observed_statistic <- column_statistics(
      matrix(residual_x, ncol = 1L),
      matrix(residual_y, ncol = 1L)
    )[[1L]]
    locked_reference <- column_statistics(
      locked_permutation_residuals,
      residual_y
    )
    projected_y <- as.vector(grid$basis %*% (
      grid$weights[, x_index] * crossprod(grid$basis, residual_y)
    ))
    numerator_variance <- sum(projected_y^2)
    expected_studentizer <- sum(
      grid$residual_leverage[, x_index] * residual_y^2
    ) - numerator_variance / n
    gamma_plugin <- numerator_variance / expected_studentizer
    observed_dimension <- grid$effective_dimension[x_index]

    data.frame(
      replicate = replicate_index,
      cv_rule = rule_name,
      normal_reject = abs(observed_statistic) > critical_value,
      locked_reject = mc_pvalue(
        locked_reference,
        observed_statistic
      ) <= nominal_level,
      observed_statistic = observed_statistic,
      gamma_plugin = gamma_plugin,
      alpha_normal_plugin = 2 * stats::pnorm(
        -critical_value / sqrt(gamma_plugin)
      ),
      observed_x_candidate = x_index,
      observed_y_candidate = y_index,
      observed_x_dimension = observed_dimension,
      observed_y_dimension = grid$effective_dimension[y_index],
      locked_abs_q95 = unname(stats::quantile(
        abs(locked_reference),
        0.95,
        type = 8
      )),
      stringsAsFactors = FALSE
    )
  })

  selection_rows <- lapply(names(selection$rules), function(rule_name) {
    indices <- selection$rules[[rule_name]]$index
    roles <- c("Observed Y", "Observed X")
    counts <- do.call(rbind, lapply(unique(roles), function(role) {
      role_indices <- indices[roles == role]
      data.frame(
        replicate = replicate_index,
        cv_rule = rule_name,
        sample_role = role,
        candidate = grid$candidate,
        effective_dimension = grid$effective_dimension,
        count = tabulate(
          role_indices,
          nbins = length(grid$candidate)
        ),
        denominator = length(role_indices),
        stringsAsFactors = FALSE
      )
    }))
    counts
  })

  list(
    diagnostics = do.call(rbind, rule_rows),
    selection = do.call(rbind, selection_rows)
  )
}

message(
  "Running paired KRR minimum-CV/+1SE locked-tuning confirmation: R=", R,
  ", B=", B, ", n=", n, ", q_penalized=", q_penalized,
  ", V=", signal_variance, ", cores=", cores, "."
)

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

raw_diagnostics <- do.call(rbind, lapply(
  replicate_results,
  `[[`,
  "diagnostics"
))
raw_selection <- do.call(rbind, lapply(
  replicate_results,
  `[[`,
  "selection"
))

calibration_map <- list(
  `Normal calibration` = "normal_reject",
  `Locked-tuning refit` = "locked_reject"
)
calibration_summary <- do.call(rbind, lapply(
  unique(raw_diagnostics$cv_rule),
  function(rule_name) {
    rule_data <- raw_diagnostics[
      raw_diagnostics$cv_rule == rule_name,
      ,
      drop = FALSE
    ]
    do.call(rbind, lapply(names(calibration_map), function(method_name) {
      rejected <- rule_data[[calibration_map[[method_name]]]]
      rejection_count <- sum(rejected)
      interval <- exact_interval(rejection_count, length(rejected))
      data.frame(
        learner = "KRR",
        cv_rule = rule_name,
        method = method_name,
        rejections = rejection_count,
        R = length(rejected),
        B = if (method_name == "Normal calibration") NA_integer_ else B,
        estimate = mean(rejected),
        lower_95 = interval[[1L]],
        upper_95 = interval[[2L]],
        stringsAsFactors = FALSE
      )
    }))
  }
))

tuning_diagnostics <- do.call(rbind, lapply(
  unique(raw_diagnostics$cv_rule),
  function(rule_name) {
    rule_data <- raw_diagnostics[
      raw_diagnostics$cv_rule == rule_name,
      ,
      drop = FALSE
    ]
    data.frame(
      learner = "KRR",
      cv_rule = rule_name,
      mean_observed_x_effective_dimension = mean(
        rule_data$observed_x_dimension
      ),
      median_observed_x_effective_dimension = stats::median(
        rule_data$observed_x_dimension
      ),
      mean_observed_y_effective_dimension = mean(
        rule_data$observed_y_dimension
      ),
      median_observed_y_effective_dimension = stats::median(
        rule_data$observed_y_dimension
      ),
      observed_x_hard_fit_rate = mean(
        rule_data$observed_x_dimension == grid$hard_residual_rank
      ),
      observed_y_hard_fit_rate = mean(
        rule_data$observed_y_dimension == grid$hard_residual_rank
      ),
      mean_gamma_plugin = mean(rule_data$gamma_plugin),
      median_gamma_plugin = stats::median(rule_data$gamma_plugin),
      mean_alpha_normal_plugin = mean(rule_data$alpha_normal_plugin),
      observed_abs_q95 = unname(stats::quantile(
        abs(rule_data$observed_statistic),
        0.95,
        type = 8
      )),
      mean_locked_abs_q95 = mean(rule_data$locked_abs_q95),
      stringsAsFactors = FALSE
    )
  }
))

selection_summary <- aggregate(
  cbind(count, denominator) ~
    cv_rule + sample_role + candidate + effective_dimension,
  data = raw_selection,
  FUN = sum
)
selection_summary$selection_rate <-
  selection_summary$count / selection_summary$denominator
selection_summary$learner <- "KRR"
selection_summary <- selection_summary[, c(
  "learner",
  "cv_rule",
  "sample_role",
  "candidate",
  "effective_dimension",
  "count",
  "denominator",
  "selection_rate"
)]

candidate_grid <- data.frame(
  learner = "KRR",
  candidate = grid$candidate,
  effective_dimension = grid$effective_dimension,
  penalty_ratio = grid$penalty_ratio,
  stringsAsFactors = FALSE
)

output_directory <- dirname(output_prefix)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  calibration_summary,
  paste0(output_prefix, "-calibration-summary.csv"),
  row.names = FALSE
)
write.csv(
  tuning_diagnostics,
  paste0(output_prefix, "-tuning-diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  selection_summary,
  paste0(output_prefix, "-selection.csv"),
  row.names = FALSE
)
write.csv(
  candidate_grid,
  paste0(output_prefix, "-candidate-grid.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    calibration_summary = calibration_summary,
    tuning_diagnostics = tuning_diagnostics,
    selection_summary = selection_summary,
    candidate_grid = candidate_grid,
    raw_diagnostics = raw_diagnostics,
    design = list(
      n = n,
      q_penalized = q_penalized,
      q_total = q_total,
      hard_residual_rank = grid$hard_residual_rank,
      signal_variance = signal_variance,
      signal_pattern = signal_pattern,
      signal_rank = signal_rank,
      confounding_strength = confounding_strength,
      orthogonality_error = geometry$orthogonality_error
    ),
    cross_validation = list(
      criterion = "exact leave-one-out mean squared prediction error",
      minimum_rule = "largest penalty among numerical minimum-MSE ties",
      one_se_rule = paste(
        "largest penalty with mean LOOCV loss no greater than",
        "the minimum plus the standard error of per-observation losses",
        "at the minimizing candidate"
      ),
      locked = paste(
        "select X and Y once on the observed data and hold the observed",
        "X-side candidate fixed for every raw-X permutation"
      )
    ),
    R = R,
    B = B,
    seed = simulation_seed,
    cell_seed = cell_seed,
    nominal_level = nominal_level
  ),
  paste0(output_prefix, "-results.rds")
)

print(calibration_summary, row.names = FALSE, digits = 5)
print(tuning_diagnostics, row.names = FALSE, digits = 5)
message("Saved KRR confirmation outputs with prefix: ", output_prefix)
