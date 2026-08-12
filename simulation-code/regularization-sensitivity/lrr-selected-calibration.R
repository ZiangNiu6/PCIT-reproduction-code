#!/usr/bin/env Rscript

# Confirm the q=178, leading-rank-16, V=7.5 minimum-CV and +1SE LRR
# calibration under observed-selection
# locked tuning. For each null dataset, exact LOOCV selects the X- and Y-side
# penalties once on the observed data. The observed X-side penalty is then
# held fixed while the complete raw X vector is permuted B times and refit.
# No cross-validation or tuning selection is performed on permuted X.

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
default_output_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)

read_positive_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

n <- 300L
q_penalized <- read_positive_integer("PCIT_LRR_Q_PENALIZED", "178")
unpenalized_rank <- 2L
q_total <- q_penalized + unpenalized_rank
hard_residual_rank <- n - q_total
signal_variance <- as.numeric(Sys.getenv(
  "PCIT_LRR_SIGNAL_VARIANCE",
  unset = "7.5"
))
signal_pattern <- Sys.getenv("PCIT_LRR_SIGNAL_PATTERN", unset = "leading")
signal_rank <- as.integer(Sys.getenv(
  "PCIT_LRR_SIGNAL_RANK",
  unset = "16"
))
confounding_strength <- 2
nominal_level <- 0.05
R <- read_positive_integer("PCIT_LRR_R", "1000")
B <- read_positive_integer("PCIT_LRR_B", "499")
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(detected_cores) != 1L ||
    !is.finite(detected_cores) ||
    detected_cores < 1L) {
  detected_cores <- 1L
}
cores <- read_positive_integer(
  "PCIT_LRR_CORES",
  as.character(min(4L, detected_cores))
)
block_size <- read_positive_integer("PCIT_LRR_BLOCK_SIZE", "10")
base_seed <- read_positive_integer("PCIT_LRR_SEED", "20260803")
seed_block_size <- read_positive_integer("PCIT_LRR_SEED_BLOCK_SIZE", "250")
seed_block_increment <- read_positive_integer(
  "PCIT_LRR_SEED_BLOCK_INCREMENT",
  "1000000"
)
output_directory <- Sys.getenv(
  "PCIT_LRR_OUTPUT",
  unset = default_output_directory
)

if (length(signal_variance) != 1L ||
    !is.finite(signal_variance) ||
    signal_variance < 0) {
  stop("PCIT_LRR_SIGNAL_VARIANCE must be one nonnegative number.")
}
if (q_penalized >= n - unpenalized_rank || q_penalized %% 2L != 0L) {
  stop(
    "PCIT_LRR_Q_PENALIZED must be even and smaller than n - ",
    "unpenalized_rank."
  )
}
if (!signal_pattern %in% c("spread", "leading")) {
  stop("PCIT_LRR_SIGNAL_PATTERN must be 'spread' or 'leading'.")
}
if (length(signal_rank) != 1L ||
    is.na(signal_rank) ||
    signal_rank < 2L ||
    signal_rank > q_penalized ||
    signal_rank %% 2L != 0L) {
  stop("PCIT_LRR_SIGNAL_RANK must be even and between 2 and q_penalized.")
}
if (B < 19L) {
  stop("PCIT_LRR_B must be at least 19.")
}

effective_dimension <- function(weights) {
  sum(weights^2)^2 / sum(weights^4)
}

find_penalty <- function(weight_function, target) {
  objective <- function(log_ratio) {
    effective_dimension(weight_function(exp(log_ratio))) - target
  }
  bracket <- c(-30, 30)
  endpoint_values <- vapply(bracket, objective, numeric(1))
  if (prod(sign(endpoint_values)) > 0) {
    stop("Could not attain effective residual dimension ", target, ".")
  }
  exp(stats::uniroot(objective, bracket, tol = 1e-12)$root)
}

make_linear_geometry <- function(n, design_seed = 20260804L) {
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

make_lrr_grid <- function(geometry) {
  residual_rank <- n - unpenalized_rank
  eigenvalues <- if (signal_pattern == "leading") {
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
  penalty_scale <- stats::median(eigenvalues)
  weight_function <- function(ratio) {
    lambda <- ratio * penalty_scale
    c(
      lambda / (eigenvalues + lambda),
      rep(1, hard_residual_rank)
    )
  }

  # Add the q-dependent hard-fit endpoint to a common interpretable path,
  # retaining only dimensions strictly between the two endpoints.
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
  interior_ratios <- vapply(
    interior_dimensions,
    function(target) find_penalty(weight_function, target),
    numeric(1)
  )
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
    stop("The LRR effective-dimension path is incorrect.")
  }
  list(
    basis = geometry$perpendicular,
    active_basis = geometry$perpendicular[
      , seq_len(q_penalized), drop = FALSE
    ],
    active_weights = weights[seq_len(q_penalized), , drop = FALSE],
    residual_diagonal = (geometry$perpendicular^2) %*% weights,
    residual_leverage = (geometry$perpendicular^2) %*% (weights^2),
    candidate_index = seq_along(penalty_ratios),
    effective_dimension = attained_dimensions,
    penalty_ratio = penalty_ratios,
    penalty_absolute = penalty_ratios * penalty_scale,
    eigenvalues = eigenvalues
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

apply_selected_candidates <- function(grid, geometry, indices, values) {
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  if (length(indices) != ncol(values)) {
    stop("Candidate indices do not match the value columns.")
  }
  residual_infinite <- values - geometry$unpenalized %*%
    crossprod(geometry$unpenalized, values)
  active_coordinates <- crossprod(grid$active_basis, values)
  result <- matrix(NA_real_, nrow = nrow(values), ncol = ncol(values))
  for (candidate in sort(unique(indices))) {
    columns <- which(indices == candidate)
    adjustment <- grid$active_weights[, candidate] - 1
    if (all(adjustment == 0)) {
      result[, columns] <- residual_infinite[, columns, drop = FALSE]
    } else {
      result[, columns] <- residual_infinite[, columns, drop = FALSE] +
        grid$active_basis %*% (
          adjustment * active_coordinates[, columns, drop = FALSE]
        )
    }
  }
  result
}

select_exact_loocv <- function(grid, geometry, values) {
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  residual_infinite <- values - geometry$unpenalized %*%
    crossprod(geometry$unpenalized, values)
  active_coordinates <- crossprod(grid$active_basis, values)
  candidate_count <- length(grid$candidate_index)
  score_matrix <- matrix(
    NA_real_, nrow = ncol(values), ncol = candidate_count
  )
  se_matrix <- score_matrix
  minimum_score <- rep(Inf, ncol(values))
  minimum_index <- rep(NA_integer_, ncol(values))
  minimum_residual <- matrix(
    NA_real_, nrow = nrow(values), ncol = ncol(values)
  )

  # Candidate order is increasing penalty/effective dimension. Strict
  # improvement gives which.min behavior (the first candidate on ties).
  for (candidate in grid$candidate_index) {
    adjustment <- grid$active_weights[, candidate] - 1
    residual <- if (all(adjustment == 0)) {
      residual_infinite
    } else {
      residual_infinite + grid$active_basis %*% (
        adjustment * active_coordinates
      )
    }
    loo_residual <- sweep(
      residual,
      1L,
      grid$residual_diagonal[, candidate],
      "/"
    )
    losses <- loo_residual^2
    score <- colMeans(losses)
    squared_loss_sum <- colSums(losses^2)
    loss_variance <- pmax(
      (squared_loss_sum - nrow(losses) * score^2) /
        (nrow(losses) - 1L),
      0
    )
    score_matrix[, candidate] <- score
    se_matrix[, candidate] <- sqrt(loss_variance / nrow(losses))

    replace <- is.na(minimum_index) | score < minimum_score
    if (any(replace)) {
      minimum_score[replace] <- score[replace]
      minimum_index[replace] <- candidate
      minimum_residual[, replace] <- residual[, replace, drop = FALSE]
    }
  }

  row_index <- seq_len(nrow(score_matrix))
  one_se_threshold <- score_matrix[cbind(row_index, minimum_index)] +
    se_matrix[cbind(row_index, minimum_index)]
  one_se_index <- max.col(
    score_matrix <= one_se_threshold,
    ties.method = "last"
  )
  one_se_residual <- apply_selected_candidates(
    grid,
    geometry,
    one_se_index,
    values
  )

  list(
    minimum = list(
      index = minimum_index,
      residual = minimum_residual,
      score = score_matrix[cbind(row_index, minimum_index)],
      se = se_matrix[cbind(row_index, minimum_index)]
    ),
    one_se = list(
      index = one_se_index,
      residual = one_se_residual,
      score = score_matrix[cbind(row_index, one_se_index)],
      threshold = one_se_threshold
    )
  )
}

column_statistics <- function(residual_x, residual_y) {
  scores <- residual_x * residual_y
  numerator <- colSums(scores)
  denominator_squared <- colSums(scores^2) - numerator^2 / nrow(scores)
  numerator / sqrt(pmax(denominator_squared, .Machine$double.xmin))
}

summarize_locked_orbit <- function(
    locked_residual_x,
    residual_y,
    selected_observed_residual_x) {
  statistic <- column_statistics(locked_residual_x, residual_y)
  selected_observed_statistic <- column_statistics(
    matrix(selected_observed_residual_x, ncol = 1L),
    matrix(residual_y, ncol = 1L)
  )[[1L]]
  if (abs(statistic[[1L]] - selected_observed_statistic) > 1e-10) {
    stop("Locked and selected observed statistics do not agree.")
  }
  observed_statistic <- statistic[[1L]]
  reference <- statistic[-1L]
  list(
    observed_statistic = observed_statistic,
    pvalue = (
      1 + sum(abs(reference) >= abs(observed_statistic))
    ) / (B + 1),
    reference_q95 = unname(stats::quantile(
      abs(reference), 0.95, type = 8
    ))
  )
}

geometry <- make_linear_geometry(n)
grid <- make_lrr_grid(geometry)
signals <- make_orthogonal_signals(grid)
mean_x <- 0.5 + confounding_strength * geometry$confounder + signals$x
mean_y <- -0.5 + confounding_strength * geometry$confounder + signals$y

# Draw all X and Y innovations first, followed by the raw-X permutation bank.
innovation_x <- matrix(NA_real_, nrow = n, ncol = R)
innovation_y <- matrix(NA_real_, nrow = n, ncol = R)
permutation_states <- vector("list", R)
seed_blocks <- split(seq_len(R), ceiling(seq_len(R) / seed_block_size))
for (seed_block_index in seq_along(seed_blocks)) {
  indices <- seed_blocks[[seed_block_index]]
  block_seed <- base_seed +
    (seed_block_index - 1L) * seed_block_increment +
    q_penalized * 1009L
  if (block_seed > .Machine$integer.max) {
    stop("The requested seed schedule exceeds R's integer range.")
  }
  set.seed(block_seed)
  innovation_x[, indices] <- matrix(
    stats::rnorm(n * length(indices)),
    nrow = n
  )
  innovation_y[, indices] <- matrix(
    stats::rnorm(n * length(indices)),
    nrow = n
  )
  for (replicate_index in indices) {
    permutation_states[[replicate_index]] <- .Random.seed
    for (draw in seq_len(B)) {
      sample.int(n)
    }
  }
}
blocks <- split(seq_len(R), ceiling(seq_len(R) / block_size))

process_block <- function(replicates) {
  replicate_count <- length(replicates)
  columns_per_orbit <- B + 1L
  observed_x <- mean_x + innovation_x[, replicates, drop = FALSE]
  observed_y <- mean_y + innovation_y[, replicates, drop = FALSE]

  # Select X and Y exactly once on observed data.
  selection <- select_exact_loocv(
    grid,
    geometry,
    cbind(observed_y, observed_x)
  )
  y_columns <- seq_len(replicate_count)
  x_columns <- replicate_count + seq_len(replicate_count)

  # Construct the raw-X permutation orbits without running LOOCV on them.
  x_orbits <- matrix(
    NA_real_,
    nrow = n,
    ncol = replicate_count * columns_per_orbit
  )
  for (local_index in seq_along(replicates)) {
    replicate_index <- replicates[[local_index]]
    x <- observed_x[, local_index]
    destination <- (local_index - 1L) * columns_per_orbit +
      seq_len(columns_per_orbit)
    x_orbits[, destination[[1L]]] <- x
    assign(
      ".Random.seed",
      permutation_states[[replicate_index]],
      envir = .GlobalEnv
    )
    for (draw in seq_len(B)) {
      x_orbits[, destination[[draw + 1L]]] <- x[sample.int(n)]
    }
  }

  rule_objects <- list(
    "Minimum CV" = selection$minimum,
    "+1SE" = selection$one_se
  )
  rule_rows <- lapply(names(rule_objects), function(rule_name) {
    selected <- rule_objects[[rule_name]]
    y_indices <- selected$index[y_columns]
    x_indices <- selected$index[x_columns]
    y_residuals <- selected$residual[, y_columns, drop = FALSE]
    x_residuals <- selected$residual[, x_columns, drop = FALSE]

    locked_indices <- rep(x_indices, each = columns_per_orbit)
    locked_residuals <- apply_selected_candidates(
      grid,
      geometry,
      locked_indices,
      x_orbits
    )

    rows <- vector("list", replicate_count)
    for (local_index in seq_along(replicates)) {
      source <- (local_index - 1L) * columns_per_orbit +
        seq_len(columns_per_orbit)
      orbit_summary <- summarize_locked_orbit(
        locked_residual_x = locked_residuals[, source, drop = FALSE],
        residual_y = y_residuals[, local_index],
        selected_observed_residual_x = x_residuals[, local_index]
      )

      projected_y <- apply_selected_candidates(
        grid,
        geometry,
        x_indices[[local_index]],
        y_residuals[, local_index]
      )[, 1L]
      numerator_variance <- sum(projected_y^2)
      expected_studentizer <- sum(
        grid$residual_leverage[, x_indices[[local_index]]] *
          y_residuals[, local_index]^2
      ) - numerator_variance / n
      if (!is.finite(expected_studentizer) || expected_studentizer <= 0) {
        stop("The selected-operator plug-in studentizer is nonpositive.")
      }
      gamma_plugin <- numerator_variance / expected_studentizer

      rows[[local_index]] <- data.frame(
        selection_rule = rule_name,
        replicate = replicates[[local_index]],
        observed_statistic = orbit_summary$observed_statistic,
        gamma_plugin = gamma_plugin,
        alpha_normal_plugin = 2 * stats::pnorm(
          -stats::qnorm(1 - nominal_level / 2) / sqrt(gamma_plugin)
        ),
        normal_reject = abs(orbit_summary$observed_statistic) >
          stats::qnorm(1 - nominal_level / 2),
        locked_pvalue = orbit_summary$pvalue,
        locked_reject = orbit_summary$pvalue <= nominal_level,
        observed_x_candidate = x_indices[[local_index]],
        observed_y_candidate = y_indices[[local_index]],
        observed_x_effective_dimension = grid$effective_dimension[
          x_indices[[local_index]]
        ],
        observed_y_effective_dimension = grid$effective_dimension[
          y_indices[[local_index]]
        ],
        locked_reference_abs_q95 = orbit_summary$reference_q95,
        stringsAsFactors = FALSE
      )
    }
    list(
      raw = do.call(rbind, rows),
      observed_x_counts = tabulate(
        x_indices,
        nbins = length(grid$candidate_index)
      ),
      observed_y_counts = tabulate(
        y_indices,
        nbins = length(grid$candidate_index)
      )
    )
  })
  names(rule_rows) <- names(rule_objects)
  rule_rows
}

message(
  "Running locked-only paired LRR minimum-CV/+1SE confirmation: n=", n,
  ", q_penalized=", q_penalized,
  ", V=", signal_variance,
  ", R=", R,
  ", B=", B,
  ", blocks=", length(blocks),
  ", cores=", cores,
  ", seed blocks=", length(seed_blocks)
)
if (.Platform$OS.type == "unix" && cores > 1L) {
  block_results <- parallel::mclapply(
    blocks,
    process_block,
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  block_results <- lapply(blocks, process_block)
}

rules <- c("Minimum CV", "+1SE")
raw <- do.call(rbind, lapply(rules, function(rule) {
  do.call(rbind, lapply(block_results, function(block) block[[rule]]$raw))
}))
raw <- raw[order(match(raw$selection_rule, rules), raw$replicate), ]
rownames(raw) <- NULL

exact_binomial_row <- function(rule, calibration, rejected) {
  rejected <- as.logical(rejected)
  rejections <- sum(rejected)
  interval <- stats::binom.test(
    rejections,
    length(rejected),
    conf.level = 0.95
  )$conf.int
  estimate <- mean(rejected)
  data.frame(
    learner = "LRR",
    selection_rule = rule,
    calibration = calibration,
    rejections = rejections,
    R = length(rejected),
    B = if (calibration == "Normal calibration") NA_integer_ else B,
    estimate = estimate,
    mc_se = sqrt(estimate * (1 - estimate) / length(rejected)),
    lower_95_exact = interval[[1L]],
    upper_95_exact = interval[[2L]],
    interval_method = "Clopper-Pearson exact",
    stringsAsFactors = FALSE
  )
}

calibration_rows <- lapply(rules, function(rule) {
  rule_raw <- raw[raw$selection_rule == rule, ]
  rbind(
    exact_binomial_row(
      rule, "Normal calibration", rule_raw$normal_reject
    ),
    exact_binomial_row(
      rule, "Locked-tuning refit", rule_raw$locked_reject
    )
  )
})
calibration_summary <- do.call(rbind, calibration_rows)

tuning_rows <- lapply(rules, function(rule) {
  rule_raw <- raw[raw$selection_rule == rule, ]
  data.frame(
    learner = "LRR",
    selection_rule = rule,
    n = n,
    q_penalized = q_penalized,
    q_total = q_total,
    hard_residual_rank = hard_residual_rank,
    signal_variance = signal_variance,
    signal_pattern = signal_pattern,
    signal_rank = signal_rank,
    R = R,
    B = B,
    mean_observed_x_effective_dimension = mean(
      rule_raw$observed_x_effective_dimension
    ),
    median_observed_x_effective_dimension = stats::median(
      rule_raw$observed_x_effective_dimension
    ),
    mean_observed_y_effective_dimension = mean(
      rule_raw$observed_y_effective_dimension
    ),
    median_observed_y_effective_dimension = stats::median(
      rule_raw$observed_y_effective_dimension
    ),
    observed_x_hard_fit_rate = mean(
      rule_raw$observed_x_effective_dimension == hard_residual_rank
    ),
    observed_y_hard_fit_rate = mean(
      rule_raw$observed_y_effective_dimension == hard_residual_rank
    ),
    observed_statistic_mean = mean(rule_raw$observed_statistic),
    observed_statistic_sd = stats::sd(rule_raw$observed_statistic),
    mean_gamma_plugin = mean(rule_raw$gamma_plugin),
    median_gamma_plugin = stats::median(rule_raw$gamma_plugin),
    mean_alpha_normal_plugin = mean(rule_raw$alpha_normal_plugin),
    observed_statistic_abs_q95 = unname(stats::quantile(
      abs(rule_raw$observed_statistic), 0.95, type = 8
    )),
    mean_locked_reference_abs_q95 = mean(
      rule_raw$locked_reference_abs_q95
    ),
    median_locked_reference_abs_q95 = stats::median(
      rule_raw$locked_reference_abs_q95
    ),
    stringsAsFactors = FALSE
  )
})
tuning_diagnostics <- do.call(rbind, tuning_rows)

selection_rows <- lapply(rules, function(rule) {
  aggregate_counts <- function(name) {
    Reduce(
      `+`,
      lapply(block_results, function(block) block[[rule]][[name]])
    )
  }
  sample_definitions <- list(
    "Observed X" = list(count = aggregate_counts("observed_x_counts")),
    "Observed Y" = list(count = aggregate_counts("observed_y_counts"))
  )
  do.call(rbind, lapply(names(sample_definitions), function(sample_type) {
    values <- sample_definitions[[sample_type]]
    data.frame(
      learner = "LRR",
      selection_rule = rule,
      sample_type = sample_type,
      candidate = grid$candidate_index,
      effective_dimension = grid$effective_dimension,
      penalty_ratio = grid$penalty_ratio,
      penalty_absolute = grid$penalty_absolute,
      count = values$count,
      total = R,
      selection_rate = values$count / R,
      stringsAsFactors = FALSE
    )
  }))
})
selection_distribution <- do.call(rbind, selection_rows)

grid_summary <- data.frame(
  learner = "LRR",
  candidate = grid$candidate_index,
  effective_dimension = grid$effective_dimension,
  penalty_ratio = grid$penalty_ratio,
  penalty_absolute = grid$penalty_absolute,
  stringsAsFactors = FALSE
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  calibration_summary,
  file.path(output_directory, "lrr-calibration-summary.csv"),
  row.names = FALSE
)
write.csv(
  tuning_diagnostics,
  file.path(output_directory, "lrr-tuning-diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  selection_distribution,
  file.path(output_directory, "lrr-selection-distribution.csv"),
  row.names = FALSE
)
write.csv(
  grid_summary,
  file.path(output_directory, "lrr-candidate-grid.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    calibration_summary = calibration_summary,
    tuning_diagnostics = tuning_diagnostics,
    selection_distribution = selection_distribution,
    candidate_grid = grid_summary,
    raw = raw,
    design = list(
      n = n,
      q_penalized = q_penalized,
      q_total = q_total,
      hard_residual_rank = hard_residual_rank,
      signal_variance = signal_variance,
      signal_pattern = signal_pattern,
      signal_rank = signal_rank,
      confounding_strength = confounding_strength,
      R = R,
      B = B,
      base_seed = base_seed,
      seed_block_size = seed_block_size,
      seed_block_increment = seed_block_increment,
      seed_block_base_values = base_seed +
        (seq_along(seed_blocks) - 1L) * seed_block_increment,
      one_se_definition = paste(
        "largest-penalty/largest-effective-dimension candidate with mean",
        "LOOCV loss <= minimum mean + SE at the minimizing candidate"
      ),
      permutation_calibration = paste(
        "select X once on observed data, hold its candidate fixed,",
        "and refit B raw-X permutations without re-running LOOCV"
      )
    )
  ),
  file.path(output_directory, "lrr-confirmation-results.rds"),
  compress = "xz"
)

print(calibration_summary, row.names = FALSE, digits = 5)
print(tuning_diagnostics, row.names = FALSE, digits = 5)
message("Saved locked-only LRR outputs under ", normalizePath(output_directory))
