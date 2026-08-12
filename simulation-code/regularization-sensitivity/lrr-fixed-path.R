#!/usr/bin/env Rscript

# Locked-only fixed-lambda LRR path experiment. Each candidate is fixed before
# seeing X or Y; no cross-validation occurs. All candidates share the same
# innovations and raw-X permutation indices within each Monte Carlo dataset.

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

read_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

read_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || !is.finite(value)) {
    stop(name, " must be one finite number.")
  }
  value
}

started_at <- proc.time()[["elapsed"]]
n <- 300L
unpenalized_rank <- 2L
q_penalized <- 178L
q_total <- q_penalized + unpenalized_rank
hard_residual_rank <- n - q_total
signal_rank <- read_int("PCIT_FIXED_SIGNAL_RANK", "16")
signal_variance <- read_num("PCIT_FIXED_V", "7.5")
R <- read_int("PCIT_FIXED_R", "1000")
B <- read_int("PCIT_FIXED_B", "499")
cores <- read_int("PCIT_FIXED_CORES", "4")
block_size <- read_int("PCIT_FIXED_BLOCK_SIZE", "5")
base_seed <- read_int("PCIT_FIXED_SEED", "20260803")
seed_block_size <- read_int("PCIT_FIXED_SEED_BLOCK_SIZE", "250")
seed_block_increment <- read_int(
  "PCIT_FIXED_SEED_BLOCK_INCREMENT",
  "1000000"
)
output_csv <- Sys.getenv(
  "PCIT_FIXED_OUTPUT",
  unset = file.path(regularization_output_directory, "lrr-fixed-path-summary.csv")
)
output_rds <- sub("[.]csv$", ".rds", output_csv)
if (identical(output_rds, output_csv)) {
  output_rds <- paste0(output_csv, ".rds")
}
nominal_level <- 0.05
critical_value <- stats::qnorm(1 - nominal_level / 2)

if (q_penalized %% 2L != 0L || q_penalized >= n - unpenalized_rank) {
  stop("The fixed q_penalized must be even and smaller than n - 2.")
}
if (signal_rank %% 2L != 0L || signal_rank > q_penalized) {
  stop("PCIT_FIXED_SIGNAL_RANK must be even and no larger than q.")
}
if (signal_variance < 0) stop("PCIT_FIXED_V must be nonnegative.")
if (B < 19L) stop("PCIT_FIXED_B must be at least 19.")

alpha_normal <- function(gamma) {
  2 * stats::pnorm(-critical_value / sqrt(gamma))
}

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
    stop("Cannot attain effective dimension ", target, ".")
  }
  exp(stats::uniroot(objective, endpoints, tol = 1e-12)$root)
}

make_linear_geometry <- function(design_seed = 20260804L) {
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
  orthogonality_error <- max(abs(
    crossprod(cbind(unpenalized, perpendicular)) - diag(n)
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

make_grid <- function(geometry) {
  residual_rank <- n - unpenalized_rank
  eigenvalues <- exp(seq(
    log(4 * n),
    log(0.25 * n),
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

  # q=178 gives exactly 20 path points: 120,130,...,290,295,298.
  target_dimensions <- c(
    seq.int(hard_residual_rank, 290L, by = 10L),
    295L,
    residual_rank
  )
  interior <- target_dimensions[
    target_dimensions > hard_residual_rank &
      target_dimensions < residual_rank
  ]
  ratios <- c(
    0,
    vapply(interior, function(target) {
      find_penalty(weight_function, target)
    }, numeric(1)),
    Inf
  )
  weights <- do.call(cbind, lapply(ratios, function(ratio) {
    if (ratio == 0) {
      c(rep(0, q_penalized), rep(1, hard_residual_rank))
    } else if (is.infinite(ratio)) {
      rep(1, residual_rank)
    } else {
      weight_function(ratio)
    }
  }))
  attained <- apply(weights, 2L, effective_dimension)
  if (length(attained) != 20L ||
      max(abs(attained - target_dimensions)) > 1e-6) {
    stop("The requested fine 20-candidate path was not attained.")
  }
  list(
    targets = attained,
    ratios = ratios,
    penalty_absolute = ratios * penalty_scale,
    weights = weights,
    active_basis = geometry$perpendicular[
      , seq_len(q_penalized), drop = FALSE
    ],
    active_weights = weights[seq_len(q_penalized), , drop = FALSE],
    A2_diagonal = (geometry$perpendicular^2) %*% (weights^2),
    eigenvalues = eigenvalues
  )
}

make_signals <- function(geometry) {
  active_x <- seq.int(1L, signal_rank, by = 2L)
  active_y <- seq.int(2L, signal_rank, by = 2L)
  scale <- sqrt(n * signal_variance)
  list(
    x = as.vector(
      scale * geometry$perpendicular[, active_x, drop = FALSE] %*%
        rep(1 / sqrt(length(active_x)), length(active_x))
    ),
    y = as.vector(
      scale * geometry$perpendicular[, active_y, drop = FALSE] %*%
        rep(1 / sqrt(length(active_y)), length(active_y))
    )
  )
}

exact_interval <- function(rejections, trials) {
  unname(stats::binom.test(
    rejections,
    trials,
    conf.level = 0.95
  )$conf.int)
}

geometry <- make_linear_geometry()
grid <- make_grid(geometry)
signals <- make_signals(geometry)
mean_x <- 0.5 + 2 * geometry$confounder + signals$x
mean_y <- -0.5 + 2 * geometry$confounder + signals$y

# Shared innovations for every candidate, using the legacy q-indexed blocks.
innovation_x <- matrix(NA_real_, nrow = n, ncol = R)
innovation_y <- matrix(NA_real_, nrow = n, ncol = R)
seed_blocks <- split(seq_len(R), ceiling(seq_len(R) / seed_block_size))
for (seed_block_index in seq_along(seed_blocks)) {
  indices <- seed_blocks[[seed_block_index]]
  block_seed <- base_seed +
    (seed_block_index - 1L) * seed_block_increment +
    q_penalized * 1009L
  set.seed(block_seed)
  innovation_x[, indices] <- matrix(
    stats::rnorm(n * length(indices)),
    nrow = n
  )
  innovation_y[, indices] <- matrix(
    stats::rnorm(n * length(indices)),
    nrow = n
  )
}

# Each replicate gets one deterministic permutation bank shared by all 20
# candidates. This avoids generating or storing candidate-specific banks.
permutation_seed <- base_seed + 500000000L +
  q_penalized * 1009L + seq_len(R) * 7919L
if (max(permutation_seed) > .Machine$integer.max) {
  stop("Permutation seed schedule exceeds R's integer range.")
}

population_rows <- lapply(seq_along(grid$targets), function(candidate) {
  weights <- grid$weights[, candidate]
  active_weights <- grid$active_weights[, candidate]
  residual_mean_y <- signals$y + grid$active_basis %*% (
    (active_weights - 1) * crossprod(grid$active_basis, signals$y)
  )
  projected_mean_y <- residual_mean_y + grid$active_basis %*% (
    (active_weights - 1) * crossprod(
      grid$active_basis,
      residual_mean_y
    )
  )
  numerator_expectation <- sum(projected_mean_y^2) + sum(weights^4)
  residual_y_second_moment <- residual_mean_y^2 +
    grid$A2_diagonal[, candidate]
  studentizer_expectation <- sum(
    grid$A2_diagonal[, candidate] * residual_y_second_moment
  ) - numerator_expectation / n
  gamma_population <- numerator_expectation / studentizer_expectation
  data.frame(
    candidate = candidate,
    gamma_population = gamma_population,
    normal_annihilation_population = alpha_normal(gamma_population)
  )
})
population_theory <- do.call(rbind, population_rows)

blocks <- split(seq_len(R), ceiling(seq_len(R) / block_size))

process_block <- function(replicates) {
  replicate_count <- length(replicates)
  orbit_size <- B + 1L
  observed_x <- mean_x + innovation_x[, replicates, drop = FALSE]
  observed_y <- mean_y + innovation_y[, replicates, drop = FALSE]
  x_orbits <- matrix(
    NA_real_,
    nrow = n,
    ncol = replicate_count * orbit_size
  )
  for (local_index in seq_along(replicates)) {
    replicate <- replicates[[local_index]]
    x <- observed_x[, local_index]
    destination <- (local_index - 1L) * orbit_size + seq_len(orbit_size)
    x_orbits[, destination[[1L]]] <- x
    set.seed(permutation_seed[[replicate]])
    for (draw in seq_len(B)) {
      x_orbits[, destination[[draw + 1L]]] <- x[sample.int(n)]
    }
  }

  residual_infinite_x <- x_orbits - geometry$unpenalized %*%
    crossprod(geometry$unpenalized, x_orbits)
  active_coordinates_x <- crossprod(grid$active_basis, x_orbits)
  residual_infinite_y <- observed_y - geometry$unpenalized %*%
    crossprod(geometry$unpenalized, observed_y)
  active_coordinates_y <- crossprod(grid$active_basis, observed_y)

  candidate_rows <- vector("list", length(grid$targets))
  for (candidate in seq_along(grid$targets)) {
    adjustment <- grid$active_weights[, candidate] - 1
    residual_x <- residual_infinite_x + grid$active_basis %*% (
      adjustment * active_coordinates_x
    )
    residual_y <- residual_infinite_y + grid$active_basis %*% (
      adjustment * active_coordinates_y
    )
    repeated_y <- residual_y[
      , rep(seq_len(replicate_count), each = orbit_size),
      drop = FALSE
    ]
    scores <- residual_x * repeated_y
    numerator <- colSums(scores)
    quadratic <- colSums(scores^2) - numerator^2 / n
    statistic <- numerator / sqrt(pmax(
      quadratic,
      .Machine$double.xmin
    ))
    statistic <- matrix(
      statistic,
      nrow = orbit_size,
      ncol = replicate_count
    )
    observed_statistic <- statistic[1L, ]
    reference <- statistic[-1L, , drop = FALSE]
    locked_pvalue <- (
      1 + colSums(sweep(
        abs(reference),
        2L,
        abs(observed_statistic),
        ">="
      ))
    ) / orbit_size
    locked_q95 <- apply(
      abs(reference),
      2L,
      stats::quantile,
      probs = 0.95,
      type = 8,
      names = FALSE
    )

    # Conditional reference dispersion. It is also the observed-theory gamma
    # under annihilation, but not when A mu_x is nonzero.
    active_coordinates_residual_y <- crossprod(
      grid$active_basis,
      residual_y
    )
    projected_y <- residual_y + grid$active_basis %*% (
      adjustment * active_coordinates_residual_y
    )
    numerator_variance <- colSums(projected_y^2)
    expected_studentizer <- colSums(
      grid$A2_diagonal[, candidate] * residual_y^2
    ) - numerator_variance / n
    if (any(!is.finite(expected_studentizer)) ||
        any(expected_studentizer <= 0)) {
      stop("A fixed-operator plug-in studentizer is nonpositive.")
    }
    gamma_plugin <- numerator_variance / expected_studentizer

    # DGP-specific conditional Gaussian approximation. Unlike the theorem
    # plug-in above, this includes the nonannihilated X-side conditional mean.
    residual_mean_x <- as.vector(
      signals$x + grid$active_basis %*% (
        adjustment * crossprod(grid$active_basis, signals$x)
      )
    )
    conditional_mean_scores <- residual_mean_x * residual_y
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

    candidate_rows[[candidate]] <- data.frame(
      candidate = candidate,
      replicate = replicates,
      observed_statistic = observed_statistic,
      normal_reject = abs(observed_statistic) > critical_value,
      locked_pvalue = locked_pvalue,
      locked_reject = locked_pvalue <= nominal_level,
      locked_reference_abs_q95 = locked_q95,
      gamma_plugin = gamma_plugin,
      alpha_normal_plugin = alpha_normal(gamma_plugin),
      normal_theory_dgp = normal_theory_dgp,
      permutation_theory_dgp = permutation_theory_dgp,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, candidate_rows)
}

message(
  "Running locked fixed-LRR path: q=", q_penalized,
  ", rank=", signal_rank,
  ", V=", signal_variance,
  ", candidates=", length(grid$targets),
  ", R=", R,
  ", B=", B,
  ", blocks=", length(blocks),
  ", cores=", cores
)
block_results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    blocks,
    process_block,
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  lapply(blocks, process_block)
}
raw <- do.call(rbind, block_results)
raw <- raw[order(raw$candidate, raw$replicate), ]
rownames(raw) <- NULL

summary_rows <- lapply(seq_along(grid$targets), function(candidate) {
  cell <- raw[raw$candidate == candidate, ]
  normal_rejections <- sum(cell$normal_reject)
  locked_rejections <- sum(cell$locked_reject)
  normal_interval <- exact_interval(normal_rejections, R)
  locked_interval <- exact_interval(locked_rejections, R)
  population <- population_theory[
    population_theory$candidate == candidate,
  ]
  data.frame(
    learner = "LRR",
    tuning_regime = "fixed lambda",
    n = n,
    R = R,
    B = B,
    q_penalized = q_penalized,
    q_total = q_total,
    hard_residual_rank = hard_residual_rank,
    signal_pattern = "leading",
    signal_rank = signal_rank,
    signal_variance = signal_variance,
    candidate = candidate,
    effective_dimension = grid$targets[[candidate]],
    penalty_ratio = grid$ratios[[candidate]],
    penalty_absolute = grid$penalty_absolute[[candidate]],
    normal_rejections = normal_rejections,
    normal_empirical = normal_rejections / R,
    normal_mc_se = sqrt(
      (normal_rejections / R) * (1 - normal_rejections / R) / R
    ),
    normal_lower_95_exact = normal_interval[[1L]],
    normal_upper_95_exact = normal_interval[[2L]],
    mean_gamma_reference_plugin = mean(cell$gamma_plugin),
    median_gamma_reference_plugin = stats::median(cell$gamma_plugin),
    normal_annihilation_plugin = mean(cell$alpha_normal_plugin),
    normal_theory_dgp = mean(cell$normal_theory_dgp),
    gamma_reference_population = population$gamma_population,
    normal_annihilation_population =
      population$normal_annihilation_population,
    permutation_rejections = locked_rejections,
    permutation_empirical = locked_rejections / R,
    permutation_mc_se = sqrt(
      (locked_rejections / R) * (1 - locked_rejections / R) / R
    ),
    permutation_lower_95_exact = locked_interval[[1L]],
    permutation_upper_95_exact = locked_interval[[2L]],
    nominal_level = nominal_level,
    permutation_theory_dgp = mean(cell$permutation_theory_dgp),
    observed_statistic_mean = mean(cell$observed_statistic),
    observed_statistic_sd = stats::sd(cell$observed_statistic),
    observed_statistic_abs_q95 = unname(stats::quantile(
      abs(cell$observed_statistic), 0.95, type = 8
    )),
    mean_locked_reference_abs_q95 = mean(
      cell$locked_reference_abs_q95
    ),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
runtime_seconds <- proc.time()[["elapsed"]] - started_at
summary$runtime_seconds <- runtime_seconds

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, output_csv, row.names = FALSE)
saveRDS(
  list(
    summary = summary,
    raw = raw,
    population_theory = population_theory,
    candidate_grid = data.frame(
      candidate = seq_along(grid$targets),
      effective_dimension = grid$targets,
      penalty_ratio = grid$ratios,
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
      innovation = "independent N(0,1)",
      null = "orthogonal odd/even leading Fourier signals",
      permutation = "complete raw X; fixed candidate reapplied",
      R = R,
      B = B,
      base_seed = base_seed,
      permutation_seed_formula = paste(
        "base_seed + 500000000 + q*1009 + replicate*7919"
      ),
      orthogonality_error = geometry$orthogonality_error
    ),
    runtime_seconds = runtime_seconds
  ),
  output_rds,
  compress = "xz"
)

print(summary[, c(
  "candidate",
  "effective_dimension",
  "normal_empirical",
  "normal_theory_dgp",
  "permutation_empirical",
  "permutation_theory_dgp"
)], row.names = FALSE, digits = 5)
message("Runtime: ", round(runtime_seconds, 2), " seconds")
message("Saved CSV: ", normalizePath(output_csv))
message("Saved RDS: ", normalizePath(output_rds))
