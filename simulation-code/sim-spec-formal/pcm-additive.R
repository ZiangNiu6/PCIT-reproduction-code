###############################################################################
# Simulation specifier: three-row GCM/PCM continuous-regression comparison
#
# DGPs (Z1 and ε are independent standard Gaussian; ξ is row-specific):
#   setting1_xindep:
#       X = ξ,                         Y = sin(2πZ1) + γ X² + ε
#   setting1 (Additive Gaussian):
#       X = sin(2πZ1) + ξ,            Y = sin(2πZ1) + γ X² + ε
#   setting2 (Heteroscedastic; ξ+1~Exp(1)):
#       X = sin(2πZ1)(1-ξ),          Y = sin(2πZ1) + γ X + ε
#
# Each row has its own gamma=0 sample-size grid.  Power uses n=600 for
# settings 1-2 and n=100 for setting1_xindep.  The two quadratic rows
# demonstrate PCM's advantage, whereas the linear heteroscedastic alternative
# makes GCM powerful.
#
# Matched-calibration method:
#   * One pcmt::perm_gcm call supplies both the observed GCM statistic and its
#     adaptive permutation p-value.  The normal p-value is computed from that
#     same observed statistic.
#   * One pcmt::perm_pcm call supplies both the observed PCM statistic and its
#     adaptive permutation p-value.  The normal p-value is computed from that
#     same observed statistic and therefore uses the same PCM direction and
#     sample split.
#
# Thus, within each test, the normal and permutation procedures differ only in
# calibration.  Random-forest seeds are fixed by nuisance-regression role, so
# each observed/permuted refit uses the same learner randomness.
###############################################################################

#####################
# 1. Parameter grid
#####################
specifier_file <- normalizePath(sys.frame(1L)$ofile)
sys.source(
  file.path(dirname(specifier_file), "pcm-additive-grid.R"),
  envir = environment()
)
parameter_grid <- as.data.frame(parameter_grid)

expected_dgps <- c("setting1_xindep", "setting1", "setting2")
if (nrow(parameter_grid) != 27L ||
    anyDuplicated(parameter_grid[c("n", "dgp", "gamma_signal")]) ||
    !setequal(parameter_grid$dgp, expected_dgps) ||
    any(parameter_grid$ground_truth != (parameter_grid$gamma_signal != 0))) {
  stop(
    "The pcm-additive parameter grid is incompatible with this specifier."
  )
}

#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  B = 1000,          # number of dataset replications per grid point
  seed = 2024
)

###################
# 3. Generate data
###################
generate_data <- function(n, dgp, gamma_signal, B_perm) {
  Z1 <- rnorm(n)
  Z <- matrix(Z1, ncol = 1)
  colnames(Z) <- "Z1"
  eps <- rnorm(n)

  if (dgp == "setting1") {
    xi <- rnorm(n)
    X <- sin(2 * pi * Z1) + xi
    Y <- sin(2 * pi * Z1) + gamma_signal * X^2 + eps
  } else if (dgp == "setting2") {
    xi <- rexp(n, 1) - 1  # ξ+1 ~ Exp(1)
    X <- sin(2 * pi * Z1) - sin(2 * pi * Z1) * xi
    Y <- sin(2 * pi * Z1) + gamma_signal * X + eps
  } else if (dgp == "setting1_xindep") {
    xi <- rnorm(n)
    X <- xi  # X independent of Z
    Y <- sin(2 * pi * Z1) + gamma_signal * X^2 + eps
  } else {
    stop("Unknown pcm-additive DGP: ", dgp)
  }

  list(X = X, Y = Y, Z = Z, n = n, B_perm = B_perm)
}

generate_data_function <- simulatr::simulatr_function(
  f = generate_data,
  arg_names = formalArgs(generate_data),
  loop = TRUE
)

######################
# 4. Method function
######################

matched_gcm_pcm_rf_f <- function(data) {
  schema_version <- "pcm_additive_matched"
  alpha <- 0.05
  adaptive_h <- 15L
  rf_num_trees <- 500L
  rf_num_threads <- 1L
  n_cores <- 1L

  # These purpose-specific constants make the result independent of whether
  # the GCM or PCM code is evaluated first.  The important stability property
  # is that the seed for a nuisance-regression role is reused for every refit
  # of that role along the relevant permutation orbit.
  seeds <- c(
    pcm_split = 53111L,
    gcm_permutation = 53112L,
    pcm_permutation = 53113L,
    gcm_YonZ = 53114L,
    gcm_XonZ = 53115L,
    pcm_YonXZ = 53116L,
    pcm_YonZ = 53117L,
    pcm_YhatonZ = 53118L,
    pcm_RonZ = 53119L
  )

  rf_args <- function(seed) {
    list(
      seed = as.integer(seed),
      num.threads = rf_num_threads,
      num.trees = rf_num_trees
    )
  }

  local_seed <- function(seed, expr) {
    R.utils::withSeed(expr, seed = as.integer(seed))
  }

  assert_probability <- function(x, label) {
    if (length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
      stop(label, " must be one finite number in [0, 1].", call. = FALSE)
    }
    invisible(TRUE)
  }

  validate_fit <- function(fit, label, alternative, B_max) {
    statistic <- as.numeric(fit$statistic)
    p_value <- as.numeric(fit$p.value)
    actual_B <- length(fit$perm_stats)

    if (length(statistic) != 1L || !is.finite(statistic)) {
      stop(label, " returned a nonfinite observed statistic.", call. = FALSE)
    }
    assert_probability(p_value, paste0(label, " permutation p-value"))
    if (!identical(fit$alternative, alternative)) {
      stop(label, " returned an unexpected alternative.", call. = FALSE)
    }
    if (actual_B < 1L || actual_B > B_max ||
        length(fit$parameter) < 1L ||
        as.integer(fit$parameter[["B"]]) != actual_B) {
      stop(label, " returned inconsistent permutation counts.", call. = FALSE)
    }
    if (is.null(fit$adaptive) ||
        !is.data.frame(fit$adaptive$trajectory) ||
        nrow(fit$adaptive$trajectory) != actual_B ||
        !identical(as.integer(fit$adaptive$h), adaptive_h) ||
        !isTRUE(all.equal(as.numeric(fit$adaptive$alpha), alpha)) ||
        !fit$adaptive$stop_reason %in% c("rejection", "futility", "cap") ||
        !identical(isTRUE(fit$adaptive$stopped_early), actual_B < B_max)) {
      stop(label, " returned inconsistent adaptive diagnostics.", call. = FALSE)
    }
    final_adaptive_p <- fit$adaptive$trajectory$p_AV[[actual_B]]
    if (!isTRUE(all.equal(p_value, final_adaptive_p, tolerance = 1e-12))) {
      stop(label, " p-value disagrees with its adaptive trajectory.", call. = FALSE)
    }

    invisible(TRUE)
  }

  B_max <- as.integer(data$B_perm)
  n <- as.integer(data$n)
  if (length(B_max) != 1L || is.na(B_max) || B_max < 1L) {
    stop("B_perm must be one positive integer.", call. = FALSE)
  }
  if (length(n) != 1L || is.na(n) || n < 4L || n != NROW(data$Z) ||
      n != NROW(data$X) || n != NROW(data$Y)) {
    stop("The generated data have inconsistent sample sizes.", call. = FALSE)
  }

  pcm_D2_indices <- local_seed(
    seeds[["pcm_split"]],
    sample.int(n, ceiling(0.5 * n))
  )
  pcm_D2_indices <- as.integer(pcm_D2_indices)
  pcm_D1_indices <- setdiff(seq_len(n), pcm_D2_indices)
  if (length(intersect(pcm_D1_indices, pcm_D2_indices)) != 0L ||
      length(c(pcm_D1_indices, pcm_D2_indices)) != n) {
    stop("The explicit PCM split is not a partition of the observations.",
         call. = FALSE)
  }

  gcm_fit <- local_seed(
    seeds[["gcm_permutation"]],
    pcmt::perm_gcm(
      Y = data$Y,
      X = data$X,
      Z = data$Z,
      reg_YonZ = "rf",
      reg_XonZ = "rf",
      args_YonZ = rf_args(seeds[["gcm_YonZ"]]),
      args_XonZ = rf_args(seeds[["gcm_XonZ"]]),
      B = B_max,
      n_cores = n_cores,
      h = adaptive_h,
      alpha = alpha,
      alternative = "two.sided"
    )
  )

  pcm_fit <- local_seed(
    seeds[["pcm_permutation"]],
    pcmt::perm_pcm(
      Y = data$Y,
      X = data$X,
      Z = data$Z,
      reg_YonXZ = "rf",
      reg_YonZ = "rf",
      reg_YhatonZ = "rf",
      reg_RonZ = "rf",
      args_YonXZ = rf_args(seeds[["pcm_YonXZ"]]),
      args_YonZ = rf_args(seeds[["pcm_YonZ"]]),
      args_YhatonZ = rf_args(seeds[["pcm_YhatonZ"]]),
      args_RonZ = rf_args(seeds[["pcm_RonZ"]]),
      frac = 0.5,
      indices = pcm_D2_indices,
      B = B_max,
      n_cores = n_cores,
      h = adaptive_h,
      alpha = alpha,
      alternative = "greater"
    )
  )

  validate_fit(gcm_fit, "perm_gcm", "two.sided", B_max)
  validate_fit(pcm_fit, "perm_pcm", "greater", B_max)
  if (!identical(as.integer(pcm_fit$D2_indices), pcm_D2_indices) ||
      !identical(as.integer(pcm_fit$D1_indices), pcm_D1_indices)) {
    stop("perm_pcm did not use the supplied sample split.", call. = FALSE)
  }

  statistics <- c(
    gcm = as.numeric(gcm_fit$statistic),
    pcm = as.numeric(pcm_fit$statistic)
  )
  p_values <- c(
    perm_gcm_rf = as.numeric(gcm_fit$p.value),
    asym_gcm_rf = 2 * stats::pnorm(abs(statistics[["gcm"]]),
                                   lower.tail = FALSE),
    perm_pcm_rf = as.numeric(pcm_fit$p.value),
    asym_pcm_rf = stats::pnorm(statistics[["pcm"]], lower.tail = FALSE)
  )
  invisible(lapply(names(p_values), function(name) {
    assert_probability(p_values[[name]], name)
  }))

  adaptive_diagnostic <- function(fit) {
    list(
      actual_B = as.integer(length(fit$perm_stats)),
      K = as.integer(fit$adaptive$K),
      t_star = as.numeric(fit$adaptive$t_star),
      stopped_early = isTRUE(fit$adaptive$stopped_early),
      stop_reason = as.character(fit$adaptive$stop_reason)
    )
  }

  # A compact, deterministic checksum is sufficient here because the supplied
  # indices are also checked against the indices returned by pcmt::perm_pcm.
  D2_checksum <- paste0(
    "n", length(pcm_D2_indices),
    "-s", sum(pcm_D2_indices),
    "-w", sum(seq_along(pcm_D2_indices) * pcm_D2_indices)
  )

  list(
    schema_version = schema_version,
    p_values = p_values,
    statistics = statistics,
    alternatives = c(gcm = "two.sided", pcm = "greater"),
    matched = c(gcm = TRUE, pcm = TRUE, all = TRUE),
    adaptive = list(
      gcm = adaptive_diagnostic(gcm_fit),
      pcm = adaptive_diagnostic(pcm_fit)
    ),
    pcm_split = list(
      D1_n = as.integer(length(pcm_D1_indices)),
      D2_n = as.integer(length(pcm_D2_indices)),
      D2_checksum = D2_checksum
    ),
    seeds = seeds,
    configuration = list(
      B_max = B_max,
      h = adaptive_h,
      alpha = alpha,
      rf_num_trees = rf_num_trees,
      rf_num_threads = rf_num_threads,
      n_cores = n_cores
    ),
    flags = c(
      matched_statistics = TRUE,
      explicit_pcm_split = TRUE,
      fixed_rf_randomness = TRUE,
      single_thread_rf = TRUE
    )
  )
}

matched_gcm_pcm_rf_spec_f <- simulatr::simulatr_function(
  f = matched_gcm_pcm_rf_f,
  arg_names = character(0),
  loop = TRUE
)

run_method_functions <- list(
  matched_gcm_pcm_rf = matched_gcm_pcm_rf_spec_f
)

############################
# 5. Assemble specifier
############################
simulatr_spec <- simulatr::simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function,
  run_method_functions
)

# ############################### check results ##################################
# check_results_nofit <- simulatr::check_simulatr_specifier_object(simulatr_spec, B_in = 1)
