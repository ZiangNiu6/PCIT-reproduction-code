###############################################################################
# Simulation specifier: five score calibrations under sparse exposure and
# variance misspecification
#
# DGP (n = 400 in the parameter grid):
#   log(mu) = log(mu0) + beta_z Z + gamma X
#
# Focused rows (n = 400):
#   1. common X: Z ~ N(0,1), Pr(X=1)=0.50, correctly specified NB with
#      dispersion 2 (theta 0.5);
#   2. rare X: the same model with Pr(X=1)=0.025;
#   3. severe variance misspecification: Z ~ Bernoulli(0.75), marginal
#      Pr(X=1)=0.14, means 20 and 30 in the two strata, true variance factors
#      Var(Y|Z=z)/E(Y|Z=z) = 5 and 1, and a Poisson working score model.
# Under MI, X is independent of Z. Under CI, X depends on Z while X and Y are
# generated independently conditional on Z at gamma = 0.
#
# Every design/scenario cell has five total signal points: the null and four
# equally spaced downregulation alternatives.
#
# All methods share one null NB fit and the same efficient-score contributions.
# They differ only in calibration:
#   1. standard parametric score with a normal cutoff;
#   2. Barry et al. working-information-studentized permutation;
#   3. empirically studentized raw-label permutation;
#   4. HGF Monte Carlo Rademacher sign flipping (no SPA);
#   5. De Santis standardized sign flipping.
###############################################################################

#####################
# 1. Parameter grid
#####################
specifier_file <- normalizePath(sys.frame(1L)$ofile)
sys.source(
  file.path(dirname(specifier_file), "nb-permutation-signflip-grid.R"),
  envir = environment()
)
parameter_grid <- as.data.frame(parameter_grid)

scenario_levels <- c("common_nb", "rare_nb", "severe_variance")
expected_signal_magnitudes <- list(
  MI = list(
    common_nb = seq(0, 0.45, length.out = 5L),
    rare_nb = seq(0, 0.85, length.out = 5L),
    severe_variance = seq(0, 0.105, length.out = 5L)
  ),
  CI = list(
    common_nb = seq(0, 0.45, length.out = 5L),
    rare_nb = seq(0, 0.85, length.out = 5L),
    severe_variance = seq(0, 0.18, length.out = 5L)
  )
)
if (nrow(parameter_grid) != 30L ||
    !setequal(unique(parameter_grid$null_design), c("MI", "CI")) ||
    !setequal(unique(parameter_grid$scenario_id), scenario_levels) ||
    any(parameter_grid$n != 400L) ||
    any(parameter_grid$nuisance_mode != "fitted") ||
    any(parameter_grid$gamma_signal > 0)) {
  stop("nb-permutation-signflip requires the focused 30-cell parameter grid")
}
for (null_design_value in c("MI", "CI")) {
  for (scenario_value in scenario_levels) {
    group <- parameter_grid[
      parameter_grid$null_design == null_design_value &
        parameter_grid$scenario_id == scenario_value,
      ,
      drop = FALSE
    ]
    expected <- expected_signal_magnitudes[[null_design_value]][[scenario_value]]
    if (nrow(group) != 5L ||
        sum(!group$ground_truth) != 1L ||
        sum(group$ground_truth) != 4L ||
        !isTRUE(all.equal(
          sort(abs(group$gamma_signal)), expected, tolerance = 1e-12
        ))) {
      stop("Unexpected null/alternative signal grid")
    }
  }
}

#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  B = 2000L,
  seed = 2024L,
  n_processors = 1L
)

###############################################################################
# 3. Shared implementation environment
###############################################################################
.nb_signflip_env <- new.env(parent = baseenv())

local({
  .failure_message <- function(prefix, error) {
    paste0(prefix, ": ", conditionMessage(error))
  }

  .fit_null_model <- function(Y, Z, nuisance_mode, working_family) {
    if (!identical(nuisance_mode, "fitted")) {
      return(list(
        success = FALSE,
        failure_reason = paste0("unknown nuisance_mode: ", nuisance_mode),
        mu_hat = rep(NA_real_, length(Y)),
        theta_hat = NA_real_
      ))
    }

    fit <- tryCatch({
      if (identical(working_family, "negative_binomial")) {
        suppressWarnings(
          MASS::glm.nb(
            Y ~ Z,
            control = stats::glm.control(maxit = 100L)
          )
        )
      } else if (identical(working_family, "poisson")) {
        design <- cbind(`(Intercept)` = 1, Z = Z)
        suppressWarnings(
          stats::glm.fit(
            x = design,
            y = Y,
            family = stats::poisson(link = "log"),
            control = stats::glm.control(maxit = 100L)
          )
        )
      } else {
        stop("unknown working_family: ", working_family)
      }
    }, error = function(e) e)

    if (inherits(fit, "error")) {
      return(list(
        success = FALSE,
        failure_reason = .failure_message("null working-model fit failed", fit),
        mu_hat = rep(NA_real_, length(Y)),
        theta_hat = NA_real_
      ))
    }

    mu_hat <- as.numeric(stats::fitted(fit))
    theta_hat <- if (identical(working_family, "negative_binomial")) {
      as.numeric(fit$theta)
    } else {
      NA_real_
    }
    if (identical(working_family, "negative_binomial") &&
        !is.null(fit$th.warn)) {
      return(list(
        success = FALSE,
        failure_reason = paste0(
          "null NB dispersion fit warning: ",
          paste(fit$th.warn, collapse = "; ")
        ),
        mu_hat = mu_hat,
        theta_hat = theta_hat
      ))
    }

    fit_valid <- isTRUE(fit$converged) &&
      length(mu_hat) == length(Y) &&
      all(is.finite(mu_hat)) && all(mu_hat > 0) &&
      (
        identical(working_family, "poisson") ||
          (length(theta_hat) == 1L && is.finite(theta_hat) && theta_hat > 0)
      )

    if (!fit_valid) {
      return(list(
        success = FALSE,
        failure_reason = "null working-model fit returned invalid estimates",
        mu_hat = mu_hat,
        theta_hat = theta_hat
      ))
    }

    list(
      success = TRUE,
      failure_reason = NA_character_,
      mu_hat = mu_hat,
      theta_hat = theta_hat
    )
  }

  .prepare_score_cache <- function(
      Y, X, Z, mu_hat, theta_hat, working_family) {
    n <- length(Y)
    D <- cbind(`(Intercept)` = 1, Z = Z)

    H <- if (identical(working_family, "negative_binomial")) {
      1 / (1 + mu_hat / theta_hat)
    } else if (identical(working_family, "poisson")) {
      rep(1, n)
    } else {
      stop("unknown working_family: ", working_family)
    }
    W <- mu_hat * H
    if (any(!is.finite(H)) || any(!is.finite(W)) || any(W <= 0)) {
      stop("invalid NB score weights")
    }

    weighted_gram <- crossprod(D, D * W)
    gram_inverse <- tryCatch(
      solve(weighted_gram),
      error = function(e) e
    )
    if (inherits(gram_inverse, "error")) {
      stop(.failure_message("weighted nuisance projection failed", gram_inverse))
    }

    projection_coef <- gram_inverse %*% crossprod(D, W * X)
    X_tilde <- as.numeric(X - D %*% projection_coef)
    outcome_factor <- H * (Y - mu_hat)
    score_contrib <- X_tilde * outcome_factor

    sqrt_W <- sqrt(W)
    weighted_design <- D * sqrt_W
    design_qr <- qr(weighted_design)
    if (design_qr$rank < ncol(D)) {
      stop("weighted nuisance design is rank deficient")
    }
    U <- qr.Q(design_qr, complete = FALSE)[, seq_len(ncol(D)), drop = FALSE]

    # b = (I - U U') sqrt(W) X = sqrt(W) X_tilde.  Using X_tilde
    # keeps the efficient-score contribution exactly common across methods.
    b <- sqrt_W * X_tilde
    b_norm_sq <- sum(b^2)
    score_norm_sq <- sum(score_contrib^2)
    degenerate <- !is.finite(b_norm_sq) || !is.finite(score_norm_sq) ||
      b_norm_sq <= .Machine$double.eps || score_norm_sq <= .Machine$double.eps

    if (degenerate) {
      A <- matrix(0, nrow = n, ncol = ncol(D))
      score_eff_n <- 0
      max_score_share <- 0
    } else {
      A <- U * as.numeric(b / sqrt(b_norm_sq))
      score_fourth_sum <- sum(score_contrib^4)
      score_eff_n <- if (score_fourth_sum > 0 && is.finite(score_fourth_sum)) {
        score_norm_sq^2 / score_fourth_sum
      } else {
        0
      }
      max_score_share <- max(score_contrib^2) / score_norm_sq
    }

    list(
      D = D,
      W = W,
      gram_inverse = gram_inverse,
      outcome_factor = outcome_factor,
      X_tilde = X_tilde,
      score_contrib = score_contrib,
      A = A,
      information_norm_sq = as.numeric(b_norm_sq),
      degenerate = degenerate,
      score_eff_n = as.numeric(score_eff_n),
      max_score_share = as.numeric(max_score_share)
    )
  }

  generate_data <- function(
      n, null_design, scenario_id, nuisance_mode,
      z_distribution, z_probability,
      x_prevalence, x_intercept, x_z_slope,
      x_probability_z0, x_probability_z1,
      mu0, beta_z, variance_model, working_family,
      nb_size, variance_factor_z0, variance_factor_z1,
      gamma_signal, B_resample) {
    if (!null_design %in% c("MI", "CI")) {
      stop("null_design must be MI or CI")
    }

    Z <- if (identical(z_distribution, "normal")) {
      stats::rnorm(n)
    } else if (identical(z_distribution, "binary")) {
      stats::rbinom(n, size = 1L, prob = z_probability)
    } else {
      stop("unknown z_distribution: ", z_distribution)
    }
    propensity <- if (identical(z_distribution, "normal")) {
      stats::plogis(x_intercept + x_z_slope * Z)
    } else {
      ifelse(Z == 0L, x_probability_z0, x_probability_z1)
    }
    X <- stats::rbinom(n, size = 1L, prob = propensity)

    mu_null <- mu0 * exp(beta_z * Z)
    mu_true <- mu_null * exp(gamma_signal * X)
    Y <- if (identical(variance_model, "constant_nb")) {
      as.integer(MASS::rnegbin(n = n, mu = mu_true, theta = nb_size))
    } else if (identical(variance_model, "stratum_quasi_poisson")) {
      variance_factor <- ifelse(
        Z == 0L,
        variance_factor_z0,
        variance_factor_z1
      )
      poisson <- abs(variance_factor - 1) <= 1e-12
      output <- integer(n)
      if (any(poisson)) {
        output[poisson] <- stats::rpois(
          sum(poisson), lambda = mu_true[poisson]
        )
      }
      if (any(!poisson)) {
        if (any(variance_factor[!poisson] <= 1)) {
          stop("non-Poisson variance factors must exceed one")
        }
        observation_size <-
          mu_true[!poisson] / (variance_factor[!poisson] - 1)
        output[!poisson] <- stats::rnbinom(
          sum(!poisson),
          mu = mu_true[!poisson],
          size = observation_size
        )
      }
      output
    } else {
      stop("unknown variance_model: ", variance_model)
    }

    # The simulator controls the data seed.  Storing method-specific seeds makes
    # resampling reproducible under any scheduling order.  The two flip methods
    # deliberately share flip_seed and therefore the same sign matrix.
    perm_seed <- sample.int(.Machine$integer.max, size = 1L)
    flip_seed <- sample.int(.Machine$integer.max, size = 1L)

    nuisance <- .fit_null_model(
      Y = Y,
      Z = Z,
      nuisance_mode = nuisance_mode,
      working_family = working_family
    )

    score_cache <- NULL
    success <- nuisance$success
    failure_reason <- nuisance$failure_reason
    if (success) {
      score_cache <- tryCatch(
        .prepare_score_cache(
          Y = Y,
          X = X,
          Z = Z,
          mu_hat = nuisance$mu_hat,
          theta_hat = nuisance$theta_hat,
          working_family = working_family
        ),
        error = function(e) e
      )
      if (inherits(score_cache, "error")) {
        success <- FALSE
        failure_reason <- .failure_message("score preparation failed", score_cache)
        score_cache <- NULL
      }
    }

    list(
      Y = Y,
      X = X,
      Z = Z,
      n = as.integer(n),
      null_design = null_design,
      scenario_id = scenario_id,
      nuisance_mode = nuisance_mode,
      working_family = working_family,
      x_prevalence = x_prevalence,
      propensity = propensity,
      mu_true = mu_true,
      mu_hat = nuisance$mu_hat,
      theta_hat = nuisance$theta_hat,
      B_resample = as.integer(B_resample),
      perm_seed = as.integer(perm_seed),
      flip_seed = as.integer(flip_seed),
      success = success,
      failure_reason = failure_reason,
      score_cache = score_cache,
      n_exposed = as.integer(sum(X == 1L)),
      n_exposed_nonzero = as.integer(sum(X == 1L & Y > 0L))
    )
  }

  .method_output <- function(
      data, p_value, statistic, success = TRUE,
      failure_reason = NA_character_, degenerate = FALSE) {
    cache <- data$score_cache
    list(
      p.value = as.numeric(p_value),
      statistic = as.numeric(statistic),
      success = as.logical(success),
      failure_reason = as.character(failure_reason),
      degenerate = as.logical(degenerate),
      theta_hat = as.numeric(data$theta_hat),
      n_exposed = as.integer(data$n_exposed),
      n_exposed_nonzero = as.integer(data$n_exposed_nonzero),
      score_eff_n = if (is.null(cache)) NA_real_ else as.numeric(cache$score_eff_n),
      max_score_share = if (is.null(cache)) NA_real_ else as.numeric(cache$max_score_share)
    )
  }

  .failed_method_output <- function(data) {
    .method_output(
      data = data,
      p_value = NA_real_,
      statistic = NA_real_,
      success = FALSE,
      failure_reason = data$failure_reason,
      degenerate = FALSE
    )
  }

  .degenerate_method_output <- function(data) {
    .method_output(
      data = data,
      p_value = 1,
      statistic = 0,
      success = TRUE,
      failure_reason = NA_character_,
      degenerate = TRUE
    )
  }

  .studentized_statistic <- function(score_matrix) {
    if (is.null(dim(score_matrix))) {
      score_matrix <- matrix(score_matrix, ncol = 1L)
    }
    n <- nrow(score_matrix)
    numerator <- colSums(score_matrix)
    centered_ss <- colSums(score_matrix^2) - numerator^2 / n
    scale <- pmax(colSums(score_matrix^2), 1)
    positive <- centered_ss > 100 * .Machine$double.eps * scale

    statistic <- numeric(length(numerator))
    statistic[positive] <- numerator[positive] / sqrt(centered_ss[positive])
    zero_denom_nonzero_num <- !positive &
      abs(numerator) > sqrt(.Machine$double.eps * scale)
    statistic[zero_denom_nonzero_num] <-
      sign(numerator[zero_denom_nonzero_num]) * Inf
    statistic
  }

  .information_statistic <- function(score_matrix, x_tilde_matrix, W) {
    if (is.null(dim(score_matrix))) {
      score_matrix <- matrix(score_matrix, ncol = 1L)
    }
    if (is.null(dim(x_tilde_matrix))) {
      x_tilde_matrix <- matrix(x_tilde_matrix, ncol = 1L)
    }
    numerator <- colSums(score_matrix)
    denominator_sq <- colSums(x_tilde_matrix^2 * W)
    scale <- pmax(colSums(x_tilde_matrix^2), 1)
    positive <- denominator_sq > 100 * .Machine$double.eps * scale

    statistic <- numeric(length(numerator))
    statistic[positive] <- numerator[positive] / sqrt(denominator_sq[positive])
    zero_denom_nonzero_num <- !positive &
      abs(numerator) > sqrt(.Machine$double.eps * scale)
    statistic[zero_denom_nonzero_num] <-
      sign(numerator[zero_denom_nonzero_num]) * Inf
    statistic
  }

  parametric_score_f <- function(data) {
    if (!isTRUE(data$success)) {
      return(.failed_method_output(data))
    }
    if (isTRUE(data$score_cache$degenerate)) {
      return(.degenerate_method_output(data))
    }

    cache <- data$score_cache
    observed <- .information_statistic(
      cache$score_contrib,
      cache$X_tilde,
      cache$W
    )[1L]
    p_value <- 2 * stats::pnorm(-abs(observed))
    .method_output(data, p_value = p_value, statistic = observed)
  }

  information_studentized_permutation_f <- function(data) {
    if (!isTRUE(data$success)) {
      return(.failed_method_output(data))
    }
    if (isTRUE(data$score_cache$degenerate)) {
      return(.degenerate_method_output(data))
    }

    cache <- data$score_cache
    observed <- .information_statistic(
      cache$score_contrib,
      cache$X_tilde,
      cache$W
    )[1L]
    B <- data$B_resample
    n <- data$n
    batch_size <- 128L
    hits <- 0L
    completed <- 0L

    # Use exactly the same permutations as the empirical-studentized method so
    # their Monte Carlo comparison is paired dataset by dataset.
    set.seed(data$perm_seed)
    while (completed < B) {
      current_batch <- min(batch_size, B - completed)
      permutation_indices <- vapply(
        seq_len(current_batch),
        function(i) sample.int(n),
        FUN.VALUE = integer(n)
      )
      X_permuted <- matrix(
        data$X[permutation_indices],
        nrow = n,
        ncol = current_batch
      )

      projection_coef <- cache$gram_inverse %*%
        crossprod(cache$D, X_permuted * cache$W)
      X_tilde_permuted <- X_permuted - cache$D %*% projection_coef
      score_permuted <- X_tilde_permuted * cache$outcome_factor
      reference <- .information_statistic(
        score_permuted,
        X_tilde_permuted,
        cache$W
      )

      hits <- hits + sum(abs(reference) >= abs(observed))
      completed <- completed + current_batch
    }

    p_value <- (1 + hits) / (B + 1)
    .method_output(data, p_value = p_value, statistic = observed)
  }

  studentized_permutation_f <- function(data) {
    if (!isTRUE(data$success)) {
      return(.failed_method_output(data))
    }
    if (isTRUE(data$score_cache$degenerate)) {
      return(.degenerate_method_output(data))
    }

    cache <- data$score_cache
    observed <- .studentized_statistic(cache$score_contrib)[1L]
    B <- data$B_resample
    n <- data$n
    batch_size <- 128L
    hits <- 0L
    completed <- 0L

    set.seed(data$perm_seed)
    while (completed < B) {
      current_batch <- min(batch_size, B - completed)
      permutation_indices <- vapply(
        seq_len(current_batch),
        function(i) sample.int(n),
        FUN.VALUE = integer(n)
      )
      X_permuted <- matrix(
        data$X[permutation_indices],
        nrow = n,
        ncol = current_batch
      )

      projection_coef <- cache$gram_inverse %*%
        crossprod(cache$D, X_permuted * cache$W)
      X_tilde_permuted <- X_permuted - cache$D %*% projection_coef
      score_permuted <- X_tilde_permuted * cache$outcome_factor
      reference <- .studentized_statistic(score_permuted)

      hits <- hits + sum(abs(reference) >= abs(observed))
      completed <- completed + current_batch
    }

    p_value <- (1 + hits) / (B + 1)
    .method_output(data, p_value = p_value, statistic = observed)
  }

  .make_rademacher_matrix <- function(B, n, seed) {
    set.seed(seed)
    matrix(
      sample(c(-1L, 1L), size = B * n, replace = TRUE),
      nrow = B,
      ncol = n,
      byrow = TRUE
    )
  }

  hgf_signflip_mc_f <- function(data) {
    if (!isTRUE(data$success)) {
      return(.failed_method_output(data))
    }
    if (isTRUE(data$score_cache$degenerate)) {
      return(.degenerate_method_output(data))
    }

    score <- data$score_cache$score_contrib
    observed <- sum(score)
    signs <- .make_rademacher_matrix(
      B = data$B_resample,
      n = data$n,
      seed = data$flip_seed
    )
    reference <- as.numeric(signs %*% score)
    p_value <- (1 + sum(abs(reference) >= abs(observed))) /
      (data$B_resample + 1)

    .method_output(data, p_value = p_value, statistic = observed)
  }

  .standardized_ratio <- function(numerator, denominator_sq) {
    denominator_sq[denominator_sq < 0 & denominator_sq > -1e-10] <- 0
    if (any(denominator_sq < 0) || any(!is.finite(denominator_sq))) {
      stop("invalid standardized-flip variance correction")
    }

    positive <- denominator_sq > 1e-12
    statistic <- numeric(length(numerator))
    statistic[positive] <- numerator[positive] / sqrt(denominator_sq[positive])
    zero_denom_nonzero_num <- !positive & abs(numerator) > sqrt(.Machine$double.eps)
    statistic[zero_denom_nonzero_num] <-
      sign(numerator[zero_denom_nonzero_num]) * Inf
    statistic
  }

  standardized_signflip_f <- function(data) {
    if (!isTRUE(data$success)) {
      return(.failed_method_output(data))
    }
    if (isTRUE(data$score_cache$degenerate)) {
      return(.degenerate_method_output(data))
    }

    cache <- data$score_cache
    observed <- sum(cache$score_contrib)
    signs <- .make_rademacher_matrix(
      B = data$B_resample,
      n = data$n,
      seed = data$flip_seed
    )

    numerator <- as.numeric(signs %*% cache$score_contrib)
    nuisance_projection <- signs %*% cache$A
    denominator_sq <- 1 - rowSums(nuisance_projection^2)
    reference <- tryCatch(
      .standardized_ratio(numerator, denominator_sq),
      error = function(e) e
    )
    if (inherits(reference, "error")) {
      return(.method_output(
        data = data,
        p_value = NA_real_,
        statistic = observed,
        success = FALSE,
        failure_reason = .failure_message(
          "standardized sign flip failed",
          reference
        )
      ))
    }

    p_value <- (1 + sum(abs(reference) >= abs(observed))) /
      (data$B_resample + 1)
    .method_output(data, p_value = p_value, statistic = observed)
  }
}, envir = .nb_signflip_env)

generate_data <- .nb_signflip_env$generate_data
parametric_score_f <- .nb_signflip_env$parametric_score_f
information_studentized_permutation_f <-
  .nb_signflip_env$information_studentized_permutation_f
studentized_permutation_f <- .nb_signflip_env$studentized_permutation_f
hgf_signflip_mc_f <- .nb_signflip_env$hgf_signflip_mc_f
standardized_signflip_f <- .nb_signflip_env$standardized_signflip_f

generate_data_function <- simulatr::simulatr_function(
  f = generate_data,
  arg_names = formalArgs(generate_data),
  packages = "MASS",
  loop = TRUE
)

studentized_permutation_spec_f <- simulatr::simulatr_function(
  f = studentized_permutation_f,
  arg_names = character(0),
  loop = TRUE
)
parametric_score_spec_f <- simulatr::simulatr_function(
  f = parametric_score_f,
  arg_names = character(0),
  loop = TRUE
)
information_studentized_permutation_spec_f <- simulatr::simulatr_function(
  f = information_studentized_permutation_f,
  arg_names = character(0),
  loop = TRUE
)
hgf_signflip_mc_spec_f <- simulatr::simulatr_function(
  f = hgf_signflip_mc_f,
  arg_names = character(0),
  loop = TRUE
)
standardized_signflip_spec_f <- simulatr::simulatr_function(
  f = standardized_signflip_f,
  arg_names = character(0),
  loop = TRUE
)

run_method_functions <- list(
  parametric_score = parametric_score_spec_f,
  information_studentized_permutation =
    information_studentized_permutation_spec_f,
  studentized_permutation = studentized_permutation_spec_f,
  hgf_signflip_mc = hgf_signflip_mc_spec_f,
  standardized_signflip = standardized_signflip_spec_f
)

############################
# 4. Assemble specifier
############################
simulatr_spec <- simulatr::simulatr_specifier(
  parameter_grid = parameter_grid,
  fixed_parameters = fixed_parameters,
  generate_data_function = generate_data_function,
  run_method_functions = run_method_functions
)

# Development checks:
# simulatr::check_simulatr_specifier_object(
#   simulatr_spec, B_in = 2, parallel = FALSE
# )
