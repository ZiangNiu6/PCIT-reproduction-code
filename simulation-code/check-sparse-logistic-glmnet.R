#!/usr/bin/env Rscript

###############################################################################
# Sparse high-dimensional X | Z check for the GCM and PCM.
#
# Data-generating process (the conditional-independence null):
#
#   Z_i ~ N_p(0, I_p)
#   X_i | Z_i ~ Bernoulli(expit(-3 + Z_i^T beta_x))
#   Y_i | Z_i ~ N(Z_i^T beta_y, 1), independently of X_i | Z_i,
#
# with n = 200, p = 50, s = 5, and every active coefficient equal to
# 1 / sqrt(s).  By default the two supports overlap in all s coordinates, but
# --overlap=k moves s-k active beta_y coordinates onto otherwise inactive
# covariates.  Hence ||beta_x|| = ||beta_y|| = 1 while the amount of marginal
# confounding is controlled by k, and X is nevertheless conditionally
# independent of Y given Z.
#
# All observed nuisance fits are five-fold cross-validated lasso fits from
# glmnet.  The CV folds and the PCM 50/50 sample split are fixed by
# purpose-specific seeds.
# For GCM, every raw-X permutation refits X | Z.  For PCM, the D2 direction is
# held fixed while its value is recomputed at the permuted D1 exposure and the
# D1 direction-on-Z projection is refit.  Within each GCM or PCM pair, the
# permutation and asymptotic p-values use exactly the same observed statistic.
#
# By default, each permutation refit repeats cross-validation.  With
# --locked-lambda, cross-validation is run only on the observed response and
# every permutation refit uses plain glmnet at that observed lambda.min.  Thus
# the two modes have identical observed fits and differ only in whether tuning
# is rerun or locked along the permutation orbit.  Because a tuning parameter
# selected at the unpermuted response is not permutation-invariant, the locked
# mode is a sensitivity check rather than a finite-sample exact permutation
# procedure.
#
# IMPORTANT: pcmt delegates reg_XonZ = "lasso" to comets::predict.lasso.  For
# a binomial glmnet fit that method uses glmnet's default link-scale prediction,
# which would form X - logit(phat), not X - phat.  The GCM is implemented here
# directly so its logistic-lasso predictions explicitly use type = "response".
# PCM is also implemented directly so its projection refit can use either
# cross-validated or locked tuning without changing the observed statistic.
#
# Calibration:
#   GCM: two-sided, adaptive permutation versus N(0, 1)
#   PCM: greater,   adaptive permutation versus N(0, 1)
#
# Adaptive permutation uses the repository's h = 15 Besag-Clifford rule with a
# cap of B = 499 by default.  Pass --fixed-permutation to use the ordinary
# fixed-B plus-one Monte Carlo p-value instead.
#
# Usage:
#   Rscript --vanilla check-sparse-logistic-glmnet.R
#   Rscript --vanilla check-sparse-logistic-glmnet.R --quick  # smoke only
#   Rscript --vanilla check-sparse-logistic-glmnet.R \
#     --reps=500 --B=499 --cores=4
#   Rscript --vanilla check-sparse-logistic-glmnet.R \
#     --locked-lambda --reps=200 --B=499 --cores=4
#   Rscript --vanilla check-sparse-logistic-glmnet.R \
#     --overlap=2 --locked-lambda --reps=200 --B=499 --cores=4
###############################################################################

suppressPackageStartupMessages({
  library(glmnet)
  library(parallel)
})

## ---------------------------------------------------------------- arguments
cli <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% cli
fixed_permutation <- "--fixed-permutation" %in% cli
locked_lambda <- "--locked-lambda" %in% cli

get_numeric_option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), cli, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  if (length(hit) > 1L) {
    stop("Option --", name, " was supplied more than once.", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(sub(paste0("^--", name, "="), "", hit)))
  if (length(value) != 1L || !is.finite(value)) {
    stop("Option --", name, " must be numeric.", call. = FALSE)
  }
  value
}

get_character_option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), cli, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  if (length(hit) > 1L) {
    stop("Option --", name, " was supplied more than once.", call. = FALSE)
  }
  sub(paste0("^--", name, "="), "", hit)
}

recognized_flags <- c("--quick", "--fixed-permutation", "--locked-lambda")
recognized_value_prefixes <- c(
  "--reps=", "--B=", "--cores=", "--overlap=", "--out="
)
unknown <- cli[
  !(cli %in% recognized_flags) &
    !vapply(cli, function(x) {
      any(startsWith(x, recognized_value_prefixes))
    }, logical(1L))
]
if (length(unknown) > 0L) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

n_reps_input <- get_numeric_option("reps", if (quick) 4L else 200L)
B_max_input <- get_numeric_option("B", if (quick) 7L else 499L)
n_cores_input <- get_numeric_option("cores", if (quick) 1L else 4L)
overlap_input <- get_numeric_option("overlap", 5L)
output_override <- get_character_option("out", "")

integer_options <- c(
  reps = n_reps_input,
  B = B_max_input,
  cores = n_cores_input
)
if (any(integer_options < 1) ||
    any(integer_options != floor(integer_options)) ||
    any(integer_options > .Machine$integer.max)) {
  stop("--reps, --B, and --cores must be positive integers.", call. = FALSE)
}
if (overlap_input < 0 || overlap_input != floor(overlap_input) ||
    overlap_input > .Machine$integer.max) {
  stop("--overlap must be a nonnegative integer.", call. = FALSE)
}
n_reps <- as.integer(n_reps_input)
B_max <- as.integer(B_max_input)
n_cores <- as.integer(n_cores_input)
SUPPORT_OVERLAP <- as.integer(overlap_input)
if (n_cores > n_reps) {
  n_cores <- n_reps
}

## --------------------------------------------------------------- fixed design
N <- 200L
P <- 50L
S <- 5L
X_INTERCEPT <- -3
Y_NOISE_SD <- 1
ALPHA <- 0.05
ADAPTIVE_H <- 15L
CV_FOLDS <- 5L
GLMNET_ALPHA <- 1
GLMNET_NLAMBDA <- 100L
LAMBDA_RULE <- "lambda.min"
SCHEMA_VERSION <- "sparse_logistic_glmnet_gcm_pcm_v3"

if (SUPPORT_OVERLAP > S) {
  stop("--overlap cannot exceed s = ", S, ".", call. = FALSE)
}
if (2L * S - SUPPORT_OVERLAP > P || N %% 2L != 0L ||
    N / 2L < CV_FOLDS) {
  stop("The fixed design is internally inconsistent.", call. = FALSE)
}

beta_x <- numeric(P)
beta_y <- numeric(P)
x_support <- seq_len(S)
y_support <- seq.int(S - SUPPORT_OVERLAP + 1L, 2L * S - SUPPORT_OVERLAP)
beta_x[x_support] <- 1 / sqrt(S)
beta_y[y_support] <- 1 / sqrt(S)
linear_predictor_correlation <- sum(beta_x * beta_y) /
  sqrt(sum(beta_x^2) * sum(beta_y^2))
if (length(intersect(x_support, y_support)) != SUPPORT_OVERLAP ||
    !isTRUE(all.equal(sum(beta_x^2), 1, tolerance = 1e-14)) ||
    !isTRUE(all.equal(sum(beta_y^2), 1, tolerance = 1e-14)) ||
    !isTRUE(all.equal(
      linear_predictor_correlation,
      SUPPORT_OVERLAP / S,
      tolerance = 1e-14
    ))) {
  stop("Failed to construct the requested sparse supports.", call. = FALSE)
}

expected_prevalence <- stats::integrate(
  function(z) stats::plogis(X_INTERCEPT + z) * stats::dnorm(z),
  lower = -Inf,
  upper = Inf,
  rel.tol = 1e-10
)$value
expected_logistic_derivative <- stats::integrate(
  function(z) {
    probability <- stats::plogis(X_INTERCEPT + z)
    probability * (1 - probability) * stats::dnorm(z)
  },
  lower = -Inf,
  upper = Inf,
  rel.tol = 1e-10
)$value
expected_marginal_xy_covariance <-
  linear_predictor_correlation * expected_logistic_derivative
expected_y_variance <- sum(beta_y^2) + Y_NOISE_SD^2
expected_marginal_xy_correlation <- expected_marginal_xy_covariance / sqrt(
  expected_prevalence * (1 - expected_prevalence) * expected_y_variance
)

## --------------------------------------------------------- deterministic RNG
with_local_seed <- function(seed, code) {
  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(R.utils::withSeed(code, seed = as.integer(seed)))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

balanced_foldid <- function(n, k, seed) {
  foldid <- with_local_seed(
    seed,
    sample(rep(seq_len(k), length.out = n), size = n, replace = FALSE)
  )
  if (!identical(sort(unique(foldid)), seq_len(k)) ||
      max(tabulate(foldid, nbins = k)) - min(tabulate(foldid, nbins = k)) > 1L) {
    stop("Failed to construct balanced CV folds.", call. = FALSE)
  }
  as.integer(foldid)
}

binomial_foldid <- function(y, k, seed) {
  y <- as.numeric(y)
  counts <- table(y)
  if (length(counts) < 2L || min(counts) < 3L) {
    return(balanced_foldid(length(y), k, seed))
  }

  foldid <- integer(length(y))
  foldid <- with_local_seed(seed, {
    for (value in sort(unique(y))) {
      indices <- which(y == value)
      indices <- sample(indices, length(indices), replace = FALSE)
      foldid[indices] <- rep(seq_len(k), length.out = length(indices))
    }
    foldid
  })

  training_has_both_classes <- vapply(seq_len(k), function(fold) {
    training_counts <- table(y[foldid != fold])
    length(training_counts) == 2L && min(training_counts) >= 2L
  }, logical(1L))
  if (!all(training_has_both_classes)) {
    stop("Failed to construct viable stratified binomial CV folds.")
  }
  as.integer(foldid)
}

pcm_D2_indices <- with_local_seed(831001L, sample.int(N, N / 2L))
pcm_D2_indices <- as.integer(pcm_D2_indices)
pcm_D1_indices <- setdiff(seq_len(N), pcm_D2_indices)

foldids <- list(
  gcm_YonZ = balanced_foldid(N, CV_FOLDS, 831011L),
  pcm_YonXZ = balanced_foldid(length(pcm_D2_indices), CV_FOLDS, 831013L),
  pcm_YhatonZ = balanced_foldid(length(pcm_D2_indices), CV_FOLDS, 831014L),
  pcm_YonZ = balanced_foldid(length(pcm_D1_indices), CV_FOLDS, 831015L),
  pcm_RonZ = balanced_foldid(length(pcm_D1_indices), CV_FOLDS, 831016L)
)

if (length(intersect(pcm_D1_indices, pcm_D2_indices)) != 0L ||
    !setequal(c(pcm_D1_indices, pcm_D2_indices), seq_len(N))) {
  stop("The fixed PCM split is not a partition.", call. = FALSE)
}

replicate_seed <- function(rep_id, stream) {
  seeds <- c(data = 832000L, gcm_permutation = 833000L,
             pcm_permutation = 834000L)
  as.integer(seeds[[stream]] + 1009L * as.integer(rep_id))
}

## --------------------------------------------------------------- glmnet fits
constant_mean_fit <- function(value, family) {
  fit <- structure(
    list(value = as.numeric(value), family = family),
    class = "constant_mean_fit"
  )
  attr(fit, "selected_lambda") <- Inf
  fit
}

fit_glmnet_mean <- function(y, x, family, foldid) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  if (nrow(x) != length(y) || length(foldid) != length(y)) {
    stop("glmnet response, design, and foldid lengths disagree.", call. = FALSE)
  }

  is_constant <- if (identical(family, "binomial")) {
    counts <- table(y)
    length(counts) < 2L || min(counts) < 3L
  } else {
    !is.finite(stats::var(y)) || stats::var(y) <= .Machine$double.eps
  }
  if (is_constant) {
    return(constant_mean_fit(mean(y), family))
  }

  type_measure <- if (identical(family, "binomial")) "deviance" else "mse"
  fit <- suppressWarnings(glmnet::cv.glmnet(
    x = x,
    y = y,
    family = family,
    alpha = GLMNET_ALPHA,
    foldid = foldid,
    type.measure = type_measure,
    standardize = TRUE,
    intercept = TRUE,
    nlambda = GLMNET_NLAMBDA,
    parallel = FALSE
  ))
  attr(fit, "selected_lambda") <- as.numeric(fit[[LAMBDA_RULE]])
  fit
}

fit_glmnet_at_lambda <- function(y, x, family, lambda) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  is_constant <- if (identical(family, "binomial")) {
    length(unique(y)) < 2L
  } else {
    !is.finite(stats::var(y)) || stats::var(y) <= .Machine$double.eps
  }
  if (is_constant || length(lambda) != 1L || !is.finite(lambda)) {
    return(constant_mean_fit(mean(y), family))
  }

  fit <- suppressWarnings(glmnet::glmnet(
    x = x,
    y = y,
    family = family,
    alpha = GLMNET_ALPHA,
    lambda = as.numeric(lambda),
    standardize = TRUE,
    intercept = TRUE
  ))
  attr(fit, "selected_lambda") <- as.numeric(lambda)
  fit
}

selected_lambda <- function(fit) {
  as.numeric(attr(fit, "selected_lambda"))
}

predict_glmnet_mean <- function(fit, newx) {
  if (inherits(fit, "constant_mean_fit")) {
    return(rep(fit$value, nrow(as.matrix(newx))))
  }
  as.numeric(stats::predict(
    fit,
    newx = as.matrix(newx),
    s = selected_lambda(fit),
    type = "response"
  ))
}

## ---------------------------------------------------------- test statistics
studentized_statistic <- function(score) {
  score <- as.numeric(score)
  score_mean <- mean(score)
  score_variance <- mean(score^2) - score_mean^2
  if (!is.finite(score_variance) || score_variance <= 0) {
    return(0)
  }
  sqrt(length(score)) * score_mean / sqrt(score_variance)
}

score_map <- function(statistic, alternative) {
  switch(
    alternative,
    greater = statistic,
    less = -statistic,
    two.sided = abs(statistic),
    stop("Unknown alternative: ", alternative, call. = FALSE)
  )
}

adaptive_pvalue <- function(t, K, h, t_star) {
  if (K < h) h / (t + h - K) else h / t_star
}

fixed_monte_carlo_cutoff <- function(
    permutation_statistics,
    alternative,
    alpha) {
  B <- length(permutation_statistics)
  maximum_exceedances <- floor(alpha * (B + 1) - 1)
  if (maximum_exceedances < 0L) {
    return(Inf)
  }
  permutation_scores <- sort(vapply(
    permutation_statistics,
    score_map,
    numeric(1L),
    alternative = alternative
  ))
  permutation_scores[[B - maximum_exceedances]]
}

permutation_reference <- function(
    observed,
    draw_statistic,
    B,
    alternative,
    adaptive,
    h,
    alpha) {
  observed_score <- score_map(observed, alternative)
  permutation_statistics <- numeric(B)

  if (!adaptive) {
    for (b in seq_len(B)) {
      permutation_statistics[[b]] <- draw_statistic()
    }
    permutation_scores <- vapply(
      permutation_statistics,
      score_map,
      numeric(1L),
      alternative = alternative
    )
    return(list(
      p_value = (1 + sum(permutation_scores >= observed_score)) / (B + 1),
      statistics = permutation_statistics,
      actual_B = as.integer(B),
      K = as.integer(sum(permutation_scores >= observed_score)),
      t_star = NA_real_,
      stopped_early = FALSE,
      stop_reason = "fixed_cap"
    ))
  }

  K <- 0L
  t_star <- Inf
  stop_reason <- "cap"
  actual_B <- B
  p_value <- 1
  for (b in seq_len(B)) {
    permutation_statistics[[b]] <- draw_statistic()
    if (score_map(permutation_statistics[[b]], alternative) >= observed_score) {
      K <- K + 1L
      if (K == h && is.infinite(t_star)) {
        t_star <- b
      }
    }
    p_value <- adaptive_pvalue(b, K, h, t_star)
    if (p_value <= alpha) {
      actual_B <- b
      stop_reason <- "rejection"
      break
    }
    if (K == h && p_value > alpha) {
      actual_B <- b
      stop_reason <- "futility"
      break
    }
  }

  list(
    p_value = as.numeric(p_value),
    statistics = permutation_statistics[seq_len(actual_B)],
    actual_B = as.integer(actual_B),
    K = as.integer(K),
    t_star = as.numeric(t_star),
    stopped_early = actual_B < B,
    stop_reason = stop_reason
  )
}

## ----------------------------------------------------------- data generation
generate_null_data <- function(rep_id) {
  with_local_seed(replicate_seed(rep_id, "data"), {
    Z <- matrix(stats::rnorm(N * P), nrow = N, ncol = P)
    colnames(Z) <- paste0("Z", seq_len(P))
    x_probability <- stats::plogis(X_INTERCEPT + c(Z %*% beta_x))
    X <- stats::rbinom(N, size = 1L, prob = x_probability)
    Y_mean <- c(Z %*% beta_y)
    Y <- Y_mean + stats::rnorm(N, sd = Y_NOISE_SD)
    list(
      X = X,
      Y = Y,
      Z = Z,
      x_probability = x_probability,
      Y_mean = Y_mean
    )
  })
}

## -------------------------------------------------------------- one dataset
run_one_replicate <- function(rep_id) {
  data <- generate_null_data(rep_id)

  # GCM: direct glmnet implementation is required for response-scale binomial
  # predictions.  The Y-side fit is held fixed; X | Z is refit after every raw
  # X permutation, exactly as in pcmt::perm_gcm.
  gcm_y_fit <- fit_glmnet_mean(
    data$Y, data$Z, family = "gaussian", foldid = foldids$gcm_YonZ
  )
  gcm_y_prediction <- predict_glmnet_mean(gcm_y_fit, data$Z)
  gcm_rY <- data$Y - gcm_y_prediction

  gcm_x_fit <- fit_glmnet_mean(
    data$X,
    data$Z,
    family = "binomial",
    foldid = binomial_foldid(data$X, CV_FOLDS, 831012L)
  )
  gcm_observed_lambda <- selected_lambda(gcm_x_fit)
  gcm_x_probability <- predict_glmnet_mean(gcm_x_fit, data$Z)
  if (any(!is.finite(gcm_x_probability)) ||
      any(gcm_x_probability < 0 | gcm_x_probability > 1)) {
    stop("The logistic glmnet fit returned an invalid response probability.")
  }
  gcm_rX <- data$X - gcm_x_probability
  gcm_statistic <- studentized_statistic(gcm_rX * gcm_rY)

  gcm_permutation <- with_local_seed(
    replicate_seed(rep_id, "gcm_permutation"),
    permutation_reference(
      observed = gcm_statistic,
      draw_statistic = function() {
        permuted_X <- data$X[sample.int(N)]
        permuted_fit <- if (locked_lambda) {
          fit_glmnet_at_lambda(
            permuted_X,
            data$Z,
            family = "binomial",
            lambda = gcm_observed_lambda
          )
        } else {
          fit_glmnet_mean(
            permuted_X,
            data$Z,
            family = "binomial",
            foldid = binomial_foldid(permuted_X, CV_FOLDS, 831012L)
          )
        }
        permuted_probability <- predict_glmnet_mean(permuted_fit, data$Z)
        if (any(!is.finite(permuted_probability)) ||
            any(permuted_probability < 0 | permuted_probability > 1)) {
          stop("A permuted logistic glmnet fit returned an invalid probability.")
        }
        studentized_statistic((permuted_X - permuted_probability) * gcm_rY)
      },
      B = B_max,
      alternative = "two.sided",
      adaptive = !fixed_permutation,
      h = ADAPTIVE_H,
      alpha = ALPHA
    )
  )

  # PCM: learn the direction on D2, then evaluate and residualize it on D1.
  # Only the final D1 direction-on-Z fit is response-adaptive under a D1 raw-X
  # permutation, so it is the only PCM nuisance refit along the orbit.
  Y_D1 <- data$Y[pcm_D1_indices]
  Y_D2 <- data$Y[pcm_D2_indices]
  X_D1 <- data$X[pcm_D1_indices]
  X_D2 <- data$X[pcm_D2_indices]
  Z_D1 <- data$Z[pcm_D1_indices, , drop = FALSE]
  Z_D2 <- data$Z[pcm_D2_indices, , drop = FALSE]
  XZ_D2 <- cbind(X = X_D2, Z_D2)

  pcm_full_fit <- fit_glmnet_mean(
    Y_D2,
    XZ_D2,
    family = "gaussian",
    foldid = foldids$pcm_YonXZ
  )
  pcm_full_prediction_D2 <- predict_glmnet_mean(pcm_full_fit, XZ_D2)
  pcm_projected_fit <- fit_glmnet_mean(
    pcm_full_prediction_D2,
    Z_D2,
    family = "gaussian",
    foldid = foldids$pcm_YhatonZ
  )
  pcm_direction <- function(exposure, covariates) {
    predict_glmnet_mean(
      pcm_full_fit,
      cbind(X = as.numeric(exposure), covariates)
    ) - predict_glmnet_mean(pcm_projected_fit, covariates)
  }

  pcm_y_fit <- fit_glmnet_mean(
    Y_D1,
    Z_D1,
    family = "gaussian",
    foldid = foldids$pcm_YonZ
  )
  pcm_rY <- Y_D1 - predict_glmnet_mean(pcm_y_fit, Z_D1)
  pcm_direction_D1 <- pcm_direction(X_D1, Z_D1)
  pcm_projection_fit <- fit_glmnet_mean(
    pcm_direction_D1,
    Z_D1,
    family = "gaussian",
    foldid = foldids$pcm_RonZ
  )
  pcm_observed_lambda <- selected_lambda(pcm_projection_fit)
  pcm_rF <- pcm_direction_D1 - predict_glmnet_mean(
    pcm_projection_fit, Z_D1
  )
  pcm_statistic <- studentized_statistic(pcm_rY * pcm_rF)

  pcm_permutation <- with_local_seed(
    replicate_seed(rep_id, "pcm_permutation"),
    permutation_reference(
      observed = pcm_statistic,
      draw_statistic = function() {
        permuted_X <- X_D1[sample.int(length(X_D1))]
        permuted_direction <- pcm_direction(permuted_X, Z_D1)
        permuted_projection_fit <- if (locked_lambda) {
          fit_glmnet_at_lambda(
            permuted_direction,
            Z_D1,
            family = "gaussian",
            lambda = pcm_observed_lambda
          )
        } else {
          fit_glmnet_mean(
            permuted_direction,
            Z_D1,
            family = "gaussian",
            foldid = foldids$pcm_RonZ
          )
        }
        permuted_rF <- permuted_direction - predict_glmnet_mean(
          permuted_projection_fit, Z_D1
        )
        studentized_statistic(pcm_rY * permuted_rF)
      },
      B = B_max,
      alternative = "greater",
      adaptive = !fixed_permutation,
      h = ADAPTIVE_H,
      alpha = ALPHA
    )
  )

  p_values <- c(
    perm_gcm_glmnet = as.numeric(gcm_permutation$p_value),
    asym_gcm_glmnet = 2 * stats::pnorm(
      abs(gcm_statistic), lower.tail = FALSE
    ),
    perm_pcm_glmnet = as.numeric(pcm_permutation$p_value),
    asym_pcm_glmnet = stats::pnorm(
      pcm_statistic, lower.tail = FALSE
    )
  )
  if (any(!is.finite(p_values)) || any(p_values < 0 | p_values > 1)) {
    stop("A method returned an invalid p-value.")
  }

  gcm_cutoff <- if (fixed_permutation) {
    fixed_monte_carlo_cutoff(
      gcm_permutation$statistics,
      alternative = "two.sided",
      alpha = ALPHA
    )
  } else {
    NA_real_
  }
  pcm_cutoff <- if (fixed_permutation) {
    fixed_monte_carlo_cutoff(
      pcm_permutation$statistics,
      alternative = "greater",
      alpha = ALPHA
    )
  } else {
    NA_real_
  }

  data.frame(
    rep = as.integer(rep_id),
    support_overlap = SUPPORT_OVERLAP,
    linear_predictor_correlation = linear_predictor_correlation,
    expected_marginal_xy_correlation = expected_marginal_xy_correlation,
    carrier_count = as.integer(sum(data$X)),
    observed_prevalence = mean(data$X),
    mean_true_probability = mean(data$x_probability),
    gcm_statistic = gcm_statistic,
    pcm_statistic = pcm_statistic,
    gcm_observed_lambda = gcm_observed_lambda,
    pcm_observed_projection_lambda = pcm_observed_lambda,
    permutation_refit = if (locked_lambda) "locked_lambda" else "repeat_cv",
    perm_gcm_glmnet = p_values[["perm_gcm_glmnet"]],
    asym_gcm_glmnet = p_values[["asym_gcm_glmnet"]],
    perm_pcm_glmnet = p_values[["perm_pcm_glmnet"]],
    asym_pcm_glmnet = p_values[["asym_pcm_glmnet"]],
    gcm_actual_B = as.integer(gcm_permutation$actual_B),
    pcm_actual_B = as.integer(pcm_permutation$actual_B),
    gcm_stop_reason = as.character(gcm_permutation$stop_reason),
    pcm_stop_reason = as.character(pcm_permutation$stop_reason),
    gcm_constant_x_fit = inherits(gcm_x_fit, "constant_mean_fit"),
    pcm_degenerate_direction =
      !is.finite(stats::var(pcm_rF)) ||
      stats::var(pcm_rF) <= .Machine$double.eps,
    gcm_permutation_cutoff = as.numeric(gcm_cutoff),
    pcm_permutation_cutoff = as.numeric(pcm_cutoff),
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------------- execute
cat(sprintf(
  paste0(
    "Sparse logistic glmnet CI-null check: n=%d, p=%d, s=%d, ",
    "intercept=%.1f, R=%d, B<=%d, %d core(s)\n"
  ),
  N, P, S, X_INTERCEPT, n_reps, B_max, n_cores
))
cat(sprintf(
  "Expected marginal exposure prevalence: %.4f (about %.1f carriers)\n",
  expected_prevalence,
  N * expected_prevalence
))
cat(sprintf(
  paste0(
    "Active supports: beta_x={%s}; beta_y={%s}; overlap=%d/%d; ",
    "linear-predictor correlation=%.3f; implied marginal Corr(X,Y)=%.4f\n"
  ),
  paste(x_support, collapse = ","),
  paste(y_support, collapse = ","),
  SUPPORT_OVERLAP,
  S,
  linear_predictor_correlation,
  expected_marginal_xy_correlation
))
cat(sprintf(
  "Calibration: %s permutation; glmnet alpha=%.0f, %d-fold CV, %s\n\n",
  if (fixed_permutation) "fixed-B" else paste0("adaptive h=", ADAPTIVE_H),
  GLMNET_ALPHA,
  CV_FOLDS,
  LAMBDA_RULE
))
cat(
  "Permutation nuisance refits: ",
  if (locked_lambda) {
    "plain glmnet at the corresponding observed lambda.min"
  } else {
    "repeat cv.glmnet"
  },
  "\n\n",
  sep = ""
)
utils::flush.console()

start_time <- proc.time()[["elapsed"]]
replicate_results <- parallel::mclapply(
  seq_len(n_reps),
  function(rep_id) {
    tryCatch(
      list(ok = TRUE, value = run_one_replicate(rep_id)),
      error = function(error) {
        list(ok = FALSE, rep = rep_id, error = conditionMessage(error))
      }
    )
  },
  mc.cores = n_cores,
  mc.set.seed = FALSE
)
elapsed_seconds <- proc.time()[["elapsed"]] - start_time

failed <- !vapply(replicate_results, function(x) isTRUE(x$ok), logical(1L))
if (any(failed)) {
  messages <- vapply(
    replicate_results[failed],
    function(x) paste0("rep ", x$rep, ": ", x$error),
    character(1L)
  )
  stop(
    sum(failed), " replicate(s) failed:\n", paste(messages, collapse = "\n"),
    call. = FALSE
  )
}
replicates <- do.call(rbind, lapply(replicate_results, `[[`, "value"))
rownames(replicates) <- NULL

## ---------------------------------------------------------------- summarize
wilson_interval <- function(rejections, total, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  proportion <- rejections / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half_width <- z / denominator * sqrt(
    proportion * (1 - proportion) / total + z^2 / (4 * total^2)
  )
  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

method_contract <- data.frame(
  method = c(
    "perm_gcm_glmnet", "asym_gcm_glmnet",
    "perm_pcm_glmnet", "asym_pcm_glmnet"
  ),
  test = c("GCM", "GCM", "PCM", "PCM"),
  calibration = c("Permutation", "Asymptotic", "Permutation", "Asymptotic"),
  alternative = c("two.sided", "two.sided", "greater", "greater"),
  stringsAsFactors = FALSE
)

summary_rows <- lapply(seq_len(nrow(method_contract)), function(index) {
  contract <- method_contract[index, , drop = FALSE]
  p_values <- replicates[[contract$method]]
  rejection_count <- sum(p_values <= ALPHA)
  interval <- wilson_interval(rejection_count, n_reps)

  reference_cutoff <- if (contract$calibration == "Asymptotic") {
    if (contract$test == "GCM") stats::qnorm(1 - ALPHA / 2) else
      stats::qnorm(1 - ALPHA)
  } else if (fixed_permutation) {
    mean(replicates[[
      if (contract$test == "GCM") "gcm_permutation_cutoff" else
        "pcm_permutation_cutoff"
    ]])
  } else {
    NA_real_
  }

  data.frame(
    support_overlap = SUPPORT_OVERLAP,
    linear_predictor_correlation = linear_predictor_correlation,
    expected_marginal_xy_correlation = expected_marginal_xy_correlation,
    permutation_refit = if (locked_lambda) "locked_lambda" else "repeat_cv",
    method = contract$method,
    test = contract$test,
    calibration = contract$calibration,
    alternative = contract$alternative,
    replications = n_reps,
    rejections = rejection_count,
    rejection_rate = rejection_count / n_reps,
    monte_carlo_se = sqrt(
      (rejection_count / n_reps) * (1 - rejection_count / n_reps) / n_reps
    ),
    lower_95_wilson = interval[["lower"]],
    upper_95_wilson = interval[["upper"]],
    mean_p_value = mean(p_values),
    reference_cutoff = reference_cutoff,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL

configuration <- list(
  schema_version = SCHEMA_VERSION,
  null = "X conditionally independent of Y given Z",
  n = N,
  p = P,
  s = S,
  x_intercept = X_INTERCEPT,
  beta_x = beta_x,
  beta_y = beta_y,
  beta_x_support = x_support,
  beta_y_support = y_support,
  support_union_size = length(union(x_support, y_support)),
  support_overlap = SUPPORT_OVERLAP,
  linear_predictor_correlation = linear_predictor_correlation,
  expected_logistic_derivative = expected_logistic_derivative,
  expected_marginal_xy_covariance = expected_marginal_xy_covariance,
  expected_marginal_xy_correlation = expected_marginal_xy_correlation,
  coefficient_value = 1 / sqrt(S),
  y_noise_sd = Y_NOISE_SD,
  expected_marginal_prevalence = expected_prevalence,
  alpha = ALPHA,
  replications = n_reps,
  permutation_cap = B_max,
  permutation_calibration = if (fixed_permutation) "fixed_B_plus_one" else
    "adaptive_Besag_Clifford",
  adaptive_h = if (fixed_permutation) NA_integer_ else ADAPTIVE_H,
  permutation_refit = if (locked_lambda) {
    "glmnet_at_corresponding_observed_lambda_min"
  } else {
    "repeat_cv_glmnet"
  },
  permutation_interpretation = if (locked_lambda) {
    paste(
      "sensitivity check: lambda selected at the unpermuted response is",
      "held fixed and is not permutation-invariant"
    )
  } else {
    "the full response-adaptive fitting rule is repeated along the orbit"
  },
  glmnet_alpha = GLMNET_ALPHA,
  glmnet_nlambda = GLMNET_NLAMBDA,
  lambda_rule = LAMBDA_RULE,
  cv_folds = CV_FOLDS,
  binomial_cv_fold_policy = paste(
    "response-stratified folds; intercept-only fallback when either class",
    "has fewer than three observations"
  ),
  pcm_D1_indices = pcm_D1_indices,
  pcm_D2_indices = pcm_D2_indices,
  foldids = foldids,
  cores = n_cores,
  elapsed_seconds = elapsed_seconds,
  package_versions = c(
    R = as.character(getRversion()),
    glmnet = as.character(utils::packageVersion("glmnet"))
  )
)

output_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
output_name <- if (locked_lambda) {
  "sparse-logistic-glmnet-check-locked-lambda"
} else {
  "sparse-logistic-glmnet-check-repeat-cv"
}
if (SUPPORT_OVERLAP != S) {
  output_name <- paste0(output_name, "-overlap", SUPPORT_OVERLAP)
}
output_dir <- if (nzchar(output_override)) {
  output_override
} else if (nzchar(output_root)) {
  file.path(output_root, "simulation-results", output_name)
} else {
  file.path(getwd(), "results", output_name)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

result <- list(
  schema_version = SCHEMA_VERSION,
  configuration = configuration,
  summary = summary,
  replicates = replicates
)
saveRDS(
  result,
  file.path(output_dir, "sparse-logistic-glmnet-check-results.rds")
)
utils::write.csv(
  summary,
  file.path(output_dir, "sparse-logistic-glmnet-check-summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  replicates,
  file.path(output_dir, "sparse-logistic-glmnet-check-replicates.csv"),
  row.names = FALSE
)

cat(sprintf("Completed in %.1f seconds.\n", elapsed_seconds))
cat(sprintf(
  "Observed mean exposure prevalence: %.4f (mean %.1f carriers)\n",
  mean(replicates$observed_prevalence),
  mean(replicates$carrier_count)
))
cat(sprintf(
  "Degenerate fits: GCM X | Z %d/%d; PCM direction %d/%d\n\n",
  sum(replicates$gcm_constant_x_fit), n_reps,
  sum(replicates$pcm_degenerate_direction), n_reps
))
print(
  summary[, c(
    "test", "calibration", "alternative", "rejections", "replications",
    "rejection_rate", "monte_carlo_se", "lower_95_wilson",
    "upper_95_wilson", "reference_cutoff"
  )],
  row.names = FALSE,
  digits = 4
)
cat("\nWritten to: ", output_dir, "\n", sep = "")
