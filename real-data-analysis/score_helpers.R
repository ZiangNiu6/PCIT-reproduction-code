## score_helpers.R
##
## Efficient-score primitives for the five Frangieh calibrations:
##   1. parametric score;
##   2. information-studentized label permutation;
##   3. empirically studentized label permutation;
##   4. HGF Monte Carlo sign flipping (without SPA); and
##   5. standardized sign flipping.
##
## Every method uses the same fitted negative-binomial null model and the same
## efficient-score contributions.  Only the reference distribution and, for
## the two permutation methods, the studentization differ.

suppressPackageStartupMessages({
  library(MASS)
})

## ---- Robust NB dispersion fallback -----------------------------------------
## Roughly one in seven of the intentionally low-count null genes
## (mean count <~ 0.5, Var/mean ~= 1) makes MASS::glm.nb's theta.ml Newton
## step diverge, which crashes the whole replicate.  The pilot -> MLE ->
## method-of-moments ladder below reproduces spacrt's C++ estimate_theta()
## (spacrt/src/nb_disp_param_precomputation.cpp, itself adapted from sceptre,
## GPL-3) in pure R, so no compiled dependency is added to this pipeline.  The
## loops use isTRUE() around the tolerance test so a NaN step exits the loop
## and flags a warning -- matching C++'s NaN-comparison-is-false semantics --
## instead of raising "missing value where TRUE/FALSE needed".

## Pearson pilot theta = n / sum((y/mu - 1)^2); the glm.nb Newton starting value.
nb_theta_pilot <- function(y, mu) {
  length(y) / sum((y / mu - 1)^2)
}

## Profile-likelihood MLE Newton iteration (score/information in theta).
nb_theta_mle <- function(t0, y, mu, limit = 200L,
                         eps = sqrt(.Machine$double.eps)) {
  it <- 0L
  del <- 1
  while ((it <- it + 1L) < limit && isTRUE(abs(del) > eps)) {
    t0 <- abs(t0)
    score <- sum(digamma(t0 + y) - digamma(t0) + log(t0) + 1 -
                   log(t0 + mu) - (y + t0) / (mu + t0))
    info <- sum(trigamma(t0) - trigamma(t0 + y) - 1 / t0 +
                  2 / (mu + t0) - (y + t0) / (mu + t0)^2)
    del <- score / info
    t0 <- t0 + del
  }
  c(estimate = t0,
    warning = as.numeric(t0 < 0 || it == limit || !is.finite(t0)))
}

## Method-of-moments fallback matching the Pearson chi-square to its residual df.
nb_theta_mm <- function(t0, y, mu, dfr, limit = 200L,
                        eps = sqrt(.Machine$double.eps)) {
  it <- 0L
  del <- 1
  while ((it <- it + 1L) < limit && isTRUE(abs(del) > eps)) {
    t0 <- abs(t0)
    del <- (sum((y - mu)^2 / (mu + mu^2 / t0)) - dfr) /
      sum((y - mu)^2 / (mu + t0)^2)
    t0 <- t0 - del
  }
  c(estimate = t0,
    warning = as.numeric(t0 < 0 || it == limit || !is.finite(t0)))
}

## MLE first; on a flagged MLE, fall back to MM; on a flagged MM, keep the pilot.
estimate_theta_r <- function(y, mu, dfr, limit = 200L,
                             eps = sqrt(.Machine$double.eps)) {
  t0 <- nb_theta_pilot(y, mu)
  tryCatch({
    mle <- nb_theta_mle(t0, y, mu, limit, eps)
    est <- mle[["estimate"]]
    if (mle[["warning"]] > 0.5) {
      mm <- nb_theta_mm(t0, y, mu, dfr, limit, eps)
      est <- mm[["estimate"]]
      if (mm[["warning"]] > 0.5) est <- t0
    }
    est
  }, error = function(e) t0)
}

## Poisson mean + Pearson/MoM dispersion, used when glm.nb fails or warns.
## The Poisson fitted values are the conditional-mean MLE under both the
## Poisson and the NB log-linear model, so they serve as mu_hat directly.
nb_fit_poisson_fallback <- function(Y, Z_design) {
  Y <- as.numeric(Y)
  pois_fit <- suppressWarnings(
    stats::glm.fit(x = Z_design, y = Y, family = stats::poisson())
  )
  mu_hat <- as.numeric(pois_fit$fitted.values)
  theta_hat <- estimate_theta_r(
    y = Y, mu = mu_hat, dfr = pois_fit$df.residual
  )
  theta_hat <- max(1e-4, theta_hat)  # spacrt's lower truncation
  if (length(mu_hat) != length(Y) || any(!is.finite(mu_hat)) ||
      any(mu_hat <= 0) || is.na(theta_hat) || theta_hat <= 0) {
    stop("negative-binomial fallback fit returned invalid estimates")
  }
  list(mu_hat = mu_hat, theta_hat = theta_hat,
       dispersion_method = "poisson_pearson_mm")
}

## ---- NB null fit: glm.nb when it converges, robust fallback otherwise -------
nb_fit <- function(Y, Z_design) {
  Z_design <- as.matrix(Z_design)
  fit <- tryCatch(
    suppressWarnings(
      MASS::glm.nb(
        Y ~ Z_design - 1,
        control = stats::glm.control(maxit = 100L)
      )
    ),
    error = function(e) NULL
  )
  if (!is.null(fit) && is.null(fit$th.warn) && isTRUE(fit$converged)) {
    mu_hat <- as.numeric(stats::fitted(fit))
    theta_hat <- as.numeric(fit$theta)
    if (length(mu_hat) == length(Y) && all(is.finite(mu_hat)) &&
        all(mu_hat > 0) && length(theta_hat) == 1L &&
        is.finite(theta_hat) && theta_hat > 0) {
      return(list(mu_hat = mu_hat, theta_hat = theta_hat,
                  dispersion_method = "glm.nb"))
    }
  }
  ## glm.nb errored, hit its dispersion iteration cap, or returned an invalid
  ## fit: use the compiled-free Poisson + Pearson-moments fallback instead.
  nb_fit_poisson_fallback(Y, Z_design)
}

## A scalar numerator divided by the square root of a non-negative quadratic
## form.  This matches the degenerate-case convention in the Section 6.1 code.
score_ratio <- function(numerator, denominator_sq, tolerance = 1e-12) {
  if (!is.finite(denominator_sq) || denominator_sq < -tolerance) {
    stop("invalid score variance")
  }
  denominator_sq <- max(0, denominator_sq)
  if (denominator_sq > tolerance) {
    return(numerator / sqrt(denominator_sq))
  }
  if (abs(numerator) <= sqrt(.Machine$double.eps)) return(0)
  sign(numerator) * Inf
}

## Population-variance form retained as a compatibility helper.  Algebraically
## it is sum(R) / sqrt(sum((R - mean(R))^2)).
studentized_stat <- function(R) {
  centered_ss <- sum((R - mean(R))^2)
  score_ratio(sum(R), centered_ss)
}

## ---- Gene-level cache: quantities independent of the tested exposure X ----
precompute_score_cache <- function(Y, Z_design) {
  Z_design <- as.matrix(Z_design)
  nb_null <- nb_fit(Y, Z_design)
  mu_hat <- as.numeric(nb_null$mu_hat)
  phi_bar <- 1 / nb_null$theta_hat

  ## H is the derivative/variance factor and W is the working Fisher weight.
  r_hat <- (as.numeric(Y) - mu_hat) / (1 + phi_bar * mu_hat)
  w_hat <- mu_hat / (1 + phi_bar * mu_hat)
  if (any(!is.finite(r_hat)) || any(!is.finite(w_hat)) || any(w_hat <= 0)) {
    stop("invalid fitted negative-binomial score quantities")
  }

  ZWZ <- crossprod(Z_design, w_hat * Z_design)
  ## A near-zero-expression gene (e.g. KRTAP4-8: 15 nonzero counts / 4148
  ## cells) drives the working Fisher weights toward zero and makes the
  ## weighted nuisance information ZWZ numerically singular, so solve() aborts.
  ## Flag the gene as degenerate and use a Moore-Penrose inverse so the caller
  ## can resample it (or, if it cannot, fall back to the asymptotic score
  ## p-value).  Non-degenerate genes are unaffected: solve() succeeds and the
  ## flag is FALSE.
  ZWZ_inv <- tryCatch(solve(ZWZ), error = function(e) NULL)
  degenerate <- is.null(ZWZ_inv)
  if (degenerate) ZWZ_inv <- MASS::ginv(ZWZ)

  weighted_design <- sqrt(w_hat) * Z_design
  weighted_qr <- qr(weighted_design)
  if (weighted_qr$rank < ncol(Z_design)) degenerate <- TRUE
  nuisance_basis <- qr.Q(weighted_qr, complete = FALSE)[,
    seq_len(ncol(Z_design)), drop = FALSE]

  list(
    Z_design = Z_design,
    theta_hat = nb_null$theta_hat,
    dispersion_method = nb_null$dispersion_method,
    r_hat = r_hat,
    w_hat = w_hat,
    ZWZ_inv = ZWZ_inv,
    nuisance_basis = nuisance_basis,
    degenerate = degenerate
  )
}

## Recompute the W-residualized exposure and efficient-score contributions.
## This is called once for the observed X and again for each label permutation.
efficient_score_components <- function(X, score_cache) {
  X <- as.numeric(X)
  if (length(X) != length(score_cache$r_hat)) {
    stop("X and score cache have incompatible lengths")
  }

  projection_coef <- score_cache$ZWZ_inv %*%
    crossprod(score_cache$Z_design, score_cache$w_hat * X)
  X_tilde <- as.numeric(X - score_cache$Z_design %*% projection_coef)
  score_contrib <- X_tilde * score_cache$r_hat
  score_sum <- sum(score_contrib)
  information_denom_sq <- sum(score_cache$w_hat * X_tilde^2)
  empirical_denom_sq <- sum((score_contrib - mean(score_contrib))^2)

  list(
    X_tilde = X_tilde,
    score_contrib = score_contrib,
    score_sum = score_sum,
    information_denom_sq = information_denom_sq,
    empirical_denom_sq = empirical_denom_sq,
    information_stat = score_ratio(score_sum, information_denom_sq),
    empirical_stat = score_ratio(score_sum, empirical_denom_sq)
  )
}

## Pair-level cache: observed score plus the low-rank correction used by the
## standardized sign-flip reference.
precompute_score_pair <- function(X, score_cache) {
  components <- efficient_score_components(X, score_cache)
  b <- sqrt(score_cache$w_hat) * components$X_tilde
  b_norm_sq <- sum(b^2)
  score_norm_sq <- sum(components$score_contrib^2)
  degenerate <- !is.finite(b_norm_sq) || !is.finite(score_norm_sq) ||
    b_norm_sq <= .Machine$double.eps ||
    score_norm_sq <= .Machine$double.eps

  if (degenerate) {
    A <- matrix(
      0,
      nrow = length(X),
      ncol = ncol(score_cache$nuisance_basis)
    )
    score_eff_n <- 0
    max_score_share <- 0
  } else {
    A <- score_cache$nuisance_basis * as.numeric(b / sqrt(b_norm_sq))
    fourth_sum <- sum(components$score_contrib^4)
    score_eff_n <- if (is.finite(fourth_sum) && fourth_sum > 0) {
      score_norm_sq^2 / fourth_sum
    } else {
      0
    }
    max_score_share <- max(components$score_contrib^2) / score_norm_sq
  }

  c(components, list(
    A = A,
    degenerate = degenerate,
    score_eff_n = as.numeric(score_eff_n),
    max_score_share = as.numeric(max_score_share)
  ))
}

## ---- Statistics used by the four resampling callbacks ----
information_score_stat <- function(X, score_cache) {
  efficient_score_components(X, score_cache)$information_stat
}

empirical_score_stat <- function(X, score_cache) {
  efficient_score_components(X, score_cache)$empirical_stat
}

hgf_signflip_stat <- function(signs, pair_cache) {
  sum(as.numeric(signs) * pair_cache$score_contrib)
}

standardized_signflip_stat <- function(signs, pair_cache) {
  signs <- as.numeric(signs)
  numerator <- sum(signs * pair_cache$score_contrib)
  nuisance_projection <- as.numeric(crossprod(signs, pair_cache$A))
  denominator_sq <- 1 - sum(nuisance_projection^2)
  if (denominator_sq < 0 && denominator_sq > -1e-10) denominator_sq <- 0
  score_ratio(numerator, denominator_sq)
}

rademacher_signs <- function(n) {
  sample(c(-1, 1), size = n, replace = TRUE)
}

## Backward-compatible wrappers used by the focused checks.
nb_score_stat_fast <- function(X, score_cache) {
  information_score_stat(X, score_cache)
}
