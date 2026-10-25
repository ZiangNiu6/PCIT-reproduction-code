# Permutation calibration for the sparse binary-exposure CI design.
#
# A raw-label permutation of binary X with m exposed observations is a
# uniformly distributed subset of m indices.  Sampling that subset directly
# avoids forming n-by-B matrices and uses one subset for all available score
# statistics.  The permutation law is unchanged, but the RNG stream (and hence
# individual p-values at a fixed seed) differs from sample.int(n) shuffles.

.sparse_score_context <- function(cache, Z, empirical = FALSE) {
  n <- length(Z)
  if (is.null(cache$D) || nrow(cache$D) != n ||
      ncol(cache$D) != 2L ||
      !isTRUE(all(cache$D[, 1L] == 1)) ||
      !isTRUE(all(cache$D[, 2L] == Z)) ||
      !is.matrix(cache$gram_inverse) ||
      !identical(dim(cache$gram_inverse), c(2L, 2L)) ||
      length(cache$W) != n || length(cache$outcome_factor) != n) {
    stop("This engine requires the intercept-plus-Z score cache.",
         call. = FALSE)
  }
  U <- as.numeric(cache$outcome_factor)
  context <- list(
    n = n,
    Z = Z,
    W = as.numeric(cache$W),
    U = U,
    inverse = cache$gram_inverse,
    total_U = sum(U),
    total_ZU = sum(Z * U),
    total_Z = sum(Z),
    total_Z2 = sum(Z^2),
    empirical = empirical
  )
  if (empirical) {
    U2 <- U^2
    context$total_U2 <- sum(U2)
    context$total_ZU2 <- sum(Z * U2)
    context$total_Z2U2 <- sum(Z^2 * U2)
  }
  context
}

.sparse_reference_statistics <- function(context, exposed_indices) {
  Z <- context$Z[exposed_indices]
  U <- context$U[exposed_indices]
  W <- context$W[exposed_indices]
  m <- length(exposed_indices)

  weighted_0 <- sum(W)
  weighted_1 <- sum(W * Z)
  inverse <- context$inverse
  a0 <- inverse[1L, 1L] * weighted_0 +
    inverse[1L, 2L] * weighted_1
  a1 <- inverse[2L, 1L] * weighted_0 +
    inverse[2L, 2L] * weighted_1

  # sum_i (X_i - a0 - a1 Z_i) U_i
  numerator <- sum(U) - a0 * context$total_U -
    a1 * context$total_ZU

  # Weighted projection identity: sum_i W_i Xtilde_i^2 =
  # sum_{i exposed} W_i - (a0, a1) dot (sum_E W, sum_E W Z).
  information_ss <- weighted_0 - a0 * weighted_0 - a1 * weighted_1
  # The original statistic also uses sum Xtilde^2 in its tolerance guard.
  unweighted_ss <- m - 2 * a0 * m - 2 * a1 * sum(Z) +
    a0^2 * context$n + 2 * a0 * a1 * context$total_Z +
    a1^2 * context$total_Z2
  information_scale <- max(unweighted_ss, 1)
  information_positive <- information_ss >
    100 * .Machine$double.eps * information_scale
  information_statistic <- if (information_positive) {
    numerator / sqrt(information_ss)
  } else if (abs(numerator) >
             sqrt(.Machine$double.eps * information_scale)) {
    sign(numerator) * Inf
  } else {
    0
  }

  empirical_statistic <- NA_real_
  if (isTRUE(context$empirical)) {
    # sum_i [(X_i - a0 - a1 Z_i) U_i]^2, expanded around X_i = 0.
    uncentered_ss <- a0^2 * context$total_U2 +
      2 * a0 * a1 * context$total_ZU2 +
      a1^2 * context$total_Z2U2 +
      sum(U^2 * (1 - 2 * a0 - 2 * a1 * Z))
    centered_ss <- uncentered_ss - numerator^2 / context$n
    empirical_scale <- max(uncentered_ss, 1)
    empirical_positive <- centered_ss >
      100 * .Machine$double.eps * empirical_scale
    empirical_statistic <- if (empirical_positive) {
      numerator / sqrt(centered_ss)
    } else if (abs(numerator) >
               sqrt(.Machine$double.eps * empirical_scale)) {
      sign(numerator) * Inf
    } else {
      0
    }
  }
  c(information = information_statistic, empirical = empirical_statistic)
}

.sparse_observed_statistics <- function(cache, empirical = FALSE) {
  if (isTRUE(cache$degenerate)) {
    return(c(information = 0, empirical = if (empirical) 0 else NA_real_))
  }
  score <- cache$score_contrib
  numerator <- sum(score)
  information_ss <- sum(cache$X_tilde^2 * cache$W)
  information_scale <- max(sum(cache$X_tilde^2), 1)
  information_statistic <- if (
    information_ss > 100 * .Machine$double.eps * information_scale
  ) {
    numerator / sqrt(information_ss)
  } else if (abs(numerator) >
             sqrt(.Machine$double.eps * information_scale)) {
    sign(numerator) * Inf
  } else {
    0
  }
  empirical_statistic <- NA_real_
  if (empirical) {
    uncentered_ss <- sum(score^2)
    centered_ss <- uncentered_ss - numerator^2 / length(score)
    empirical_scale <- max(uncentered_ss, 1)
    empirical_statistic <- if (
      centered_ss > 100 * .Machine$double.eps * empirical_scale
    ) {
      numerator / sqrt(centered_ss)
    } else if (abs(numerator) >
               sqrt(.Machine$double.eps * empirical_scale)) {
      sign(numerator) * Inf
    } else {
      0
    }
  }
  c(information = information_statistic, empirical = empirical_statistic)
}

# data needs X, Z, n, perm_seed and (when available) score_cache. The optional
# poisson_cache has the same cache fields as prepare_poisson_only_cache in
# search-nb-sparse-ci-null.R. Passing NULL for either cache skips that method.
# Returned names match the search script's replicate columns.
sparse_binary_permutation_p_values <- function(
    data, permutations, poisson_cache = NULL) {
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 1 || permutations != as.integer(permutations)) {
    stop("permutations must be a positive integer.", call. = FALSE)
  }
  permutations <- as.integer(permutations)
  n <- as.integer(data$n)
  if (length(n) != 1L || is.na(n) || n < 1L ||
      length(data$X) != n || length(data$Z) != n ||
      anyNA(data$X) || any(!(data$X %in% c(0, 1))) ||
      length(data$perm_seed) != 1L || is.na(data$perm_seed)) {
    stop("Invalid binary-exposure data or permutation seed.", call. = FALSE)
  }
  if (is.null(data$score_cache) && is.null(poisson_cache)) {
    stop("At least one score cache is required.", call. = FALSE)
  }

  observed_nb <- if (is.null(data$score_cache)) {
    c(information = NA_real_, empirical = NA_real_)
  } else {
    .sparse_observed_statistics(data$score_cache, empirical = TRUE)
  }
  observed_poisson <- if (is.null(poisson_cache)) {
    c(information = NA_real_, empirical = NA_real_)
  } else {
    .sparse_observed_statistics(poisson_cache, empirical = FALSE)
  }
  nb_active <- !is.null(data$score_cache) &&
    !isTRUE(data$score_cache$degenerate)
  poisson_active <- !is.null(poisson_cache) &&
    !isTRUE(poisson_cache$degenerate)
  nb_context <- if (nb_active) {
    .sparse_score_context(data$score_cache, data$Z, empirical = TRUE)
  } else {
    NULL
  }
  poisson_context <- if (poisson_active) {
    .sparse_score_context(poisson_cache, data$Z, empirical = FALSE)
  } else {
    NULL
  }

  nb_information_hits <- nb_empirical_hits <- poisson_hits <- 0L
  if (nb_active || poisson_active) {
    set.seed(as.integer(data$perm_seed))
    m <- sum(data$X)
    for (b in seq_len(permutations)) {
      exposed_indices <- sample.int(n, m)
      if (nb_active) {
        reference <- .sparse_reference_statistics(
          nb_context, exposed_indices
        )
        nb_information_hits <- nb_information_hits +
          (abs(reference[["information"]]) >=
             abs(observed_nb[["information"]]))
        nb_empirical_hits <- nb_empirical_hits +
          (abs(reference[["empirical"]]) >=
             abs(observed_nb[["empirical"]]))
      }
      if (poisson_active) {
        reference <- .sparse_reference_statistics(
          poisson_context, exposed_indices
        )
        poisson_hits <- poisson_hits +
          (abs(reference[["information"]]) >=
             abs(observed_poisson[["information"]]))
      }
    }
  }
  if (!is.null(data$score_cache) && !nb_active) {
    nb_information_hits <- nb_empirical_hits <- permutations
  }
  if (!is.null(poisson_cache) && !poisson_active) {
    poisson_hits <- permutations
  }
  list(
    information_statistic = observed_nb[["information"]],
    empirical_statistic = observed_nb[["empirical"]],
    information_hits = if (is.null(data$score_cache)) NA_integer_ else
      nb_information_hits,
    empirical_hits = if (is.null(data$score_cache)) NA_integer_ else
      nb_empirical_hits,
    p_barry = if (is.null(data$score_cache)) NA_real_ else
      (1 + nb_information_hits) / (permutations + 1),
    p_ours = if (is.null(data$score_cache)) NA_real_ else
      (1 + nb_empirical_hits) / (permutations + 1),
    poisson_only_statistic = observed_poisson[["information"]],
    poisson_only_hits = if (is.null(poisson_cache)) NA_integer_ else
      poisson_hits,
    p_poisson_only = if (is.null(poisson_cache)) NA_real_ else
      (1 + poisson_hits) / (permutations + 1)
  )
}

# Compare the sparse algebra and hit counts with a direct dense calculation
# under *the same sampled subsets*. Returns the maximum statistic difference.
# It is intended as a one-time check after changing the score implementation.
validate_sparse_binary_permutation <- function(seed = 37L, permutations = 99L) {
  set.seed(seed)
  n <- 250L
  Z <- stats::rnorm(n)
  X <- integer(n)
  X[sample.int(n, 12L)] <- 1L
  mu <- 4 * exp(0.6 * Z)
  Y <- stats::rnbinom(n, size = 0.5, mu = mu)
  D <- cbind(`(Intercept)` = 1, Z = Z)
  make_cache <- function(W, U) {
    inverse <- solve(crossprod(D, D * W))
    a <- inverse %*% crossprod(D, W * X)
    x_tilde <- as.numeric(X - D %*% a)
    score <- x_tilde * U
    list(
      D = D, W = W, gram_inverse = inverse, outcome_factor = U,
      X_tilde = x_tilde, score_contrib = score,
      degenerate = FALSE
    )
  }
  nb_cache <- make_cache(mu / (1 + mu / 0.5),
                         (Y - mu) / (1 + mu / 0.5))
  poisson_cache <- make_cache(mu, Y - mu)
  data <- list(
    n = n, X = X, Z = Z, perm_seed = sample.int(.Machine$integer.max, 1L),
    score_cache = nb_cache
  )
  contexts <- list(
    nb = .sparse_score_context(nb_cache, Z, TRUE),
    poisson = .sparse_score_context(poisson_cache, Z, FALSE)
  )
  caches <- list(nb = nb_cache, poisson = poisson_cache)
  max_difference <- 0
  dense_hits <- c(nb_information = 0L, nb_empirical = 0L,
                  poisson_information = 0L)
  observed_nb <- .sparse_observed_statistics(nb_cache, TRUE)
  observed_poisson <- .sparse_observed_statistics(poisson_cache, FALSE)
  set.seed(data$perm_seed)
  for (b in seq_len(permutations)) {
    indices <- sample.int(n, sum(X))
    for (name in names(caches)) {
      cache <- caches[[name]]
      x_perm <- integer(n)
      x_perm[indices] <- 1L
      a <- cache$gram_inverse %*% crossprod(D, cache$W * x_perm)
      x_tilde <- as.numeric(x_perm - D %*% a)
      score <- x_tilde * cache$outcome_factor
      numerator <- sum(score)
      den_i <- sum(x_tilde^2 * cache$W)
      scale_i <- max(sum(x_tilde^2), 1)
      dense_i <- if (den_i > 100 * .Machine$double.eps * scale_i) {
        numerator / sqrt(den_i)
      } else if (abs(numerator) > sqrt(.Machine$double.eps * scale_i)) {
        sign(numerator) * Inf
      } else {
        0
      }
      sparse <- .sparse_reference_statistics(contexts[[name]], indices)
      max_difference <- max(max_difference,
                            abs(dense_i - sparse[["information"]]))
      if (name == "nb") {
        den_e_raw <- sum(score^2)
        den_e <- den_e_raw - numerator^2 / n
        scale_e <- max(den_e_raw, 1)
        dense_e <- if (den_e > 100 * .Machine$double.eps * scale_e) {
          numerator / sqrt(den_e)
        } else if (abs(numerator) > sqrt(.Machine$double.eps * scale_e)) {
          sign(numerator) * Inf
        } else {
          0
        }
        max_difference <- max(max_difference,
                              abs(dense_e - sparse[["empirical"]]))
        dense_hits[["nb_information"]] <- dense_hits[["nb_information"]] +
          (abs(dense_i) >= abs(observed_nb[["information"]]))
        dense_hits[["nb_empirical"]] <- dense_hits[["nb_empirical"]] +
          (abs(dense_e) >= abs(observed_nb[["empirical"]]))
      } else {
        dense_hits[["poisson_information"]] <-
          dense_hits[["poisson_information"]] +
          (abs(dense_i) >= abs(observed_poisson[["information"]]))
      }
    }
  }
  sparse <- sparse_binary_permutation_p_values(data, permutations,
                                               poisson_cache)
  stopifnot(
    is.finite(max_difference), max_difference < 1e-10,
    identical(as.integer(sparse$information_hits),
              as.integer(dense_hits[["nb_information"]])),
    identical(as.integer(sparse$empirical_hits),
              as.integer(dense_hits[["nb_empirical"]])),
    identical(as.integer(sparse$poisson_only_hits),
              as.integer(dense_hits[["poisson_information"]]))
  )
  invisible(list(max_statistic_difference = max_difference,
                 hit_counts = dense_hits))
}
