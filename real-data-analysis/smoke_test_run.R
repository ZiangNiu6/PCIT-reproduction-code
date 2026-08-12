#!/usr/bin/env Rscript
## smoke_test_run.R
##
## Fast, cluster-independent contract test for the five Perturb-seq score
## calibrations. It uses a small synthetic NB dataset and checks the shared
## algebra; no ondisc data or adaptive-FDR package is required.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "method_config.R"))
source(file.path(script_dir, "score_helpers.R"))

expect_near <- function(x, y, tolerance, label) {
  error <- max(abs(x - y))
  if (!is.finite(error) || error > tolerance) {
    stop(sprintf("%s failed: max error %.3g > %.3g", label, error, tolerance))
  }
}

## The paired-draw seed map must be collision-free over distinct coordinates
## in the configured block, while deliberately agreeing across paired methods.
seed_grid <- expand.grid(
  rep_index = 1:2,
  q = c(0.01, 0.05),
  stream = 1:2,
  iteration = 1:3,
  pair_idx = 1:4
)
seed_grid$seed <- mapply(
  perturb_resample_seed,
  rep_index = seed_grid$rep_index,
  q = seed_grid$q,
  iteration = seed_grid$iteration,
  pair_idx = seed_grid$pair_idx,
  stream = seed_grid$stream,
  MoreArgs = list(
    max_reps = 20L,
    q_levels = c(0.01, 0.05),
    n_pairs = 1000L,
    B_max = 25000L
  )
)
stopifnot(!anyDuplicated(seed_grid$seed))

## Synthetic-null response genes must exclude every positive-control target,
## even when a legacy preprocessing object includes them in its candidate list.
test_gene_list <- c("POS1", "POS2", "NULL1", "NULL2", "NULL3")
test_null_pool <- perturb_null_gene_pool(
  gene_list = test_gene_list,
  negative_gene_ids = test_gene_list,
  positive_targets = c("POS1", "POS2"),
  expected_size = 3L
)
stopifnot(
  identical(test_null_pool, c("NULL1", "NULL2", "NULL3")),
  !any(test_null_pool %in% c("POS1", "POS2"))
)

set.seed(20260810)
n <- 120L
z <- rnorm(n)
batch <- factor(sample(letters[1:3], n, replace = TRUE))
Z_design <- model.matrix(~ z + batch)
X <- rbinom(n, 1L, plogis(-1 + 0.45 * z))
mu <- exp(1.1 + 0.35 * z)
Y <- rnbinom(n, mu = mu, size = 1.4)

gene_cache <- precompute_score_cache(Y, Z_design)
pair_cache <- precompute_score_pair(X, gene_cache)

## Weighted residualization and null-score identities.
expect_near(
  as.numeric(crossprod(Z_design, gene_cache$w_hat * pair_cache$X_tilde)),
  rep(0, ncol(Z_design)),
  1e-8,
  "weighted X residualization"
)
score_equation_scale <- max(1, sum(abs(gene_cache$r_hat)))
expect_near(
  as.numeric(crossprod(Z_design, gene_cache$r_hat)) / score_equation_scale,
  rep(0, ncol(Z_design)),
  1e-5,
  "NB null score equation"
)
expect_near(
  sum(X * gene_cache$r_hat),
  pair_cache$score_sum,
  1e-6 * score_equation_scale,
  "common efficient-score numerator"
)

## Both score denominators and the former fast formula agree.
expect_near(
  pair_cache$information_denom_sq,
  sum(gene_cache$w_hat * pair_cache$X_tilde^2),
  1e-10,
  "information denominator"
)
expect_near(
  pair_cache$empirical_denom_sq,
  sum((pair_cache$score_contrib - mean(pair_cache$score_contrib))^2),
  1e-10,
  "empirical denominator"
)
expect_near(
  nb_score_stat_fast(X, gene_cache),
  pair_cache$information_stat,
  1e-12,
  "cached information statistic"
)

## The QR correction has the standardized-flip geometry used in Section 6.1.
Q <- gene_cache$nuisance_basis
expect_near(crossprod(Q), diag(ncol(Q)), 1e-10, "orthonormal nuisance basis")
b <- sqrt(gene_cache$w_hat) * pair_cache$X_tilde
expect_near(as.numeric(crossprod(Q, b)), rep(0, ncol(Q)), 1e-8,
            "score direction orthogonal to nuisance basis")
expect_near(colSums(pair_cache$A), rep(0, ncol(pair_cache$A)), 1e-8,
            "identity-flip correction")

## Paired permutation draws preserve exposure abundance and return both
## studentizations of the same permuted efficient score.
B <- 199L
permutations <- replicate(B, sample.int(n), simplify = FALSE)
T_info <- vapply(permutations, function(index) {
  stopifnot(sum(X[index]) == sum(X))
  information_score_stat(X[index], gene_cache)
}, numeric(1))
T_empirical <- vapply(permutations, function(index) {
  empirical_score_stat(X[index], gene_cache)
}, numeric(1))
stopifnot(all(is.finite(T_info)), all(is.finite(T_empirical)))

## Paired signs obey the projection bound and both references are odd in g.
sign_draws <- replicate(B, rademacher_signs(n), simplify = FALSE)
for (signs in sign_draws) {
  denominator_sq <- 1 - sum(as.numeric(crossprod(signs, pair_cache$A))^2)
  if (denominator_sq < -1e-10 || denominator_sq > 1 + 1e-10) {
    stop("standardized sign-flip denominator left [0, 1]")
  }
  expect_near(
    hgf_signflip_stat(-signs, pair_cache),
    -hgf_signflip_stat(signs, pair_cache),
    1e-10,
    "HGF sign symmetry"
  )
  expect_near(
    standardized_signflip_stat(-signs, pair_cache),
    -standardized_signflip_stat(signs, pair_cache),
    1e-10,
    "standardized-flip sign symmetry"
  )
}

## One fixed-B p-value from each calibration verifies the end-to-end contract.
T_hgf <- vapply(
  sign_draws, hgf_signflip_stat, numeric(1), pair_cache = pair_cache
)
T_standardized <- vapply(
  sign_draws, standardized_signflip_stat, numeric(1), pair_cache = pair_cache
)
p_values <- c(
  parametric_score = 2 * pnorm(-abs(pair_cache$information_stat)),
  information_studentized_permutation =
    (1 + sum(abs(T_info) >= abs(pair_cache$information_stat))) / (B + 1),
  studentized_permutation =
    (1 + sum(abs(T_empirical) >= abs(pair_cache$empirical_stat))) / (B + 1),
  hgf_signflip_mc =
    (1 + sum(abs(T_hgf) >= abs(pair_cache$score_sum))) / (B + 1),
  standardized_signflip =
    (1 + sum(abs(T_standardized) >= abs(pair_cache$score_sum))) / (B + 1)
)
stopifnot(
  identical(names(p_values), PERTURB_METHOD_LEVELS),
  all(is.finite(p_values)),
  all(p_values > 0 & p_values <= 1)
)

## A constant exposure is a genuine degenerate score and must be harmless.
degenerate <- precompute_score_pair(rep(0, n), gene_cache)
stopifnot(
  isTRUE(degenerate$degenerate),
  degenerate$information_stat == 0,
  degenerate$empirical_stat == 0,
  hgf_signflip_stat(rep(1, n), degenerate) == 0,
  standardized_signflip_stat(rep(1, n), degenerate) == 0
)

cat("Five-method Perturb-seq score smoke test passed.\n")
print(round(p_values, 4))
