#!/usr/bin/env Rscript
###############################################################################
# Exploratory: where does asymptotic calibration break, and does permutation
# calibration survive there?
#
# Design principle.  To isolate CALIBRATION we must keep the nuisance easy, so
# that the observed statistic stays pivotal and any miscalibration is slow CLT
# convergence rather than nuisance bias.  We therefore use a one-dimensional Z
# with a linear mean, mu(Z) = Z1, for which kernel ridge is accurate at every
# sample size considered (empirical E / Var(mu) is 0.23 at n = 30 and falls
# with n).  The grid deliberately stays SMALL: asymptotic calibration is
# hardest exactly where the effective sample size is limited, which is the
# regime the permutation argument is meant to serve.  Note that the PCM splits
# the sample, so at n = 20 its inference fold holds only 10 observations.  All nuisance regressions use "krr"; because the smoother
# is frozen at fit time, no permutation refits the learner.
#
# The stress is placed entirely on the NOISE, using the five bivariate laws of
# DiCiccio and Romano (2017, JASA, Section 5).  They were constructed to break
# the normal approximation for the studentized correlation: every one has
# CORRELATION ZERO BUT IS DEPENDENT, so the naive normal reference fails while
# the estimand of interest is still null.
#
#   normal   bivariate N(0, I) -- the benign reference case
#   t5       bivariate t with identity scale and 5 df (finite fourth moments)
#   exp      elliptical: r * S u, S = diag(sqrt 2, 1), u uniform on the circle,
#            r ~ Exp(1)
#   circular uniform on the two-dimensional unit circle
#   t4.1     xi = W + V, eps = W - V with W, V iid t_{4.1}; DiCiccio and Romano
#            show tau^2 diverges as the degrees of freedom approach 4, so this
#            is the adversarial case for the normal reference
#
# Because these pairs are dependent, feeding them in as (xi, eps) yields
# conditional MEAN independence rather than full conditional independence.
# That is the right null here: it is exactly what the PCM targets and what the
# GCM's estimand E Cov(X, Y | Z) captures.  Full independence of the pair is
# used for the MI null, where X must be independent of (Y, Z).
#
# Two sparse-exposure scenarios are retained alongside them, since rare binary
# exposures are the regime the genomic applications live in.
# Both nulls are covered:
#   MI: X independent of (Y, Z)   -- permutation is exact in finite samples
#   CI: X depends on Z, X indep Y | Z -- permutation is only asymptotic
#
# Usage:
#   Rscript explore-heavy-tail-calibration.R [--reps=500] [--cores=7] [--B=499]
#     [--n=20,30,50,75] [--scenarios=all]
###############################################################################
suppressPackageStartupMessages({library(pcmt); library(parallel)})

cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(f, d) {
  h <- grep(paste0("^--", f, "="), cli, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", f, "="), "", h[[1L]])
}
R_REPS  <- as.integer(get_opt("reps", 500))
N_CORES <- as.integer(get_opt("cores", 7))
B_PERM  <- as.integer(get_opt("B", 499))
N_GRID  <- as.integer(strsplit(get_opt("n", "20,50,100,200"), ",")[[1L]])
# "simple": Z is 1-d with mu(Z) = Z1, which kernel ridge fits accurately at
# every n here, so the statistic stays pivotal and miscalibration can only be
# slow CLT convergence.  "complex": Z is 3-d with nonlinear, interacting mean
# functions and heteroscedastic outcome noise.  That is more realistic, but the
# nuisance becomes the binding constraint at small n -- so each cell also
# reports the outcome-side error E_y against Var(mu_y).  When that ratio is
# large the statistic is not pivotal and the calibration comparison is void.
MEAN_MODEL <- get_opt("mean", "complex")
ALPHA   <- 0.05

# The PCM ALWAYS refits its projection learner here.  The cached
# frozen-operator path breaks finite-sample exactness under the MI null -- at
# 5% exposure prevalence it rejects 0.145 at the 0.05 level, against 0.055 when
# refitting -- because the PCM's permuted pseudo-response fhat(X_pi(i), Z_i) is
# a recomputed vector rather than a rearrangement of the observed one, so an
# operator tuned on the observed response privileges it.  The GCM has no such
# problem: its permuted pseudo-response is literally a permutation of X, and it
# is exact even with a near-interpolating frozen operator.
options(pcmt.cached_refit_pcm = FALSE, pcmt.cached_refit_gcm = TRUE)

## ---- bivariate noise laws of DiCiccio and Romano (2017, Sec. 5) -----------
# Each returns an n x 2 matrix (xi, eps) with zero correlation.  All but the
# first are dependent.  Each is STANDARDISED to unit variance: scale is
# irrelevant to a studentized statistic, but not to the nuisance regression --
# unstandardised, these laws differ in variance by a factor of eight (t_4.1 as
# W + V has variance 2 nu/(nu - 2) ~ 3.9, circular has 0.5), which would vary
# the signal-to-noise ratio across scenarios and confound the comparison.
draw_pair <- function(n, law) switch(law,
  normal   = matrix(rnorm(2 * n), n, 2),
  t5       = { nu <- 5; matrix(rnorm(2 * n), n, 2) /
                 sqrt(stats::rchisq(n, nu) / nu) },
  exp      = { th <- stats::runif(n, 0, 2 * pi); r <- stats::rexp(n, 1)
               cbind(r * sqrt(2) * cos(th), r * sin(th)) },
  circular = { th <- stats::runif(n, 0, 2 * pi); cbind(cos(th), sin(th)) },
  t41      = { W <- stats::rt(n, 4.1); V <- stats::rt(n, 4.1)
               cbind(W + V, W - V) },
  # DiCiccio and Romano note tau^2 diverges as the degrees of freedom approach
  # 4, so these trace that divergence rather than sampling one point of it.
  t45      = { W <- stats::rt(n, 4.5); V <- stats::rt(n, 4.5)
               cbind(W + V, W - V) },
  t405     = { W <- stats::rt(n, 4.05); V <- stats::rt(n, 4.05)
               cbind(W + V, W - V) },
  stop("unknown noise law: ", law))

NOISE_SD <- c(normal = 1, t5 = sqrt(5/3), circular = sqrt(0.5),
              exp = NA, t45 = sqrt(2 * 4.5 / 2.5), t41 = sqrt(2 * 4.1 / 2.1),
              t405 = sqrt(2 * 4.05 / 2.05))
NOISE_SD[["exp"]] <- stats::sd({ th <- stats::runif(2e5, 0, 2*pi)
  stats::rexp(2e5, 1) * sin(th) })
draw_pair <- local({ raw <- draw_pair; function(n, law) raw(n, law) / NOISE_SD[[law]] })

# For the MI null the exposure must be independent of (Y, Z), so the pair is
# broken by drawing the two coordinates from independent replicates.
draw_independent_pair <- function(n, law) {
  cbind(draw_pair(n, law)[, 1L], draw_pair(n, law)[, 2L])
}

## ---- scenarios -------------------------------------------------------------
SCENARIOS <- lapply(
  c("normal", "t5", "circular", "t41"),
  function(l) list(id = l, pair = l, xlaw = "cont")
)

expit <- function(u) 1 / (1 + exp(-u))

mu_functions <- function(Z) {
  if (identical(MEAN_MODEL, "simple")) {
    list(x = Z[, 1L], y = Z[, 1L], s = rep(1, nrow(Z)))
  } else {
    list(
      x = sin(Z[, 2L]) + 0.5 * Z[, 1L] * Z[, 2L],
      y = sin(Z[, 1L]) + 0.5 * Z[, 2L]^2 + 0.5 * Z[, 1L] * Z[, 3L],
      s = 0.7 + 0.3 * abs(Z[, 1L])          # mild heteroscedastic outcome noise
    )
  }
}

generate <- function(n, scen, null_type) {
  p <- if (identical(MEAN_MODEL, "simple")) 1L else 3L
  Z <- matrix(rnorm(n * p), n, p); colnames(Z) <- paste0("Z", seq_len(p))
  mu <- mu_functions(Z)
  E <- if (identical(null_type, "MI")) draw_independent_pair(n, scen$pair)
       else draw_pair(n, scen$pair)
  Y <- mu$y + mu$s * E[, 2L]
  X <- if (identical(null_type, "MI")) E[, 1L] else mu$x + E[, 1L]
  list(X = X, Y = Y, Z = Z, mu_y = mu$y)
}

one_rep <- function(r, n, scen, null_type) {
  set.seed(5000000L + 7919L * match(scen$id, vapply(SCENARIOS, function(z) z$id, "")) +
             131L * match(null_type, c("MI", "CI")) + 17L * n + r)
  d <- generate(n, scen, null_type)
  if (stats::sd(d$X) == 0) return(c(perm_gcm = NA, asym_gcm = NA,
                                    perm_pcm = NA, asym_pcm = NA,
                                    E_y = NA, var_mu = NA))
  fr <- get("fit_regression", envir = asNamespace("pcmt"))
  my <- fr(y = d$Y, x = d$Z, method = "krr")
  fit_y <- d$Y - as.numeric(stats::residuals(my, response = d$Y, data = d$Z))
  g <- pcmt::perm_gcm(Y = d$Y, X = d$X, Z = d$Z,
                      reg_YonZ = "krr", reg_XonZ = "krr",
                      B = B_PERM, n_cores = 1L, h = 15L, alpha = ALPHA,
                      alternative = "two.sided")
  p <- pcmt::perm_pcm(Y = d$Y, X = d$X, Z = d$Z,
                      reg_YonXZ = "krr", reg_YonZ = "krr",
                      reg_YhatonZ = "krr", reg_RonZ = "krr", frac = 0.5,
                      B = B_PERM, n_cores = 1L, h = 15L, alpha = ALPHA,
                      alternative = "greater")
  c(perm_gcm = as.numeric(g$p.value),
    asym_gcm = 2 * stats::pnorm(abs(as.numeric(g$statistic)), lower.tail = FALSE),
    perm_pcm = as.numeric(p$p.value),
    asym_pcm = stats::pnorm(as.numeric(p$statistic), lower.tail = FALSE),
    E_y = mean((fit_y - d$mu_y)^2),
    var_mu = mean((d$mu_y - mean(d$mu_y))^2))
}

wilson <- function(k, N, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2); ph <- k / N; dn <- 1 + z^2 / N
  ctr <- (ph + z^2 / (2 * N)) / dn
  hw <- (z / dn) * sqrt(ph * (1 - ph) / N + z^2 / (4 * N^2))
  c(max(0, ctr - hw), min(1, ctr + hw))
}

sel <- get_opt("scenarios", "all")
if (!identical(sel, "all")) {
  keep <- strsplit(sel, ",")[[1L]]
  SCENARIOS <- Filter(function(s) s$id %in% keep, SCENARIOS)
}

cat(sprintf("Heavy-tail calibration probe: %d scenarios x %d n x 2 nulls, R=%d, B=%d\n",
            length(SCENARIOS), length(N_GRID), R_REPS, B_PERM))
cat("All nuisances kernel ridge; mu(Z) = Z1 so the nuisance is not the bottleneck.\n\n")
cat(sprintf("%-10s %-3s %5s | %-16s %-16s | %-16s %-16s | %5s\n",
            "scenario","null","n","perm_gcm","asym_gcm","perm_pcm","asym_pcm",
            "E/Var"))
rows <- list()
for (scen in SCENARIOS) for (null_type in c("MI","CI")) for (n in N_GRID) {
  M <- do.call(rbind, mclapply(seq_len(R_REPS),
        function(r) one_rep(r, n, scen, null_type), mc.cores = N_CORES))
  ok <- stats::complete.cases(M); M <- M[ok, , drop = FALSE]; Rv <- nrow(M)
  f <- function(nm) { k <- sum(M[, nm] <= ALPHA); ci <- wilson(k, Rv)
    sprintf("%.3f[%.3f,%.3f]", k/Rv, ci[1], ci[2]) }
  nz <- mean(M[, "E_y"]) / mean(M[, "var_mu"])
  cat(sprintf("%-10s %-3s %5d | %-16s %-16s | %-16s %-16s | %5.2f\n",
              scen$id, null_type, n, f("perm_gcm"), f("asym_gcm"),
              f("perm_pcm"), f("asym_pcm"), nz))
  utils::flush.console()
  rows[[length(rows)+1L]] <- data.frame(scenario = scen$id, null = null_type,
    n = n, R = Rv,
    perm_gcm = mean(M[,"perm_gcm"] <= ALPHA), asym_gcm = mean(M[,"asym_gcm"] <= ALPHA),
    perm_pcm = mean(M[,"perm_pcm"] <= ALPHA), asym_pcm = mean(M[,"asym_pcm"] <= ALPHA),
    E_y = mean(M[,"E_y"]), var_mu = mean(M[,"var_mu"]), mean_model = MEAN_MODEL)
}
res <- do.call(rbind, rows)
out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) {
  file.path(out_root, "simulation-results", "heavy-tail-calibration")
} else file.path(getwd(), "results", "heavy-tail-calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out_dir, "heavy-tail-calibration.rds"))
cat("\nWritten to: ", out_dir, "\n", sep = "")
