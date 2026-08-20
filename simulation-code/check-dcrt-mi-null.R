#!/usr/bin/env Rscript

###############################################################################
# MI-null stress test for the dCRT of Niu et al. (2024, AoS).
#
# Under the marginal-independence null X |= (Y, Z), marginal permutation is
# exact in finite samples, so any departure from the nominal level by a
# model-X procedure is a pure calibration artifact.  This script searches the
# MI null for configurations where the dCRT is inflated or deflated.
#
# The mechanism under test is condition (23) of Niu et al.: the dCRT resamples
# Xtilde ~ Lhat(X | Z) and therefore needs the second moment of Lhat to be
# accurate.  With an in-sample random-forest fit and X |= Z, the truth is
# Var(X | Z) = Var(X), but the in-sample residual second moment sigma^2 is far
# smaller because the forest interpolates.  Two errors then act in opposite
# directions:
#
#   (i)  sigma^2 understates Var(X | Z)     -> resampling reference too narrow
#                                           -> inflation;
#   (ii) in-sample residuals are negatively correlated across observations, so
#        the true sd of the statistic is below the independent-draw resampling
#        sd                                -> deflation.
#
# Which one wins is an empirical question, so the script reports the realized
# sd of the observed statistic alongside the mean resampling sd.
#
# Design (all cells satisfy the MI null exactly):
#   Z = (Z1, ..., Zp) iid N(0, 1);  Y = sin(2 pi Z1) + eps;  X |= (Y, Z).
# Only Z1 drives Y; the remaining p - 1 coordinates are noise covariates that
# increase how badly the X-on-Z forest overfits.
#
# Methods (all share one pair of nuisance fits, so they differ only in
# calibration):
#   asym_gcm_rf  studentized statistic, N(0, 1) reference
#   ndcrt_rf     studentized statistic, model-X resampling  (Appendix A)
#   dcrt_rf      unnormalized statistic, model-X resampling (Algorithm 1)
#   perm_gcm_rf  studentized statistic, marginal permutation -- exact here;
#                run only with --with-perm because it refits per draw
#
# Usage:
#   Rscript --vanilla check-dcrt-mi-null.R [--reps=N] [--cores=K] [--M=M]
#                                          [--with-perm] [--quick]
###############################################################################

suppressPackageStartupMessages({
  library(pcmt)
  library(parallel)
})

cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  hit <- grep(paste0("^--", flag, "="), cli, value = TRUE)
  if (length(hit) == 0L) return(default)
  as.numeric(sub(paste0("^--", flag, "="), "", hit[[1L]]))
}
quick <- "--quick" %in% cli
with_perm <- "--with-perm" %in% cli
n_reps <- as.integer(get_opt("reps", if (quick) 200 else 2000))
n_cores <- as.integer(get_opt("cores", 4))
M_resamp <- as.integer(get_opt("M", 499))
B_perm <- as.integer(get_opt("B", 499))
alpha <- 0.05

fit_regression <- get("fit_regression", envir = asNamespace("pcmt"))

studentized_stat <- function(R) {
  n <- length(R)
  m <- mean(R)
  v <- mean(R^2) - m^2
  if (v <= 0) return(0)
  sqrt(n) * m / sqrt(v)
}
# Sided resampling p-values, matching pcmt:::perm_pvalue.
resample_pvalue <- function(T_obs, T_resamp) {
  M <- length(T_resamp)
  c(two   = (1 + sum(abs(T_resamp) >= abs(T_obs))) / (1 + M),
    right = (1 + sum(T_resamp >= T_obs)) / (1 + M),
    left  = (1 + sum(T_resamp <= T_obs)) / (1 + M))
}

## ------------------------------------------------------------------- design
# X is drawn from its marginal law, independently of (Y, Z), so the MI null
# holds by construction whatever draw_X does.
draw_X <- function(n, law) {
  switch(law,
    gauss = rnorm(n),
    t3 = rt(n, df = 3) / sqrt(3),            # unit variance, heavy tails
    bern10 = rbinom(n, 1L, 0.10),            # sparse binary
    bern02 = rbinom(n, 1L, 0.02),            # very sparse binary
    stop("unknown X law: ", law)
  )
}

generate_mi <- function(n, p, law) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(Z) <- paste0("Z", seq_len(p))
  Y <- sin(2 * pi * Z[, 1L]) + rnorm(n)
  list(X = draw_X(n, law), Y = Y, Z = Z)
}

rf_args <- function(seed) {
  list(seed = as.integer(seed), num.threads = 1L, num.trees = 500L)
}

one_replicate <- function(rep_id, n, p, law) {
  set.seed(7000000L + rep_id)
  d <- generate_mi(n, p, law)

  mY <- fit_regression(y = d$Y, x = d$Z, method = "rf", args = rf_args(53114L))
  rY <- as.numeric(residuals(mY, response = d$Y, data = d$Z))
  mX <- fit_regression(y = d$X, x = d$Z, method = "rf", args = rf_args(53115L))
  rX <- as.numeric(residuals(mX, response = d$X, data = d$Z))

  sigma_hat <- sqrt(mean(rX^2))
  T_dcrt <- sum(rX * rY) / sqrt(n)
  T_stud <- studentized_stat(rX * rY)

  set.seed(53120L + rep_id)
  E <- matrix(rnorm(n * M_resamp), nrow = n, ncol = M_resamp)
  R_res <- (sigma_hat * E) * rY
  cs <- colSums(R_res)
  m_res <- cs / n
  v_res <- colMeans(R_res^2) - m_res^2
  v_res[v_res <= 0] <- Inf

  perm_p <- NA_real_
  if (with_perm) {
    perm_p <- R.utils::withSeed(
      pcmt::perm_gcm(Y = d$Y, X = d$X, Z = d$Z,
                     reg_YonZ = "rf", reg_XonZ = "rf",
                     args_YonZ = rf_args(53114L), args_XonZ = rf_args(53115L),
                     B = B_perm, n_cores = 1L, h = 15L, alpha = alpha,
                     alternative = "two.sided")$p.value,
      seed = 53112L
    )
  }

  nd <- resample_pvalue(T_stud, sqrt(n) * m_res / sqrt(v_res))
  dc <- resample_pvalue(T_dcrt, cs / sqrt(n))

  c(asym_gcm_rf_two   = 2 * stats::pnorm(abs(T_stud), lower.tail = FALSE),
    asym_gcm_rf_right = stats::pnorm(T_stud, lower.tail = FALSE),
    asym_gcm_rf_left  = stats::pnorm(T_stud),
    stats::setNames(nd, paste0("ndcrt_rf_", names(nd))),
    stats::setNames(dc, paste0("dcrt_rf_", names(dc))),
    perm_gcm_rf_two = perm_p,
    # diagnostics for the second-moment condition
    T_dcrt = T_dcrt,
    resamp_sd = stats::sd(cs / sqrt(n)),
    sigma_hat = sigma_hat)
}

## --------------------------------------------------------------------- grid
# --n, --p and --law take comma-separated lists, so a single extreme cell can
# be re-run on its own with the (expensive) permutation benchmark switched on.
get_list <- function(flag, default, as_int = TRUE) {
  hit <- grep(paste0("^--", flag, "="), cli, value = TRUE)
  if (length(hit) == 0L) return(default)
  v <- strsplit(sub(paste0("^--", flag, "="), "", hit[[1L]]), ",")[[1L]]
  if (as_int) as.integer(v) else v
}
grid <- expand.grid(
  n = get_list("n", if (quick) c(100L) else c(50L, 100L, 200L)),
  p = get_list("p", if (quick) c(1L, 20L) else c(1L, 5L, 20L)),
  law = get_list("law",
                 if (quick) c("gauss") else c("gauss", "t3", "bern10", "bern02"),
                 as_int = FALSE),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

cat(sprintf("MI-null dCRT stress test: %d cells x %d reps, M=%d, %d cores%s\n\n",
            nrow(grid), n_reps, M_resamp, n_cores,
            if (with_perm) ", with permutation GCM" else ""))

methods <- as.vector(t(outer(c("asym_gcm_rf", "ndcrt_rf", "dcrt_rf"),
                             c("two", "right", "left"), paste, sep = "_")))
methods <- c(methods, "perm_gcm_rf_two")
rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  t0 <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(n_reps),
                             function(r) one_replicate(r, g$n, g$p, g$law),
                             mc.cores = n_cores)
  bad <- !vapply(reps, is.numeric, logical(1L))
  if (any(bad)) stop(sprintf("%d replicate(s) failed in cell %d", sum(bad), i))
  R <- do.call(rbind, reps)

  rej <- colMeans(R[, methods, drop = FALSE] <= alpha)
  rows[[i]] <- data.frame(
    n = g$n, p = g$p, law = g$law,
    as.list(rej),
    # true null sd of the unnormalized statistic vs the sd the dCRT resamples
    true_sd = stats::sd(R[, "T_dcrt"]),
    resamp_sd = mean(R[, "resamp_sd"]),
    sigma_hat = mean(R[, "sigma_hat"]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rows[[i]]$sd_ratio <- rows[[i]]$resamp_sd / rows[[i]]$true_sd

  fmt3 <- function(m) sprintf("%.3f/%.3f/%.3f", rej[[paste0(m, "_two")]],
                              rej[[paste0(m, "_right")]], rej[[paste0(m, "_left")]])
  cat(sprintf("[%2d/%2d] n=%3d p=%2d %-7s (%.0fs)  two/right/left -- asym %s | ndcrt %s | dcrt %s%s  sd_ratio %.2f\n",
              i, nrow(grid), g$n, g$p, g$law, proc.time()[["elapsed"]] - t0,
              fmt3("asym_gcm_rf"), fmt3("ndcrt_rf"), fmt3("dcrt_rf"),
              if (with_perm) sprintf(" | perm %.3f", rej[["perm_gcm_rf_two"]]) else "",
              rows[[i]]$sd_ratio))
  utils::flush.console()
}

results <- do.call(rbind, rows)
mcse <- sqrt(alpha * (1 - alpha) / n_reps)

out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) {
  file.path(out_root, "simulation-results", "dcrt-check")
} else {
  file.path(getwd(), "results", "dcrt-check")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag <- if (with_perm) "-perm" else ""
saveRDS(results, file.path(out_dir, paste0("dcrt-mi-null-results", tag, ".rds")))
utils::write.csv(results,
                 file.path(out_dir, paste0("dcrt-mi-null-summary", tag, ".csv")),
                 row.names = FALSE)

cat(sprintf("\n=== MI-null rejection at alpha = %.2f (MCSE %.4f) ===\n",
            alpha, mcse))
sided_cols <- as.vector(t(outer(c("asym_gcm_rf", "ndcrt_rf", "dcrt_rf"),
                                c("two", "right", "left"), paste, sep = "_")))
show <- results[, c("n", "p", "law", sided_cols,
                    if (with_perm) "perm_gcm_rf_two", "sigma_hat", "sd_ratio")]
print(show, row.names = FALSE, digits = 3)

rank_cols <- c("n", "p", "law", "dcrt_rf_two", "dcrt_rf_right", "dcrt_rf_left",
               "ndcrt_rf_two", "ndcrt_rf_right", "ndcrt_rf_left", "sd_ratio")
cat("\nMost inflated cells (by two-sided dCRT):\n")
print(head(results[order(-results$dcrt_rf_two), rank_cols], 5),
      row.names = FALSE, digits = 3)
cat("\nMost deflated cells (by two-sided ndCRT):\n")
print(head(results[order(results$ndcrt_rf_two), rank_cols], 5),
      row.names = FALSE, digits = 3)
cat("\nWritten to: ", out_dir, "\n", sep = "")
