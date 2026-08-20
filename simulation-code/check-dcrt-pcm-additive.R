#!/usr/bin/env Rscript

###############################################################################
# Check simulation: does the dCRT of Niu et al. (2024, AoS) belong in the
# Section 5.2 comparison?
#
# This is an exploratory check, not a manuscript reproduction.  It reuses the
# Section 5.2 data-generating process (sim-spec-formal/pcm-additive.R) and adds
# two model-X procedures to the GCM family already compared there:
#
#   dcrt_rf   The \hat{dCRT} of Niu et al. (2024), Algorithm 1.  Statistic
#               T = n^{-1/2} sum_i (X_i - mu_x(Z_i)) (Y_i - mu_y(Z_i)),
#             with both conditional means fit in sample.  Calibration draws
#             Xtilde_i ~ Lhat(X | Z_i) and recomputes T without refitting
#             mu_x, exactly as the algorithm prescribes.
#   ndcrt_rf  The \hat{ndCRT} of Niu et al. (2024), Appendix A.  Same
#             resampling scheme, but the statistic is divided by the empirical
#             sd of the residual products, recomputed on each resample.  Its
#             observed value therefore coincides with the GCM statistic.
#
# The estimated conditional law is Lhat(X | Z) = N(mu_x(Z), sigma^2) with a
# single sigma^2 = mean of the squared X-residuals, matching the Gaussian
# model-X working law used in the simulations of Niu et al. (2024).
#
# All four GCM-family procedures share one set of random-forest nuisance fits,
# so they differ only in calibration:
#   asym_gcm_rf   studentized statistic, N(0,1) reference
#   perm_gcm_rf   studentized statistic, marginal-permutation reference (ours)
#   ndcrt_rf      studentized statistic, model-X resampling reference
#   dcrt_rf       unnormalized statistic, model-X resampling reference
#
# Usage:
#   Rscript --vanilla check-dcrt-pcm-additive.R [--reps=N] [--cores=K]
#                                               [--M=M] [--quick]
###############################################################################

suppressPackageStartupMessages({
  library(pcmt)
  library(parallel)
})

## ---------------------------------------------------------------- arguments
cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  hit <- grep(paste0("^--", flag, "="), cli, value = TRUE)
  if (length(hit) == 0L) return(default)
  as.numeric(sub(paste0("^--", flag, "="), "", hit[[1L]]))
}
quick <- "--quick" %in% cli
n_reps <- as.integer(get_opt("reps", if (quick) 40 else 200))
n_cores <- as.integer(get_opt("cores", max(1L, parallel::detectCores() - 1L)))
M_resamp <- as.integer(get_opt("M", 499))
B_perm <- as.integer(get_opt("B", 499))

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
script_dir <- if (length(script_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_argument)))
} else {
  getwd()
}
specifier_file <- file.path(script_dir, "sim-spec-formal", "pcm-additive.R")
source(specifier_file, local = FALSE)   # supplies generate_data(), grid

## ------------------------------------------------------------- dCRT methods
# Mirrors pcmt:::studentized_stat so that the observed ndCRT statistic is
# numerically identical to the observed GCM statistic.
studentized_stat <- function(R) {
  n <- length(R)
  m <- mean(R)
  v <- mean(R^2) - m^2
  if (v <= 0) return(0)
  sqrt(n) * m / sqrt(v)
}

# Two-sided resampling p-value on the magnitude scale, matching
# pcmt:::perm_pvalue(alternative = "two.sided").  Under the Gaussian working
# law the resampling distribution of the unnormalized statistic is exactly
# symmetric about zero, so the magnitude rule is the natural two-sided form.
resample_pvalue <- function(T_obs, T_resamp) {
  M <- length(T_resamp)
  c(two   = (1 + sum(abs(T_resamp) >= abs(T_obs))) / (1 + M),
    right = (1 + sum(T_resamp >= T_obs)) / (1 + M),
    left  = (1 + sum(T_resamp <= T_obs)) / (1 + M))
}

# rX, rY are the observed nuisance residuals X - mu_x(Z) and Y - mu_y(Z).
# Because the resampled statistic only involves Xtilde_i - mu_x(Z_i), and that
# difference is sigma * N(0, 1) under the Gaussian working law, the whole
# resampling loop reduces to draws of standard normals -- mu_x is never refit,
# exactly as Algorithm 1 specifies.
dcrt_pvalues <- function(rX, rY, M) {
  n <- length(rX)
  sigma_hat <- sqrt(mean(rX^2))

  T_dcrt_obs <- sum(rX * rY) / sqrt(n)
  T_ndcrt_obs <- studentized_stat(rX * rY)

  E <- matrix(stats::rnorm(n * M), nrow = n, ncol = M)
  R_resamp <- (sigma_hat * E) * rY            # n x M residual products
  colsum <- colSums(R_resamp)
  m_resamp <- colsum / n
  v_resamp <- colMeans(R_resamp^2) - m_resamp^2
  v_resamp[v_resamp <= 0] <- Inf

  T_dcrt_resamp <- colsum / sqrt(n)
  T_ndcrt_resamp <- sqrt(n) * m_resamp / sqrt(v_resamp)

  d <- resample_pvalue(T_dcrt_obs, T_dcrt_resamp)
  nd <- resample_pvalue(T_ndcrt_obs, T_ndcrt_resamp)
  c(stats::setNames(d, paste0("dcrt_rf_", names(d))),
    stats::setNames(nd, paste0("ndcrt_rf_", names(nd))))
}

## ------------------------------------------------------------ one replicate
# Seeds follow the Section 5.2 convention: one fixed seed per nuisance role, so
# every refit of a role reuses the same learner randomness.
seeds <- c(gcm_permutation = 53112L, gcm_YonZ = 53114L, gcm_XonZ = 53115L,
           dcrt_resampling = 53120L)
rf_args <- function(seed) {
  list(seed = as.integer(seed), num.threads = 1L, num.trees = 500L)
}

one_replicate <- function(rep_id, n, dgp, gamma_signal) {
  set.seed(1000000L + rep_id)
  data <- generate_data(n = n, dgp = dgp, gamma_signal = gamma_signal,
                        B_perm = B_perm)

  gcm_fit <- R.utils::withSeed(
    pcmt::perm_gcm(
      Y = data$Y, X = data$X, Z = data$Z,
      reg_YonZ = "rf", reg_XonZ = "rf",
      args_YonZ = rf_args(seeds[["gcm_YonZ"]]),
      args_XonZ = rf_args(seeds[["gcm_XonZ"]]),
      B = B_perm, n_cores = 1L, h = 15L, alpha = 0.05,
      alternative = "two.sided"
    ),
    seed = seeds[["gcm_permutation"]]
  )

  T_obs <- as.numeric(gcm_fit$statistic)
  dcrt_p <- R.utils::withSeed(
    dcrt_pvalues(c(gcm_fit$rX), c(gcm_fit$rY), M_resamp),
    seed = seeds[["dcrt_resampling"]] + rep_id
  )

  c(perm_gcm_rf_two = as.numeric(gcm_fit$p.value),
    asym_gcm_rf_two   = 2 * stats::pnorm(abs(T_obs), lower.tail = FALSE),
    asym_gcm_rf_right = stats::pnorm(T_obs, lower.tail = FALSE),
    asym_gcm_rf_left  = stats::pnorm(T_obs),
    dcrt_p)
}

## ------------------------------------------------------------- check grid
# A subset of the Section 5.2 grid: the smallest and largest null cell of each
# row, plus the largest signal on that row's power curve.  Permutation GCM
# refits a random forest per draw, so the full 27-cell grid is out of reach for
# a check run; these nine cells cover both ends of each row's sample-size range.
check_cells <- do.call(rbind, lapply(unique(parameter_grid$dgp), function(d) {
  rows <- parameter_grid[parameter_grid$dgp == d, , drop = FALSE]
  nulls <- rows[rows$gamma_signal == 0, , drop = FALSE]
  alts <- rows[rows$gamma_signal > 0, , drop = FALSE]
  alts <- alts[order(alts$gamma_signal), , drop = FALSE]
  keep <- rbind(
    nulls[which.min(nulls$n), , drop = FALSE],
    nulls[which.max(nulls$n), , drop = FALSE],
    alts[nrow(alts), , drop = FALSE]
  )
  keep[, c("n", "dgp", "gamma_signal")]
}))
if (quick) {
  check_cells <- check_cells[check_cells$gamma_signal == 0, , drop = FALSE]
}
rownames(check_cells) <- NULL

cat(sprintf(
  "dCRT check simulation: %d cells x %d reps, M=%d resamples, B=%d perms, %d cores\n",
  nrow(check_cells), n_reps, M_resamp, B_perm, n_cores
))
print(check_cells)
cat("\n")

## ------------------------------------------------------------------- run it
alpha <- 0.05
method_order <- c("asym_gcm_rf", "perm_gcm_rf", "ndcrt_rf", "dcrt_rf")

summary_rows <- list()
pvalue_rows <- list()
for (i in seq_len(nrow(check_cells))) {
  cell <- check_cells[i, ]
  t0 <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(
    seq_len(n_reps),
    function(r) {
      tryCatch(
        one_replicate(r, cell$n, cell$dgp, cell$gamma_signal),
        error = function(e) {
          warning(conditionMessage(e), call. = FALSE)
          stats::setNames(rep(NA_real_, 4L), method_order)
        }
      )
    },
    mc.cores = n_cores
  )
  failed <- !vapply(reps, function(x) is.numeric(x) && length(x) == 4L,
                    logical(1L))
  if (any(failed)) {
    stop(sprintf("%d replicate(s) failed in cell %d.", sum(failed), i))
  }
  P <- do.call(rbind, reps)[, method_order, drop = FALSE]
  elapsed <- proc.time()[["elapsed"]] - t0

  rej <- colMeans(P <= alpha, na.rm = TRUE)
  mcse <- sqrt(rej * (1 - rej) / n_reps)
  summary_rows[[i]] <- data.frame(
    dgp = cell$dgp, n = cell$n, gamma = cell$gamma_signal,
    method = method_order,
    rejection = as.numeric(rej), mcse = as.numeric(mcse),
    stringsAsFactors = FALSE
  )
  pvalue_rows[[i]] <- data.frame(
    dgp = cell$dgp, n = cell$n, gamma = cell$gamma_signal,
    rep = rep(seq_len(n_reps), times = length(method_order)),
    method = rep(method_order, each = n_reps),
    p_value = as.numeric(P),
    stringsAsFactors = FALSE
  )

  cat(sprintf("[%d/%d] %-16s n=%3d gamma=%.4g  (%.0fs)\n",
              i, nrow(check_cells), cell$dgp, cell$n, cell$gamma_signal,
              elapsed))
  cat(paste0(sprintf("        %-12s %.3f (+/- %.3f)", method_order, rej, mcse),
             collapse = "\n"), "\n\n")
  utils::flush.console()
}

results <- do.call(rbind, summary_rows)
pvalues <- do.call(rbind, pvalue_rows)

## -------------------------------------------------------------------- store
out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) {
  file.path(out_root, "simulation-results", "dcrt-check")
} else {
  file.path(dirname(script_dir), "..", "results", "dcrt-check")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(results = results, pvalues = pvalues,
             config = list(n_reps = n_reps, M = M_resamp, B = B_perm,
                           alpha = alpha, cells = check_cells)),
        file.path(out_dir, "dcrt-check-results.rds"))
utils::write.csv(results, file.path(out_dir, "dcrt-check-summary.csv"),
                 row.names = FALSE)

## ------------------------------------------------------------------- report
cat("=== Rejection rate at alpha = 0.05 ===\n")
wide <- stats::reshape(
  results[, c("dgp", "n", "gamma", "method", "rejection")],
  idvar = c("dgp", "n", "gamma"), timevar = "method", direction = "wide"
)
names(wide) <- sub("^rejection\\.", "", names(wide))
print(wide[, c("dgp", "n", "gamma", method_order)], row.names = FALSE,
      digits = 3)
cat("\nWritten to: ", out_dir, "\n", sep = "")
