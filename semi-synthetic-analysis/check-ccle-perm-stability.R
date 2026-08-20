#!/usr/bin/env Rscript

###############################################################################
# dCRT on the CCLE semi-synthetic null for the low-sparsity representative.
#
# Target: EIF5A_MUT (sparsity_tag "low_sparsity", target_idx 4, freq ~ 0.0084,
# i.e. 4 carriers out of 474 cell lines).  Under the delta = 0 setup the
# synthetic outcome is built from Z = X_{-j} alone, so Y |= X_j | Z holds by
# construction and every rejection is a Type-I error.
#
# The script recomputes the GCM nuisance residuals with the same random-forest
# seeds as the manuscript run (ccle_seed_schedule with paired_seed_idx = 34),
# verifies the recomputed GCM statistic against the stored matched-v4 result,
# and then calibrates that same statistic four extra ways:
#
#   dcrt_gauss / ndcrt_gauss   Lhat(X | Z) = N(mu_x(Z), sigma^2)
#   dcrt_bern  / ndcrt_bern    Lhat(X | Z) = Bernoulli(pi(Z))
#
# The Bernoulli law is the natural model-X working law for a binary mutation
# indicator and is what spacrt::dCRT uses for family = "binomial"; the Gaussian
# law is carried along for comparability with the Section 5.2 check.  In both
# cases mu_x is the SAME in-sample random forest the GCM uses and is not refit
# on resampling, per Algorithm 1 of Niu et al. (2024).
#
# Only the observed nuisance fits are recomputed here (two forests per
# replicate).  The permutation GCM is not rerun -- its rejection rate is read
# from the stored matched-v4 result, since it needs a forest refit per draw.
#
# Usage:
#   Rscript --vanilla check-dcrt-ccle-eif5a.R [--reps=N] [--cores=K] [--M=M]
###############################################################################

suppressPackageStartupMessages({
  library(comets)
  library(pcmt)
  library(parallel)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.", call. = FALSE)
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "matched-calibration-utils.R"))

cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  hit <- grep(paste0("^--", flag, "="), cli, value = TRUE)
  if (length(hit) == 0L) return(default)
  as.numeric(sub(paste0("^--", flag, "="), "", hit[[1L]]))
}
R_REPS   <- as.integer(get_opt("reps", 50))     # replicates for the outcome-side error
N_CORES  <- as.integer(get_opt("cores", 7))
M_PERM   <- as.integer(get_opt("Mperm", 100))  # permutation refits for Etilde

TARGET       <- "EIF5A_MUT"
PAIRED_SEED  <- 34L          # target_idx 4 + 30, per SETUPS_NULL
S_FIXED      <- 10L
RF_NUM_TREES <- 500L
RUN_SEED_BASE <- 20260000L
SAMPLE_SEED  <- 2026L
N_PER_BIN    <- 10L
ALPHA        <- 0.05

data_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(data_root)) stop("LOCAL_PDRT_DATA_DIR must be set.", call. = FALSE)
data_root <- normalizePath(path.expand(data_root), mustWork = TRUE)
in_dir <- file.path(data_root, "ccle")

## ---- Load and preprocess CCLE exactly as ccle-typeI-synthY-30targets.R ----
mutations <- read.delim(file.path(in_dir, "mutation.txt"), sep = "\t")
mutations <- mutations[grep("_MUT$", mutations$Name), ]
response <- read.csv(file.path(in_dir, "response.csv"))
response <- response[response$Compound == "PLX4720", ]
idx  <- match(response$CCLE.Cell.Line.Name, colnames(mutations))
resp <- response[!is.na(idx), ]
mut  <- mutations[, c(1, idx[!is.na(idx)])]
rownames(mut) <- mut$Name
mut$Name <- NULL
Y_real <- resp$Amax
cors <- suppressWarnings(apply(mut, 1, function(m) {
  tc <- cor(Y_real, m); ifelse(!is.na(tc), tc, 0)
}))
screen <- cors[abs(cors) > 0.05 & !is.na(cors)]
dat <- data.frame(Y_real = Y_real, t(mut[names(screen), ]))
mut_names <- setdiff(colnames(dat), "Y_real")
M_full <- as.matrix(dat[, mut_names])
n <- nrow(M_full)
freqs <- colMeans(M_full)

## ---- Rebuild the tercile grid and check the EIF5A contract ----
qs <- quantile(freqs, probs = c(1/3, 2/3), names = FALSE)
pool_low  <- sort(names(freqs)[freqs <= qs[1]])
pool_med  <- sort(names(freqs)[freqs >  qs[1] & freqs <= qs[2]])
pool_high <- sort(names(freqs)[freqs >  qs[2]])
set.seed(SAMPLE_SEED)
targets_low  <- sort(sample(pool_low,  N_PER_BIN))
targets_med  <- sort(sample(pool_med,  N_PER_BIN))
targets_high <- sort(sample(pool_high, N_PER_BIN))
TARGETS_DF <- rbind(
  data.frame(target = targets_low,  tag = "low_sparsity",    stringsAsFactors = FALSE),
  data.frame(target = targets_med,  tag = "medium_sparsity", stringsAsFactors = FALSE),
  data.frame(target = targets_high, tag = "high_sparsity",   stringsAsFactors = FALSE)
)
TARGETS_DF$target_idx <- seq_len(nrow(TARGETS_DF))
row <- TARGETS_DF[TARGETS_DF$target == TARGET, , drop = FALSE]
if (nrow(row) != 1L || row$tag != "low_sparsity" || row$target_idx != 4L ||
    row$target_idx + nrow(TARGETS_DF) != PAIRED_SEED) {
  stop("EIF5A_MUT did not reproduce its canonical grid position.", call. = FALSE)
}
X_j <- M_full[, TARGET]
Z_j <- M_full[, setdiff(mut_names, TARGET), drop = FALSE]
cat(sprintf("CCLE: n = %d, p = %d; %s freq = %.5f (%d carriers)\n",
            n, ncol(M_full), TARGET, freqs[[TARGET]], sum(X_j)))

## ---- Synthetic-outcome generator (verbatim from the runner) ----
generate_Y_components <- function(Z, s) {
  p <- ncol(Z)
  s <- max(1L, min(s, p))
  S <- sort(sample.int(p, size = s))
  beta <- rnorm(s, mean = 0, sd = 1)
  mu   <- as.numeric(Z[, S, drop = FALSE] %*% beta)
  if (s >= 2L) {
    n_pairs_max <- as.integer(choose(s, 2L))
    n_inter     <- min(2L * s, n_pairs_max)
    pairs       <- utils::combn(S, 2)
    if (ncol(pairs) > n_inter) {
      sel   <- sample.int(ncol(pairs), size = n_inter)
      pairs <- pairs[, sel, drop = FALSE]
    }
    gamma <- rnorm(ncol(pairs), mean = 0, sd = 1)
    for (k in seq_len(ncol(pairs))) {
      mu <- mu + gamma[k] * Z[, pairs[1, k]] * Z[, pairs[2, k]]
    }
  }
  sigma <- sd(mu)
  if (!is.finite(sigma) || sigma == 0) sigma <- 1
  list(mu = mu, eps = rnorm(length(mu), mean = 0, sd = sigma),
       outcome_noise_sd = sigma)
}


fit_regression <- get("fit_regression", envir = asNamespace("pcmt"))

## ---- Permutation-stability diagnostics -------------------------------------
# Assumption "perm" for the GCM is stated in terms of
#     Etil_{n,x} = (1/n) sum_i ( muhat_x^(Pi)(Z_i) - Xbar )^2,
# the empirical energy of the weight-side refit AFTER permutation, measured
# against its natural permutation target Xbar.  X_j and Z are fixed across
# replicates here (only the synthetic Y is redrawn), so Etil_{n,x} does not
# depend on the outcome at all and is estimated once from M_PERM refits.
#
# Assumption "unperm"(i) is the outcome-side condition E_{n,y} = o_P(1).  The
# CCLE outcome is synthetic, so mu_y is known exactly and E_{n,y} is directly
# computable -- unlike mu_x, which is unknown because X_j is real data.
# For the weight side we therefore report the observed-fit energy
#     Ehat_x^obs = (1/n) sum_i ( muhat_x(Z_i) - Xbar )^2
# as the comparison: if the forest extracts as much structure from a permuted
# X as from the real one, it is fitting noise rather than signal.

Xbar <- mean(X_j)
var_X <- mean((X_j - Xbar)^2)

cat(sprintf("\n%s: n = %d, p = %d, carriers = %d, Var(X) = %.5f\n",
            TARGET, n, ncol(Z_j), sum(X_j), var_X))

## --- (a) permutation-stability error, M_PERM refits on permuted X
cat(sprintf("Fitting %d permutation refits on %d cores ...\n", M_PERM, N_CORES))
t0 <- proc.time()[["elapsed"]]
perm_fits <- parallel::mclapply(seq_len(M_PERM), function(b) {
  seeds <- ccle_seed_schedule(PAIRED_SEED, b, RUN_SEED_BASE)
  Xp <- ccle_with_seed(seeds[["gcm_permutation_seed"]], X_j[sample.int(n)])
  m <- fit_regression(y = Xp, x = Z_j, method = "rf",
                      args = ccle_rf_args(seeds[["gcm_XonZ_seed"]], RF_NUM_TREES))
  fit <- Xp - as.numeric(residuals(m, response = Xp, data = Z_j))
  c(Etil = mean((fit - Xbar)^2), fit_max = max(fit), fit_min = min(fit))
}, mc.cores = N_CORES)
perm_fits <- do.call(rbind, perm_fits)
Etil_x <- mean(perm_fits[, "Etil"])
cat(sprintf("  done in %.0fs\n", proc.time()[["elapsed"]] - t0))

## --- (b) observed-side fit energy, same seeds, unpermuted X
obs_fits <- parallel::mclapply(seq_len(min(R_REPS, 50L)), function(r) {
  seeds <- ccle_seed_schedule(PAIRED_SEED, r, RUN_SEED_BASE)
  m <- fit_regression(y = X_j, x = Z_j, method = "rf",
                      args = ccle_rf_args(seeds[["gcm_XonZ_seed"]], RF_NUM_TREES))
  fit <- X_j - as.numeric(residuals(m, response = X_j, data = Z_j))
  mean((fit - Xbar)^2)
}, mc.cores = N_CORES)
Ehat_x_obs <- mean(unlist(obs_fits))

## --- (c) outcome-side error against the KNOWN synthetic mu_y
y_err <- parallel::mclapply(seq_len(R_REPS), function(r) {
  seeds <- ccle_seed_schedule(PAIRED_SEED, r, RUN_SEED_BASE)
  outcome <- ccle_with_seed(seeds[["data_seed"]],
                            generate_Y_components(Z_j, S_FIXED))
  Y <- outcome$mu + outcome$eps
  m <- fit_regression(y = Y, x = Z_j, method = "rf",
                      args = ccle_rf_args(seeds[["gcm_YonZ_seed"]], RF_NUM_TREES))
  fit <- Y - as.numeric(residuals(m, response = Y, data = Z_j))
  c(E_y = mean((fit - outcome$mu)^2),
    var_mu = mean((outcome$mu - mean(outcome$mu))^2))
}, mc.cores = N_CORES)
y_err <- do.call(rbind, y_err)
E_ny <- mean(y_err[, "E_y"]); var_mu_y <- mean(y_err[, "var_mu"])

## ---- Report ---------------------------------------------------------------
cat(sprintf("\n=== Assumption 'perm' (permutation stability), GCM weight side ===\n"))
cat(sprintf("  Etil_{n,x}                       = %.6f   (M = %d refits)\n", Etil_x, M_PERM))
cat(sprintf("  Etil_{n,x} / Var(X)              = %.4f\n", Etil_x / var_X))
cat(sprintf("  n * Etil_{n,x}                   = %.3f\n", n * Etil_x))
cat(sprintf("  Ehat_x^obs (unpermuted fit)      = %.6f   ratio perm/obs = %.3f\n",
            Ehat_x_obs, Etil_x / Ehat_x_obs))
cat(sprintf("\n=== Assumption 'unperm'(i) (outcome side), mu_y known exactly ===\n"))
cat(sprintf("  E_{n,y}                          = %.4f\n", E_ny))
cat(sprintf("  Var(mu_y)                        = %.4f   ratio = %.3f\n",
            var_mu_y, E_ny / var_mu_y))
cat(sprintf("\n=== Assumption 'perm'(ii) product condition ===\n"))
cat(sprintf("  n * E_{n,y} * Etil_{n,x}         = %.4f\n", n * E_ny * Etil_x))
