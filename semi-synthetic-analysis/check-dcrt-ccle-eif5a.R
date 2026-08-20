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
R_REPS   <- as.integer(get_opt("reps", 500))
N_CORES  <- as.integer(get_opt("cores", 4))
M_RESAMP <- as.integer(get_opt("M", 499))

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

## ---- dCRT machinery ----
fit_regression <- get("fit_regression", envir = asNamespace("pcmt"))
studentized_stat <- function(Rv) {
  nn <- length(Rv); m <- mean(Rv); v <- mean(Rv^2) - m^2
  if (v <= 0) return(0)
  sqrt(nn) * m / sqrt(v)
}
# Sided resampling p-values, following pcmt:::perm_pvalue exactly:
#   two.sided = (1 + #{|T~| >= |T_obs|}) / (1 + M)
#   greater   = (1 + #{ T~  >=  T_obs }) / (1 + M)
#   less      = (1 + #{ T~  <=  T_obs }) / (1 + M)
resample_pvalue <- function(T_obs, T_res) {
  M <- length(T_res)
  c(two   = (1 + sum(abs(T_res) >= abs(T_obs))) / (1 + M),
    right = (1 + sum(T_res >= T_obs)) / (1 + M),
    left  = (1 + sum(T_res <= T_obs)) / (1 + M))
}
# Xtilde is drawn from Lhat(X | Z); mu_x is NOT refit (Algorithm 1).
dcrt_from_draws <- function(resid_mat, rY, T_dcrt_obs, T_stud_obs) {
  Rm <- resid_mat * rY
  cs <- colSums(Rm)
  nn <- nrow(Rm)
  m <- cs / nn
  v <- colMeans(Rm^2) - m^2
  v[v <= 0] <- Inf
  d <- resample_pvalue(T_dcrt_obs, cs / sqrt(nn))
  nd <- resample_pvalue(T_stud_obs, sqrt(nn) * m / sqrt(v))
  stats::setNames(c(d, nd), c(paste0("dcrt_", names(d)), paste0("ndcrt_", names(nd))))
}

one_rep <- function(r) {
  seeds <- ccle_seed_schedule(PAIRED_SEED, r, RUN_SEED_BASE)
  outcome <- ccle_with_seed(seeds[["data_seed"]],
                            generate_Y_components(Z_j, S_FIXED))
  Y <- outcome$mu + outcome$eps          # effect_size = 0: the null

  mY <- fit_regression(y = Y, x = Z_j, method = "rf",
                       args = ccle_rf_args(seeds[["gcm_YonZ_seed"]], RF_NUM_TREES))
  rY <- as.numeric(residuals(mY, response = Y, data = Z_j))
  mX <- fit_regression(y = X_j, x = Z_j, method = "rf",
                       args = ccle_rf_args(seeds[["gcm_XonZ_seed"]], RF_NUM_TREES))
  rX <- as.numeric(residuals(mX, response = X_j, data = Z_j))

  T_stud <- studentized_stat(rX * rY)
  T_dcrt <- sum(rX * rY) / sqrt(n)
  pi_hat <- pmin(pmax(X_j - rX, 0), 1)       # fitted P(X = 1 | Z), clipped
  sigma_hat <- sqrt(mean(rX^2))

  set.seed(seeds[["gcm_permutation_seed"]])
  g <- dcrt_from_draws(matrix(rnorm(n * M_RESAMP, sd = sigma_hat), n, M_RESAMP),
                       rY, T_dcrt, T_stud)
  b <- dcrt_from_draws(matrix(rbinom(n * M_RESAMP, 1L, pi_hat), n, M_RESAMP) - pi_hat,
                       rY, T_dcrt, T_stud)

  c(asym_gcm_two   = 2 * stats::pnorm(abs(T_stud), lower.tail = FALSE),
    asym_gcm_right = stats::pnorm(T_stud, lower.tail = FALSE),
    asym_gcm_left  = stats::pnorm(T_stud),
    stats::setNames(g, paste0(sub("_", "_gauss_", names(g)))),
    stats::setNames(b, paste0(sub("_", "_bern_", names(b)))),
    T_stud = T_stud, sigma_hat = sigma_hat, pi_hat_max = max(pi_hat),
    n_carriers_resampled = mean(colSums(matrix(rbinom(n * 20L, 1L, pi_hat), n, 20L))))
}

cat(sprintf("Running %d replicates on %d cores, M = %d resamples ...\n",
            R_REPS, N_CORES, M_RESAMP))
t0 <- proc.time()[["elapsed"]]
reps <- parallel::mclapply(seq_len(R_REPS), one_rep, mc.cores = N_CORES)
bad <- !vapply(reps, is.numeric, logical(1L))
if (any(bad)) stop(sprintf("%d replicate(s) failed.", sum(bad)))
Rm <- do.call(rbind, reps)
cat(sprintf("Done in %.0fs\n\n", proc.time()[["elapsed"]] - t0))

## ---- Cross-check against the stored matched-v4 result ----
stored_path <- file.path(in_dir, sprintf(
  "ccle-typeI-synthY-30targets-matched-v4-results-%s-delta0-s10.rds", TARGET))
stored <- readRDS(stored_path)[[1L]]
# The stored matched-v4 result was produced on the Linux cluster.  ranger seeds
# its mt19937_64 through std::uniform_int_distribution, whose implementation
# differs between libstdc++ and libc++, so a forest grown here from the same
# seed is not bitwise identical to the stored one.  The synthetic outcome and
# the seed schedule DO reproduce exactly (verified separately), so the stored
# run is a valid distributional reference but not a bitwise one.  What matters
# for this comparison is internal matching: every method below is calibrated
# from one and the same pair of nuisance fits.
k <- min(R_REPS, nrow(stored$statistics))
cmp <- data.frame(mine = Rm[seq_len(k), "T_stud"],
                  stored = stored$statistics[seq_len(k), "gcm"])
cat(sprintf("GCM statistic, recomputed vs stored (%d reps):\n", k))
cat(sprintf("  sd:   mine %.4f   stored %.4f\n", sd(cmp$mine), sd(cmp$stored)))
cat(sprintf("  mean: mine %.4f   stored %.4f\n", mean(cmp$mine), mean(cmp$stored)))
cat(sprintf("  max |deviation| %.3g  (nonzero: cross-platform ranger, see note)\n",
            max(abs(cmp$mine - cmp$stored))))

# asym_pcm needs no rerun: PCM's asymptotic p-value is a deterministic function
# of the stored statistic, so all three sidedness variants are read off directly.
T_pcm <- stored$statistics[, "pcm"]
asym_pcm <- c(two   = mean(2 * stats::pnorm(abs(T_pcm), lower.tail = FALSE) <= ALPHA),
              right = mean(stats::pnorm(T_pcm, lower.tail = FALSE) <= ALPHA),
              left  = mean(stats::pnorm(T_pcm) <= ALPHA))

methods <- c("asym_gcm", "dcrt_gauss", "ndcrt_gauss", "dcrt_bern", "ndcrt_bern")
sides <- c("two", "right", "left")
rej <- vapply(methods, function(m)
  colMeans(Rm[, paste0(m, "_", sides), drop = FALSE] <= ALPHA), numeric(3L))
rownames(rej) <- sides
rej <- cbind(rej, asym_pcm = asym_pcm[sides])
mcse <- sqrt(rej * (1 - rej) / R_REPS)

cat(sprintf("\n=== %s, delta = 0 (CI null), rejection at alpha = %.2f ===\n",
            TARGET, ALPHA))
cat(sprintf("%-14s %9s %9s %9s\n", "method", "two.sided", "greater", "less"))
for (m in colnames(rej)) {
  cat(sprintf("%-14s %9s %9s %9s\n", m,
              sprintf("%.3f", rej["two", m]),
              sprintf("%.3f", rej["right", m]),
              sprintf("%.3f", rej["left", m])))
}
cat(sprintf("(MCSE at 0.05 is %.3f on %d reps; asym_pcm is from the stored run)\n",
            sqrt(0.05 * 0.95 / R_REPS), R_REPS))
cat(sprintf("\nPCM statistic under the null: mean %.3f, sd %.3f\n",
            mean(T_pcm), sd(T_pcm)))
cat(sprintf("GCM statistic under the null: mean %.3f, sd %.3f\n",
            mean(Rm[, "T_stud"]), sd(Rm[, "T_stud"])))
cat(sprintf("\n  %-14s %.3f   [stored matched-v4, two.sided, %d reps]\n", "perm_gcm",
            stored$rejection_rate[["perm_gcm"]], stored$R))
cat(sprintf("  %-14s %.3f   [stored matched-v4, two.sided]\n", "perm_pcm",
            stored$rejection_rate[["perm_pcm"]]))
cat(sprintf("\nDiagnostics: mean sigma_hat = %.4f, mean max pi_hat = %.3f,\n",
            mean(Rm[, "sigma_hat"]), mean(Rm[, "pi_hat_max"])))
cat(sprintf("  observed carriers = %d, mean carriers per Bernoulli resample = %.2f\n",
            sum(X_j), mean(Rm[, "n_carriers_resampled"])))

out_dir <- file.path(data_root, "simulation-results", "dcrt-check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(raw = Rm, rejection = rej, asym_pcm = asym_pcm, target = TARGET,
             stored_rejection = stored$rejection_rate,
             config = list(R = R_REPS, M = M_RESAMP, alpha = ALPHA)),
        file.path(out_dir, "dcrt-ccle-eif5a-null.rds"))
cat("\nWritten to: ", out_dir, "\n", sep = "")
