#!/usr/bin/env Rscript
###############################################################################
# Does raising the covariate dimension repair the Section 5.2 random forest?
#
# Section 5.2 uses a one-dimensional Z.  With a single covariate ranger must
# take mtry = 1, so there is no feature subsampling, the forest degenerates to
# bagged unpruned trees, and its nuisance error does not vanish: E_x is flat at
# roughly 0.34 from n = 400 to n = 6400 against Var(mu_x) ~ 0.50.  Assumption
# "unperm"(ii) then fails, the asymptotic theory does not cover the simulation,
# and its asymptotic tests sit near 0.10 at every n.
#
# This script holds the DGP fixed and varies ONLY the covariate dimension: the
# signal is sin(2 pi Z1) exactly as in Section 5.2, with d - 1 further standard
# normal covariates that enter nothing.  At d = 3 ranger can take mtry < p, so
# the forest is no longer degenerate.  Comparing d = 1 with d = 3 therefore
# isolates the forest's behaviour from the difficulty of the target.
#
# ONLY the asymptotic tests are computed here.  The question is whether they
# recover once the nuisance is no longer degenerate, and that needs no
# permutation reference -- dropping it removes every per-permutation forest
# refit and makes the run roughly a hundred times cheaper, which buys tighter
# Monte Carlo error and a range of n instead of a single point.
#
# Each cell reports the nuisance errors and the null spread of the GCM
# statistic: an sd far from 1 is precisely what the miscalibration is.
#
# Usage: Rscript explore-rf-dimension.R [--reps=2000] [--cores=7]
#          [--n=100,200,400] [--d=1,3]
###############################################################################
suppressPackageStartupMessages({library(pcmt); library(parallel)})

cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(f, d) { h <- grep(paste0("^--", f, "="), cli, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", f, "="), "", h[[1L]]) }
R_REPS  <- as.integer(get_opt("reps", 2000))
N_CORES <- as.integer(get_opt("cores", 7))
N_GRID  <- as.integer(strsplit(get_opt("n", "100,200,400"), ",")[[1L]])
D_GRID  <- as.integer(strsplit(get_opt("d", "1,3"), ",")[[1L]])
ALPHA   <- 0.05
RF <- function(seed) list(seed = as.integer(seed), num.threads = 1L,
                          num.trees = 500L)

generate <- function(n, dgp, d) {
  Z <- matrix(rnorm(n * d), n, d); colnames(Z) <- paste0("Z", seq_len(d))
  s <- sin(2 * pi * Z[, 1L]); eps <- rnorm(n)
  if (dgp == "setting1") { X <- s + rnorm(n); Y <- s + eps
  } else if (dgp == "setting2") { xi <- rexp(n, 1) - 1
    X <- s - s * xi; Y <- s + eps
  } else { X <- rnorm(n); Y <- s + eps }          # setting1_xindep
  list(X = X, Y = Y, Z = Z,
       mu_x = if (dgp == "setting1_xindep") rep(0, n) else s, mu_y = s)
}

DGPS <- c("setting1_xindep", "setting1", "setting2")

fit_regression <- get("fit_regression", envir = asNamespace("pcmt"))
studentized_stat <- function(v) {
  m <- mean(v); s2 <- mean(v^2) - m^2
  if (s2 <= 0) return(0)
  sqrt(length(v)) * m / sqrt(s2)
}

one_rep <- function(r, n, dgp, d) {
  set.seed(6000000L + 1013L * match(dgp, DGPS) + 37L * d + 7L * n + r)
  dd <- generate(n, dgp, d)
  # --- asymptotic GCM: two nuisance fits, one studentized statistic
  mY <- fit_regression(y = dd$Y, x = dd$Z, method = "rf", args = RF(53114L))
  rY <- as.numeric(stats::residuals(mY, response = dd$Y, data = dd$Z))
  mX <- fit_regression(y = dd$X, x = dd$Z, method = "rf", args = RF(53115L))
  rX <- as.numeric(stats::residuals(mX, response = dd$X, data = dd$Z))
  Tg <- studentized_stat(rX * rY)
  # --- asymptotic PCM.  pcmt::asym_pcm is restricted to lrm/gols nuisances,
  # so the observed statistic is taken from perm_pcm with B = 1 (the
  # permutation reference is discarded) and the normal reference applied to it.
  # This is exactly how the Section 5.2 specifier forms asym_pcm, so the two
  # are directly comparable.
  ap <- pcmt::perm_pcm(Y = dd$Y, X = dd$X, Z = dd$Z, reg_YonXZ = "rf",
        reg_YonZ = "rf", reg_YhatonZ = "rf", reg_RonZ = "rf",
        args_YonXZ = RF(53116L), args_YonZ = RF(53117L),
        args_YhatonZ = RF(53118L), args_RonZ = RF(53119L), frac = 0.5,
        B = 1L, n_cores = 1L, h = NULL, alternative = "greater")
  c(asym_gcm = 2 * stats::pnorm(abs(Tg), lower.tail = FALSE),
    asym_pcm = stats::pnorm(as.numeric(ap$statistic), lower.tail = FALSE),
    stat_gcm = Tg, stat_pcm = as.numeric(ap$statistic),
    E_x = mean((dd$X - rX - dd$mu_x)^2),
    E_y = mean((dd$Y - rY - dd$mu_y)^2),
    var_mu_x = mean((dd$mu_x - mean(dd$mu_x))^2))
}

wilson <- function(k, N) { z <- 1.959964; ph <- k/N; dn <- 1 + z^2/N
  ctr <- (ph + z^2/(2*N))/dn; hw <- (z/dn)*sqrt(ph*(1-ph)/N + z^2/(4*N^2))
  c(max(0, ctr-hw), min(1, ctr+hw)) }

cat(sprintf("RF dimension probe, ASYMPTOTIC TESTS ONLY: null cells, R=%d\n", R_REPS))
cat("Signal is sin(2*pi*Z1) throughout; d-1 covariates are inert.\n\n")
cat(sprintf("%-15s %2s %4s | %-18s %-18s | %6s %6s %6s %6s\n",
            "design","d","n","asym_gcm","asym_pcm",
            "E_x","Var_mx","sdTg","sdTp"))
rows <- list()
for (dgp in DGPS) for (d in D_GRID) for (n in N_GRID) {
  M <- do.call(rbind, mclapply(seq_len(R_REPS),
        function(r) one_rep(r, n, dgp, d), mc.cores = N_CORES))
  f <- function(nm) { k <- sum(M[, nm] <= ALPHA); ci <- wilson(k, R_REPS)
    sprintf("%.3f[%.3f,%.3f]", k/R_REPS, ci[1], ci[2]) }
  cat(sprintf("%-15s %2d %4d | %-18s %-18s | %6.3f %6.3f %6.2f %6.2f\n",
      dgp, d, n, f("asym_gcm"), f("asym_pcm"),
      mean(M[,"E_x"]), mean(M[,"var_mu_x"]),
      stats::sd(M[,"stat_gcm"]), stats::sd(M[,"stat_pcm"])))
  utils::flush.console()
  rows[[length(rows)+1L]] <- data.frame(dgp=dgp, d=d, n=n, R=R_REPS,
    asym_gcm=mean(M[,"asym_gcm"]<=ALPHA), asym_pcm=mean(M[,"asym_pcm"]<=ALPHA),
    E_x=mean(M[,"E_x"]), E_y=mean(M[,"E_y"]), var_mu_x=mean(M[,"var_mu_x"]),
    sd_stat_gcm=stats::sd(M[,"stat_gcm"]), sd_stat_pcm=stats::sd(M[,"stat_pcm"]))
}
res <- do.call(rbind, rows)
out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) file.path(out_root,"simulation-results","rf-dimension") else file.path(getwd(),"results","rf-dimension")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out_dir, "rf-dimension.rds"))
cat("\nWritten to: ", out_dir, "\n", sep = "")
