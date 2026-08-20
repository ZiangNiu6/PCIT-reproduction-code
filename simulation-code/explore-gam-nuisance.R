#!/usr/bin/env Rscript
###############################################################################
# Section 5.2 exactly as specified (d = 1), with a spline nuisance instead of a
# random forest.
#
# Motivation.  The published Section 5.2 uses ranger on a one-dimensional Z.
# With a single covariate mtry is forced to 1, the forest degenerates to bagged
# unpruned trees, and its error does not vanish (E_x flat at ~0.34 from n = 400
# to 6400 against Var(mu_x) ~ 0.50).  Assumption "unperm"(ii) then fails, so the
# asymptotic theory does not cover the simulation -- which is why its asymptotic
# tests sit near 0.10 at every n.  A penalised spline is the natural learner for
# a smooth one-dimensional signal, and is far more accurate here.
#
# Basis dimension matters more than the choice of learner.  mgcv's default of
# 10 basis functions cannot represent sin(2 pi Z1): it gives E_x ~ 0.30-0.41,
# no better than the forest and no improvement with n.  Lundborg et al. use
# floor((N-1)/d) basis functions, which works but costs 27 s per fit at
# n = 600 -- impossible inside a permutation loop.  k = 20 matches that
# accuracy (E_x 0.141, 0.090, 0.037 at n = 100, 200, 600) at 22 ms.
#
# The question: once the nuisance is competent, do the asymptotic tests recover,
# and does permutation calibration retain an edge?
#
# Usage: Rscript explore-gam-nuisance.R [--reps=300] [--cores=7] [--k=20]
###############################################################################
suppressPackageStartupMessages({library(pcmt); library(parallel)})

cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(f, d) { h <- grep(paste0("^--", f, "="), cli, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", f, "="), "", h[[1L]]) }
R_REPS  <- as.integer(get_opt("reps", 300))
N_CORES <- as.integer(get_opt("cores", 7))
B_PERM  <- as.integer(get_opt("B", 499))
GAM_K   <- as.integer(get_opt("k", 20))
ALPHA   <- 0.05

RF  <- function(seed) list(seed = as.integer(seed), num.threads = 1L,
                           num.trees = 500L)
# mgcv requires the basis dimension to be below the number of unique covariate
# values, and the PCM halves the sample again, so k is capped at n/4: that
# leaves a factor-of-two margin inside the PCM's inference fold.
GAM <- function(seed, n) list(family = stats::gaussian(),
                              k = max(3L, min(GAM_K, floor(n / 4))))

# Section 5.2's designs verbatim: Z is one-dimensional throughout.
generate <- function(n, dgp) {
  Z1 <- rnorm(n); Z <- matrix(Z1, ncol = 1L); colnames(Z) <- "Z1"
  s <- sin(2 * pi * Z1); eps <- rnorm(n)
  if (dgp == "setting1") { X <- s + rnorm(n); Y <- s + eps
  } else if (dgp == "setting2") { xi <- rexp(n, 1) - 1
    X <- s - s * xi; Y <- s + eps
  } else { X <- rnorm(n); Y <- s + eps }
  list(X = X, Y = Y, Z = Z,
       mu_x = if (dgp == "setting1_xindep") rep(0, n) else s, mu_y = s)
}

# --designs and --n override the default grid.  Section 5.2 runs setting1 and
# setting2 only at n = 200-600; the permutation advantage, where it survives a
# competent nuisance at all, appears well below that, so those rows are worth
# probing at sample sizes the published grid never visits.
designs_opt <- get_opt("designs", "")
n_opt <- get_opt("n", "")
CELLS <- if (nzchar(designs_opt) && nzchar(n_opt)) {
  expand.grid(dgp = strsplit(designs_opt, ",")[[1L]],
              n = as.integer(strsplit(n_opt, ",")[[1L]]),
              stringsAsFactors = FALSE)[, c("dgp", "n")]
} else {
  rbind(
    expand.grid(dgp = "setting1_xindep", n = c(20L, 60L, 100L),
                stringsAsFactors = FALSE),
    expand.grid(dgp = "setting1", n = c(200L, 400L, 600L),
                stringsAsFactors = FALSE),
    expand.grid(dgp = "setting2", n = c(200L, 400L, 600L),
                stringsAsFactors = FALSE))
}
CELLS <- CELLS[order(CELLS$dgp, CELLS$n), ]
DGPS <- c("setting1_xindep", "setting1", "setting2")

one_rep <- function(r, n, dgp, learner) {
  set.seed(7000000L + 1013L * match(dgp, DGPS) + 7L * n +
             97L * match(learner, c("rf", "gam")) + r)
  dd <- generate(n, dgp)
  A <- if (learner == "rf") function(sd) RF(sd) else function(sd) GAM(sd, n)
  g <- pcmt::perm_gcm(Y = dd$Y, X = dd$X, Z = dd$Z,
        reg_YonZ = learner, reg_XonZ = learner,
        args_YonZ = A(53114L), args_XonZ = A(53115L),
        B = B_PERM, n_cores = 1L, h = 15L, alpha = ALPHA,
        alternative = "two.sided")
  p <- pcmt::perm_pcm(Y = dd$Y, X = dd$X, Z = dd$Z,
        reg_YonXZ = learner, reg_YonZ = learner,
        reg_YhatonZ = learner, reg_RonZ = learner,
        args_YonXZ = A(53116L), args_YonZ = A(53117L),
        args_YhatonZ = A(53118L), args_RonZ = A(53119L), frac = 0.5,
        B = B_PERM, n_cores = 1L, h = 15L, alpha = ALPHA,
        alternative = "greater")
  c(perm_gcm = as.numeric(g$p.value),
    asym_gcm = 2 * stats::pnorm(abs(as.numeric(g$statistic)), lower.tail = FALSE),
    perm_pcm = as.numeric(p$p.value),
    asym_pcm = stats::pnorm(as.numeric(p$statistic), lower.tail = FALSE),
    stat_gcm = as.numeric(g$statistic),
    E_x = mean((c(dd$X) - c(g$rX) - dd$mu_x)^2),
    E_y = mean((c(dd$Y) - c(g$rY) - dd$mu_y)^2))
}

wilson <- function(k, N) { z <- 1.959964; ph <- k/N; dn <- 1 + z^2/N
  ctr <- (ph + z^2/(2*N))/dn; hw <- (z/dn)*sqrt(ph*(1-ph)/N + z^2/(4*N^2))
  c(max(0, ctr-hw), min(1, ctr+hw)) }

cat(sprintf("Section 5.2 designs (d=1), null cells.  R=%d, B=%d, GAM k=min(%d, n/4)\n\n",
            R_REPS, B_PERM, GAM_K))
cat(sprintf("%-15s %4s %-4s | %-15s %-15s | %-15s %-15s | %6s %6s\n",
            "design","n","reg","perm_gcm","asym_gcm","perm_pcm","asym_pcm",
            "E_x","sd(T)"))
rows <- list()
for (i in seq_len(nrow(CELLS))) for (learner in c("rf", "gam")) {
  dgp <- CELLS$dgp[i]; n <- CELLS$n[i]
  M <- do.call(rbind, mclapply(seq_len(R_REPS),
        function(r) one_rep(r, n, dgp, learner), mc.cores = N_CORES))
  f <- function(nm) { k <- sum(M[, nm] <= ALPHA); ci <- wilson(k, R_REPS)
    sprintf("%.3f[%.3f,%.3f]", k/R_REPS, ci[1], ci[2]) }
  cat(sprintf("%-15s %4d %-4s | %-15s %-15s | %-15s %-15s | %6.3f %6.2f\n",
      dgp, n, learner, f("perm_gcm"), f("asym_gcm"), f("perm_pcm"),
      f("asym_pcm"), mean(M[,"E_x"]), stats::sd(M[,"stat_gcm"])))
  utils::flush.console()
  rows[[length(rows)+1L]] <- data.frame(dgp=dgp, n=n, reg=learner, R=R_REPS,
    perm_gcm=mean(M[,"perm_gcm"]<=ALPHA), asym_gcm=mean(M[,"asym_gcm"]<=ALPHA),
    perm_pcm=mean(M[,"perm_pcm"]<=ALPHA), asym_pcm=mean(M[,"asym_pcm"]<=ALPHA),
    E_x=mean(M[,"E_x"]), E_y=mean(M[,"E_y"]),
    sd_stat=stats::sd(M[,"stat_gcm"]))
}
res <- do.call(rbind, rows)
out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) file.path(out_root,"simulation-results","gam-nuisance") else file.path(getwd(),"results","gam-nuisance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(out_dir, "gam-nuisance.rds"))
cat("\nWritten to: ", out_dir, "\n", sep = "")
