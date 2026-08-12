#!/usr/bin/env Rscript
###############################################################################
# CCLE Semi-Synthetic Type I Error -- 30 targets x s in {5, 10}
#
# Companion to ccle-typeI-synthY.R.  Whereas that script ran 3 hand-picked
# targets at a single sparsity level (s = 5), this one sweeps a broader grid:
#
#   * 10 targets sampled from each frequency tercile of the screened mutations
#     (low / medium / high sparsity)              ->  30 targets
#   * active-set size s in {5, 10}                ->   2 s values
#   * 30 x 2 = 60 setups, each with R = 500 reps
#
# The runner takes a single CLI integer arg in 1..60, the setup index into the
# (target, s) grid (see SETUPS below), and writes one RDS file per setup:
#
#   $LOCAL_PDRT_DATA_DIR/ccle/
#     ccle-typeI-synthY-30targets-results-{TARGET}-s{S}.rds
#
# After all 60 jobs finish, run combine-ccle-results.R to merge.
###############################################################################

suppressPackageStartupMessages({
  library(comets)
  library(pcmt)
})

data_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(data_root)) {
  stop(
    "LOCAL_PDRT_DATA_DIR must point to the external data/results directory.",
    call. = FALSE
  )
}
data_root <- path.expand(data_root)
if (!dir.exists(data_root)) {
  stop("LOCAL_PDRT_DATA_DIR does not exist: ", data_root, call. = FALSE)
}
data_root <- normalizePath(data_root, mustWork = TRUE)

# ---- Config ----
R         <- 500L
B_PERM    <- 499L
H_THRESH  <- 15L
PCM_REPS  <- 10L
ALPHA     <- 0.05
N_CORES   <- 1L
S_VALUES  <- c(5L, 10L)
N_PER_BIN <- 10L
SAMPLE_SEED <- 2026L      # seed for tercile sampling (must be stable across jobs)
RUN_SEED_BASE <- 20260000L  # offset added to setup index for the per-rep seed

# ---- 1. Load + preprocess CCLE ----
# Raw inputs and generated results share this external CCLE directory. Keeping
# both outside the repository prevents large replicate files from becoming
# source-controlled artifacts.
in_dir <- file.path(data_root, "ccle")
bpath <- in_dir
input_files <- file.path(in_dir, c("mutation.txt", "response.csv"))
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing CCLE input file(s):\n  ",
    paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}

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
M_full    <- as.matrix(dat[, mut_names])
n         <- nrow(M_full)
p_total   <- ncol(M_full)
freqs     <- colMeans(M_full)

cat(sprintf("Loaded CCLE: n = %d cell lines, p = %d screened mutations\n",
            n, p_total))

# ---- 2. Build the 30-target grid (deterministic given SAMPLE_SEED) ----
qs <- quantile(freqs, probs = c(1/3, 2/3), names = FALSE)
pool_low  <- sort(names(freqs)[freqs <= qs[1]])
pool_med  <- sort(names(freqs)[freqs >  qs[1] & freqs <= qs[2]])
pool_high <- sort(names(freqs)[freqs >  qs[2]])

if (any(c(length(pool_low), length(pool_med), length(pool_high)) < N_PER_BIN))
  stop("A frequency bin has fewer than ", N_PER_BIN, " mutations: ",
       sprintf("low=%d, med=%d, high=%d",
               length(pool_low), length(pool_med), length(pool_high)))

set.seed(SAMPLE_SEED)
targets_low  <- sort(sample(pool_low,  N_PER_BIN))
targets_med  <- sort(sample(pool_med,  N_PER_BIN))
targets_high <- sort(sample(pool_high, N_PER_BIN))

TARGETS_DF <- rbind(
  data.frame(target = targets_low,  tag = "low_sparsity",    stringsAsFactors = FALSE),
  data.frame(target = targets_med,  tag = "medium_sparsity", stringsAsFactors = FALSE),
  data.frame(target = targets_high, tag = "high_sparsity",   stringsAsFactors = FALSE)
)
TARGETS_DF$freq <- freqs[TARGETS_DF$target]

## Setup grid: (target, s); ordered (all targets at s=5) then (all at s=10).
SETUPS <- rbind(
  cbind(TARGETS_DF, s = S_VALUES[1L]),
  cbind(TARGETS_DF, s = S_VALUES[2L])
)
rownames(SETUPS) <- NULL
N_SETUPS <- nrow(SETUPS)            # 60

cat(sprintf("Tercile cutoffs: q1/3 = %.4f, q2/3 = %.4f\n", qs[1], qs[2]))
cat(sprintf("Pool sizes: low=%d, medium=%d, high=%d\n",
            length(pool_low), length(pool_med), length(pool_high)))
cat(sprintf("Sampled %d targets per bin; %d total setups (target x s).\n",
            N_PER_BIN, N_SETUPS))

# ---- 3. CLI: setup index 1..N_SETUPS ----
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) < 1L)
  stop("Usage: Rscript ccle-typeI-synthY-30targets.R <setup_idx in 1..", N_SETUPS, ">")
setup_idx <- suppressWarnings(as.integer(cli_args[[1L]]))
if (is.na(setup_idx) || setup_idx < 1L || setup_idx > N_SETUPS)
  stop("setup_idx must be an integer in 1..", N_SETUPS, " (got: ", cli_args[[1L]], ")")

target <- SETUPS$target[setup_idx]
tag    <- SETUPS$tag[setup_idx]
s_j    <- as.integer(SETUPS$s[setup_idx])
freq_j <- SETUPS$freq[setup_idx]

stopifnot(target %in% mut_names)
cat(sprintf("\n[Setup %d/%d]  target=%s  tag=%s  freq=%.4f  s=%d  R=%d\n",
            setup_idx, N_SETUPS, target, tag, freq_j, s_j, R))

# ---- 4. Synthetic Y generator (identical to ccle-typeI-synthY.R) ----
generate_Y_synth <- function(Z, s) {
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
  mu + rnorm(length(mu), mean = 0, sd = sigma)
}

# ---- 5. One-rep dispatcher (5 methods, same as parent script) ----
method_names <- c("asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive")

run_one_rep <- function(Y, X, Z) {
  out <- setNames(rep(NA_real_, length(method_names)), method_names)

  out["asym_gcm"] <- tryCatch({
    comets::gcm(Y = Y, X = X, Z = Z,
                reg_YonZ = "rf", reg_XonZ = "rf",
                alternative = "two.sided")$p.value
  }, error = function(e) NA_real_)

  out["asym_pcm"] <- tryCatch({
    comets::pcm(Y = Y, X = X, Z = Z, rep = PCM_REPS,
                reg_YonXZ = "rf", reg_YonZ = "rf",
                reg_YhatonZ = "rf", reg_VonXZ = "rf",
                reg_RonZ = "rf")$p.value
  }, error = function(e) NA_real_)

  out["perm_gcm"] <- tryCatch({
    pcmt::perm_gcm(Y = Y, X = X, Z = Z,
                   reg_YonZ = "rf", reg_XonZ = "rf",
                   B = B_PERM, h = H_THRESH, n_cores = N_CORES,
                   alternative = "two.sided")$p.value
  }, error = function(e) NA_real_)

  out["perm_pcm"] <- tryCatch({
    pcmt::perm_pcm(Y = Y, X = X, Z = Z,
                   reg_YonXZ = "rf", reg_YonZ = "rf",
                   reg_YhatonZ = "rf", reg_RonZ = "rf",
                   B = B_PERM, h = H_THRESH, n_cores = N_CORES,
                   alternative = "two.sided")$p.value
  }, error = function(e) NA_real_)

  out["naive"] <- tryCatch({
    suppressWarnings(cor.test(Y, X))$p.value
  }, error = function(e) NA_real_)

  out
}

# ---- 6. Run R replicates for this single setup ----
X_j <- M_full[, target]
Z_j <- M_full[, setdiff(mut_names, target), drop = FALSE]

set.seed(RUN_SEED_BASE + setup_idx)

pvals <- matrix(NA_real_, nrow = R, ncol = length(method_names),
                dimnames = list(NULL, method_names))

t0 <- Sys.time()
for (r in seq_len(R)) {
  Y_syn      <- generate_Y_synth(Z_j, s_j)
  pvals[r, ] <- run_one_rep(Y_syn, X_j, Z_j)

  if (r %% 25L == 0L || r == R) {
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- el / r * (R - r)
    cat(sprintf("  rep %4d/%d  elapsed=%6.0fs  ETA=%6.0fs  rej_so_far=%s\n",
                r, R, el, eta,
                paste(sprintf("%s=%.2f",
                              method_names,
                              colMeans(pvals[seq_len(r), , drop = FALSE] <= ALPHA,
                                       na.rm = TRUE)),
                      collapse = " ")))
  }
}

rej     <- colMeans(pvals <= ALPHA, na.rm = TRUE)
na_rate <- colMeans(is.na(pvals))
cat("  Final rejection rates:\n"); print(round(rej, 3))

result <- list(
  setup_idx      = setup_idx,
  target         = target,
  sparsity_tag   = tag,
  freq           = freq_j,
  s              = s_j,
  R              = R,
  pvals          = pvals,
  rejection_rate = rej,
  na_rate        = na_rate
)

out_file <- file.path(bpath,
  sprintf("ccle-typeI-synthY-30targets-results-%s-s%d.rds", target, s_j))
saveRDS(list(result), out_file)
cat("\nSaved:\n   ", out_file, "\n")
cat("Run combine-ccle-results.R after all ", N_SETUPS, " jobs finish.\n")
