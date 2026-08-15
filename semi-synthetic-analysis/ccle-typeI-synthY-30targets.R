#!/usr/bin/env Rscript
###############################################################################
# CCLE Semi-Synthetic Representative Power Analysis -- s = 10
#
# Companion to ccle-typeI-synthY.R.  Whereas that script ran 3 hand-picked
# targets, this one combines the full null sweep with paired alternatives on
# one preselected representative per sparsity group:
#
#   * 10 targets sampled from each frequency tercile of the screened mutations
#     (low / medium / high sparsity) at delta = 0 ->  30 null setups
#   * RIOK3_MUT effects in {1, 2, 3, 4}           ->   4 alternatives
#   * PLK1_MUT and ATP11A_MUT effects in
#     {.5, 1, 1.5, 2}                             ->   8 alternatives
#   * 30 + 12 = 42 setups, each with R = 500 reps
#
# The runner takes one CLI integer in 1..42. Setups 1:30 are the complete null
# sweep in the old target order. Setups 31:42 are ordered by within-target grid
# index, then RIOK3_MUT / PLK1_MUT / ATP11A_MUT. All effects for a target use
# the old matched-v2 s=10 seed index, pairing their generated outcomes, splits,
# permutation streams, and RF streams. Legacy files are left untouched:
#
#   $LOCAL_PDRT_DATA_DIR/ccle/
#     ccle-typeI-synthY-30targets-matched-v4-results-{TARGET}-{EFFECT_ID}-s10.rds
#
# After all 42 jobs finish, run combine-ccle-results.R to merge.
###############################################################################

suppressPackageStartupMessages({
  library(comets)
  library(pcmt)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not determine this script's path", call. = FALSE)
}
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "matched-calibration-utils.R"))

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
B_MAX     <- 499L
H_THRESH  <- 15L
ALPHA     <- 0.05
N_CORES   <- 1L
ALTERNATIVE <- "two.sided"
PCM_FRAC  <- 0.5
RF_NUM_TREES <- 500L
S_FIXED   <- 10L
EFFECT_SCALE <- "noise_sd_equals_sd_mu"
EFFECT_LOOKUP <- data.frame(
  effect_id = c(
    "delta0", "delta0p5", "delta1", "delta1p5", "delta2", "delta3", "delta4"
  ),
  effect_size = c(0, 0.5, 1, 1.5, 2, 3, 4),
  stringsAsFactors = FALSE
)
REPRESENTATIVE_TARGETS <- data.frame(
  representative_order = seq_len(3L),
  target = c("RIOK3_MUT", "PLK1_MUT", "ATP11A_MUT"),
  expected_tag = c("low_sparsity", "medium_sparsity", "high_sparsity"),
  expected_target_idx = c(9L, 15L, 22L),
  stringsAsFactors = FALSE
)
REPRESENTATIVE_EFFECT_SIZES <- rbind(
  RIOK3_MUT = c(0, 1, 2, 3, 4),
  PLK1_MUT = c(0, 0.5, 1, 1.5, 2),
  ATP11A_MUT = c(0, 0.5, 1, 1.5, 2)
)
colnames(REPRESENTATIVE_EFFECT_SIZES) <- paste0("grid", seq_len(5L))
N_PER_BIN <- 10L
SAMPLE_SEED <- 2026L      # seed for tercile sampling (must be stable across jobs)
RUN_SEED_BASE <- 20260000L # base for deterministic setup/rep/component seeds

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
TARGETS_DF$target_idx <- seq_len(nrow(TARGETS_DF))

representative_idx <- match(REPRESENTATIVE_TARGETS$target, TARGETS_DF$target)
if (anyNA(representative_idx)) {
  stop(
    "The preselected representative target(s) are absent from the sampled grid: ",
    paste(REPRESENTATIVE_TARGETS$target[is.na(representative_idx)], collapse = ", "),
    call. = FALSE
  )
}
REPRESENTATIVES_DF <- TARGETS_DF[representative_idx, , drop = FALSE]
REPRESENTATIVES_DF$representative_order <-
  REPRESENTATIVE_TARGETS$representative_order
if (!identical(REPRESENTATIVES_DF$target, REPRESENTATIVE_TARGETS$target) ||
    !identical(REPRESENTATIVES_DF$tag, REPRESENTATIVE_TARGETS$expected_tag) ||
    !identical(
      as.integer(REPRESENTATIVES_DF$target_idx),
      REPRESENTATIVE_TARGETS$expected_target_idx
    )) {
  stop(
    paste(
      "The preselected representatives do not match the expected sparsity",
      "groups and canonical target indices."
    ),
    call. = FALSE
  )
}

effect_id_for_size <- function(effect_size) {
  idx <- match(effect_size, EFFECT_LOOKUP$effect_size)
  if (anyNA(idx)) {
    stop("An effect size is absent from EFFECT_LOOKUP.", call. = FALSE)
  }
  EFFECT_LOOKUP$effect_id[idx]
}

## Preserve all 30 old s=10 null streams first. Then append nonzero effects in
## within-target grid-index-major order, with representatives displayed from
## high sparsity (rare) to low sparsity (common). Paired seed indices 31:60 are
## target-specific and therefore identical to the matched-v2 null schedule.
SETUPS_NULL <- transform(
  TARGETS_DF,
  grid_idx = 1L,
  effect_id = "delta0",
  effect_size = 0,
  hypothesis = "null",
  effect_scale = EFFECT_SCALE,
  s = S_FIXED,
  paired_seed_idx = target_idx + nrow(TARGETS_DF)
)
SETUPS_ALT <- do.call(rbind, lapply(2:5, function(grid_idx) {
  effect_size <- as.numeric(
    REPRESENTATIVE_EFFECT_SIZES[REPRESENTATIVES_DF$target, grid_idx]
  )
  transform(
    REPRESENTATIVES_DF,
    grid_idx = as.integer(grid_idx),
    effect_id = effect_id_for_size(effect_size),
    effect_size = effect_size,
    hypothesis = "alternative",
    effect_scale = EFFECT_SCALE,
    s = S_FIXED,
    paired_seed_idx = target_idx + nrow(TARGETS_DF)
  )
}))
SETUPS <- rbind(
  SETUPS_NULL,
  SETUPS_ALT[, names(SETUPS_NULL), drop = FALSE]
)
rownames(SETUPS) <- NULL
N_SETUPS <- nrow(SETUPS)            # 30 nulls + 3 representatives x 4 alternatives
stopifnot(
  N_SETUPS == 42L,
  identical(as.integer(SETUPS$grid_idx), c(rep(1L, 30L), rep(2:5, each = 3L))),
  identical(SETUPS$target[1:30], TARGETS_DF$target),
  identical(SETUPS$effect_id[1:30], rep("delta0", 30L)),
  identical(as.numeric(SETUPS$effect_size[1:30]), rep(0, 30L)),
  identical(SETUPS$hypothesis[1:30], rep("null", 30L)),
  identical(
    SETUPS$target[31:42], rep(REPRESENTATIVE_TARGETS$target, 4L)
  ),
  identical(
    SETUPS$tag[31:42], rep(REPRESENTATIVE_TARGETS$expected_tag, 4L)
  ),
  identical(
    as.integer(SETUPS$target_idx[31:42]), rep(representative_idx, 4L)
  ),
  identical(as.integer(representative_idx), c(9L, 15L, 22L)),
  identical(
    as.numeric(SETUPS$effect_size[31:42]),
    as.numeric(REPRESENTATIVE_EFFECT_SIZES[, 2:5, drop = FALSE])
  ),
  identical(
    SETUPS$effect_id[31:42], effect_id_for_size(SETUPS$effect_size[31:42])
  ),
  identical(SETUPS$hypothesis[31:42], rep("alternative", 12L)),
  identical(as.integer(SETUPS$paired_seed_idx[1:30]), 31:60),
  identical(
    as.integer(SETUPS$paired_seed_idx[31:42]),
    rep(representative_idx + 30L, 4L)
  ),
  identical(
    as.integer(SETUPS$paired_seed_idx[seq(31L, 42L, by = 3L)]),
    rep(39L, 4L)
  ),
  all(as.integer(SETUPS$s) == S_FIXED),
  all(SETUPS$effect_scale == EFFECT_SCALE)
)

cat(sprintf("Tercile cutoffs: q1/3 = %.4f, q2/3 = %.4f\n", qs[1], qs[2]))
cat(sprintf("Pool sizes: low=%d, medium=%d, high=%d\n",
            length(pool_low), length(pool_med), length(pool_high)))
cat(sprintf("Sampled %d targets per bin; %d total null/representative setups.\n",
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
target_idx <- as.integer(SETUPS$target_idx[setup_idx])
grid_idx <- as.integer(SETUPS$grid_idx[setup_idx])
effect_id <- SETUPS$effect_id[setup_idx]
hypothesis <- SETUPS$hypothesis[setup_idx]
effect_size <- as.numeric(SETUPS$effect_size[setup_idx])
effect_scale <- SETUPS$effect_scale[setup_idx]
paired_seed_idx <- as.integer(SETUPS$paired_seed_idx[setup_idx])

stopifnot(target %in% mut_names)
cat(sprintf(
  paste0(
    "\n[Setup %d/%d]  target_idx=%d  target=%s  tag=%s  freq=%.4f  ",
    "s=%d  grid_idx=%d  effect_id=%s  hypothesis=%s  ",
    "effect_size=%.2f  paired_seed_idx=%d  R=%d\n"
  ),
  setup_idx, N_SETUPS, target_idx, target, tag, freq_j, s_j, grid_idx,
  effect_id, hypothesis, effect_size, paired_seed_idx, R
))

# ---- 4. Paired synthetic-outcome generator ----
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
  list(
    mu = mu,
    eps = rnorm(length(mu), mean = 0, sd = sigma),
    outcome_noise_sd = sigma
  )
}

# ---- 5. One-rep dispatcher: one statistic, two calibrations per method ----
method_names <- c("asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive")

run_one_rep <- function(Y, X, Z, seeds, pcm_D2_indices,
                        test_order = c("gcm", "pcm")) {
  matched <- ccle_run_matched_calibration(
    Y = Y,
    X = X,
    Z = Z,
    seeds = seeds,
    pcm_D2_indices = pcm_D2_indices,
    B_max = B_MAX,
    adaptive_h = H_THRESH,
    alpha = ALPHA,
    n_cores = N_CORES,
    alternative = ALTERNATIVE,
    pcm_frac = PCM_FRAC,
    rf_num_trees = RF_NUM_TREES,
    test_order = test_order
  )

  out <- setNames(rep(NA_real_, length(method_names)), method_names)
  out[names(matched$pvalues)] <- matched$pvalues
  out[["naive"]] <- tryCatch({
    suppressWarnings(cor.test(Y, X))$p.value
  }, error = function(e) NA_real_)

  list(
    pvalues = out,
    statistics = matched$statistics,
    diagnostics = matched$diagnostics
  )
}

# ---- 6. Run R replicates for this single setup ----
X_j <- M_full[, target]
Z_j <- M_full[, setdiff(mut_names, target), drop = FALSE]

pvals <- matrix(NA_real_, nrow = R, ncol = length(method_names),
                dimnames = list(NULL, method_names))
statistics <- matrix(
  NA_real_,
  nrow = R,
  ncol = 2L,
  dimnames = list(NULL, c("gcm", "pcm"))
)

seed_names <- names(ccle_seed_schedule(paired_seed_idx, 1L, RUN_SEED_BASE))
seed_metadata <- matrix(
  NA_integer_,
  nrow = R,
  ncol = length(seed_names),
  dimnames = list(NULL, seed_names)
)
pcm_D2_indices <- vector("list", R)
pcm_D1_n <- rep(NA_integer_, R)
pcm_D2_n <- rep(NA_integer_, R)
outcome_noise_sd <- rep(NA_real_, R)
gcm_actual_B <- rep(NA_integer_, R)
pcm_actual_B <- rep(NA_integer_, R)
gcm_stop_reason <- rep(NA_character_, R)
pcm_stop_reason <- rep(NA_character_, R)
gcm_error <- rep(NA_character_, R)
pcm_error <- rep(NA_character_, R)
gcm_matched_statistic <- rep(FALSE, R)
pcm_matched_statistic <- rep(FALSE, R)
gcm_loss_count <- rep(NA_integer_, R)
pcm_loss_count <- rep(NA_integer_, R)
gcm_t_star <- rep(NA_real_, R)
pcm_t_star <- rep(NA_real_, R)
gcm_stopped_early <- rep(NA, R)
pcm_stopped_early <- rep(NA, R)
gcm_adaptive_backend <- rep(NA_character_, R)
pcm_adaptive_backend <- rep(NA_character_, R)

t0 <- Sys.time()
for (r in seq_len(R)) {
  seeds <- ccle_seed_schedule(paired_seed_idx, r, RUN_SEED_BASE)
  seed_metadata[r, ] <- seeds

  # Every random component has its own local stream.  Changing method order or
  # the number of permutations consumed by one method cannot affect another.
  # Every effect for a target reuses every stream; its Y values differ only by
  # the deterministic carrier effect below.
  outcome <- ccle_with_seed(
    seeds[["data_seed"]],
    generate_Y_components(Z_j, s_j)
  )
  outcome_noise_sd[[r]] <- outcome$outcome_noise_sd
  Y_syn <- outcome$mu + outcome$eps +
    effect_size * outcome$outcome_noise_sd * X_j
  split_D2 <- ccle_make_pcm_split(
    n = n,
    split_seed = seeds[["pcm_split_seed"]],
    frac = PCM_FRAC
  )
  pcm_D2_indices[[r]] <- split_D2
  pcm_D2_n[[r]] <- length(split_D2)
  pcm_D1_n[[r]] <- n - length(split_D2)

  rep_result <- run_one_rep(
    Y = Y_syn,
    X = X_j,
    Z = Z_j,
    seeds = seeds,
    pcm_D2_indices = split_D2
  )
  rep_complete <-
    all(is.finite(rep_result$pvalues)) &&
    all(is.finite(rep_result$statistics)) &&
    isTRUE(rep_result$diagnostics$gcm$matched_statistic) &&
    isTRUE(rep_result$diagnostics$pcm$matched_statistic) &&
    is.na(rep_result$diagnostics$gcm$error) &&
    is.na(rep_result$diagnostics$pcm$error)
  if (!rep_complete) {
    stop(
      sprintf(
        paste0(
          "Matched-calibration replicate %d failed; no partial v4 result ",
          "will be saved. GCM error: %s; PCM error: %s"
        ),
        r,
        ifelse(
          is.na(rep_result$diagnostics$gcm$error),
          "none",
          rep_result$diagnostics$gcm$error
        ),
        ifelse(
          is.na(rep_result$diagnostics$pcm$error),
          "none",
          rep_result$diagnostics$pcm$error
        )
      ),
      call. = FALSE
    )
  }
  pvals[r, ] <- rep_result$pvalues
  statistics[r, ] <- rep_result$statistics

  gcm_actual_B[[r]] <- rep_result$diagnostics$gcm$actual_B
  pcm_actual_B[[r]] <- rep_result$diagnostics$pcm$actual_B
  gcm_stop_reason[[r]] <- rep_result$diagnostics$gcm$stop_reason
  pcm_stop_reason[[r]] <- rep_result$diagnostics$pcm$stop_reason
  gcm_error[[r]] <- rep_result$diagnostics$gcm$error
  pcm_error[[r]] <- rep_result$diagnostics$pcm$error
  gcm_matched_statistic[[r]] <- rep_result$diagnostics$gcm$matched_statistic
  pcm_matched_statistic[[r]] <- rep_result$diagnostics$pcm$matched_statistic
  gcm_loss_count[[r]] <- rep_result$diagnostics$gcm$loss_count
  pcm_loss_count[[r]] <- rep_result$diagnostics$pcm$loss_count
  gcm_t_star[[r]] <- rep_result$diagnostics$gcm$t_star
  pcm_t_star[[r]] <- rep_result$diagnostics$pcm$t_star
  gcm_stopped_early[[r]] <- rep_result$diagnostics$gcm$stopped_early
  pcm_stopped_early[[r]] <- rep_result$diagnostics$pcm$stopped_early
  gcm_adaptive_backend[[r]] <- rep_result$diagnostics$gcm$adaptive_backend
  pcm_adaptive_backend[[r]] <- rep_result$diagnostics$pcm$adaptive_backend

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

replication_metadata <- data.frame(
  rep = seq_len(R),
  as.data.frame(seed_metadata, optional = TRUE),
  outcome_noise_sd = outcome_noise_sd,
  pcm_D1_n = pcm_D1_n,
  pcm_D2_n = pcm_D2_n,
  gcm_actual_B = gcm_actual_B,
  pcm_actual_B = pcm_actual_B,
  gcm_stop_reason = gcm_stop_reason,
  pcm_stop_reason = pcm_stop_reason,
  gcm_error = gcm_error,
  pcm_error = pcm_error,
  gcm_matched_statistic = gcm_matched_statistic,
  pcm_matched_statistic = pcm_matched_statistic,
  gcm_loss_count = gcm_loss_count,
  pcm_loss_count = pcm_loss_count,
  gcm_t_star = gcm_t_star,
  pcm_t_star = pcm_t_star,
  gcm_stopped_early = gcm_stopped_early,
  pcm_stopped_early = pcm_stopped_early,
  gcm_adaptive_backend = gcm_adaptive_backend,
  pcm_adaptive_backend = pcm_adaptive_backend,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
replication_metadata$pcm_D2_indices <- I(pcm_D2_indices)
all_replicates_matched <-
  all(replication_metadata$gcm_matched_statistic) &&
  all(replication_metadata$pcm_matched_statistic) &&
  all(is.na(replication_metadata$gcm_error)) &&
  all(is.na(replication_metadata$pcm_error)) &&
  all(is.finite(pvals)) &&
  all(is.finite(statistics))
if (!all_replicates_matched) {
  stop("Refusing to save an incomplete or unmatched v4 result", call. = FALSE)
}

package_versions <- c(
  R = R.version.string,
  pcmt = as.character(utils::packageVersion("pcmt")),
  comets = as.character(utils::packageVersion("comets")),
  ranger = as.character(utils::packageVersion("ranger"))
)

result <- list(
  schema_version  = "ccle_synthy_matched_v4",
  simulation_design = "representative_effect_grid_s10",
  calibration_design = "shared_observed_statistic",
  setup_idx      = setup_idx,
  target_idx     = target_idx,
  paired_seed_idx = paired_seed_idx,
  target         = target,
  sparsity_tag   = tag,
  freq           = freq_j,
  s              = s_j,
  grid_idx       = grid_idx,
  effect_id      = effect_id,
  hypothesis     = hypothesis,
  effect_size    = effect_size,
  effect_scale   = effect_scale,
  R              = R,
  B_requested    = B_MAX,
  B_max          = B_MAX,
  adaptive_h     = H_THRESH,
  adaptive_alpha = ALPHA,
  alternative    = ALTERNATIVE,
  pvals          = pvals,
  statistics     = statistics,
  replication_metadata = replication_metadata,
  matched_calibration = all_replicates_matched,
  method_definition = list(
    alternative = ALTERNATIVE,
    regression = "rf",
    B_requested = B_MAX,
    B_max = B_MAX,
    h = H_THRESH,
    alpha = ALPHA,
    permutation_mode = "adaptive_anytime_valid",
    n_cores = N_CORES,
    pcm_frac = PCM_FRAC,
    rf_num_trees = RF_NUM_TREES,
    rf_num_threads = 1L,
    role_specific_rf_seeds = TRUE,
    gcm = paste(
      "One pcmt::perm_gcm fit; fit$statistic is normal-calibrated and",
      "fit$p.value is permutation-calibrated."
    ),
    pcm = paste(
      "One pcmt::perm_pcm fit with explicit D2 indices; fit$statistic",
      "is normal-calibrated and fit$p.value is permutation-calibrated."
    ),
    seed_schedule = list(
      base = RUN_SEED_BASE,
      formula = "base + paired_seed_idx * 100000 + rep * 100 + component_offset",
      schedule_index = "paired_seed_idx = target_idx + 30",
      setup_order = "null_all_targets_then_grid_index_major_representatives",
      paired_across_effect_sizes = TRUE,
      component_offsets = CCLE_SEED_COMPONENT_OFFSETS
    )
  ),
  package_versions = package_versions,
  rejection_rate = rej,
  na_rate        = na_rate
)

out_file <- file.path(bpath,
  sprintf(
    "ccle-typeI-synthY-30targets-matched-v4-results-%s-%s-s10.rds",
    target,
    effect_id
  ))
tmp_file <- tempfile(
  pattern = paste0(basename(out_file), ".tmp-"),
  tmpdir = dirname(out_file)
)
on.exit(unlink(tmp_file), add = TRUE)
saveRDS(list(result), tmp_file)
if (!file.rename(tmp_file, out_file)) {
  stop("Failed to atomically move v4 result into place: ", out_file,
       call. = FALSE)
}
cat("\nSaved:\n   ", out_file, "\n")
cat("Run combine-ccle-results.R after all ", N_SETUPS, " jobs finish.\n")
