#!/usr/bin/env Rscript
###############################################################################
# Merge the paired matched-calibration v4 representative-target effect grid
# from ccle-typeI-synthY-30targets.R into a single versioned combined file.
#
# Inputs and outputs live in $LOCAL_PDRT_DATA_DIR/ccle:
#   ccle-typeI-synthY-30targets-matched-v4-results-{TARGET}-{EFFECT_ID}-s10.rds
#   ccle-typeI-synthY-30targets-matched-v4-results.rds
#   ccle-typeI-synthY-30targets-matched-v4-rejection-rates.csv
#
# To know which setups are expected, this script re-derives the same 30
# deterministic targets as the runner: it samples 10 targets from each
# frequency tercile of the screened mutations using SAMPLE_SEED = 2026. All 30
# targets have a delta0 job. Power is evaluated only for three pre-specified
# representatives, in high-/medium-/low-sparsity display order: EIF5A_MUT on
# {0, 1, 2, 3, 4}, and PLK1_MUT / ATP11A_MUT on {0, .5, 1, 1.5, 2}. Setup IDs
# 31:42 are grid-index-major across those three representatives.
#
# Usage: Rscript combine-ccle-results.R
###############################################################################

FIXED_S     <- 10L
N_PER_BIN   <- 10L
SAMPLE_SEED <- 2026L
RUN_SEED_BASE <- 20260000L
EXPECTED_R  <- 500L
EXPECTED_B_MAX <- 499L
EXPECTED_H <- 15L
EXPECTED_ALPHA <- 0.05
EXPECTED_PCM_FRAC <- 0.5
EXPECTED_RF_NUM_TREES <- 500L
SCHEMA_VERSION     <- "ccle_synthy_matched_v4"
SIMULATION_DESIGN  <- "representative_effect_grid_s10"
CALIBRATION_DESIGN <- "shared_observed_statistic"
EFFECT_SCALE       <- "noise_sd_equals_sd_mu"
SETUP_ORDER <- "null_all_targets_then_grid_index_major_representatives"
REPRESENTATIVE_EXPECTATIONS <- data.frame(
  representative_order = seq_len(3L),
  target = c("EIF5A_MUT", "PLK1_MUT", "ATP11A_MUT"),
  sparsity_tag = c("low_sparsity", "medium_sparsity", "high_sparsity"),
  stringsAsFactors = FALSE
)
REPRESENTATIVE_EFFECT_SIZES <- list(
  EIF5A_MUT = c(0, 1, 2, 3, 4),
  PLK1_MUT = c(0, 0.5, 1, 1.5, 2),
  ATP11A_MUT = c(0, 0.5, 1, 1.5, 2)
)
effect_id_from_size <- function(x) {
  lookup <- c(
    `0` = "delta0", `0.5` = "delta0p5", `1` = "delta1",
    `1.5` = "delta1p5", `2` = "delta2", `3` = "delta3",
    `4` = "delta4"
  )
  ids <- unname(lookup[as.character(x)])
  if (anyNA(ids)) {
    stop("No canonical effect ID for effect size: ",
         paste(x[is.na(ids)], collapse = ", "), call. = FALSE)
  }
  ids
}
EXPECTED_METHODS <- c(
  "asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive"
)
EXPECTED_STATISTICS <- c("gcm", "pcm")
EXPECTED_SEED_COMPONENT_OFFSETS <- c(
  data_seed = 1L,
  pcm_split_seed = 2L,
  gcm_permutation_seed = 3L,
  pcm_permutation_seed = 4L,
  gcm_YonZ_seed = 5L,
  gcm_XonZ_seed = 6L,
  pcm_YonXZ_seed = 7L,
  pcm_YonZ_seed = 8L,
  pcm_YhatonZ_seed = 9L,
  pcm_RonZ_seed = 10L
)

cli_args <- commandArgs(trailingOnly = TRUE)
status_mode <- FALSE
status_path <- NULL
status_indices <- integer()
if (length(cli_args)) {
  if (length(cli_args) != 3L ||
      !identical(cli_args[[1L]], "--output-status")) {
    stop(
      "Usage: Rscript combine-ccle-results.R ",
      "[--output-status <path> <n>]",
      call. = FALSE
    )
  }
  index_text <- cli_args[[3L]]
  if (!grepl("^[1-9][0-9]*$", index_text)) {
    stop("Status setup count/index must be a canonical integer in 1..42",
         call. = FALSE)
  }
  index_value <- as.integer(index_text)
  if (is.na(index_value) || index_value > 42L) {
    stop("Status setup count/index must be a canonical integer in 1..42",
         call. = FALSE)
  }
  status_mode <- TRUE
  status_path <- path.expand(cli_args[[2L]])
  if (!dir.exists(dirname(status_path))) {
    stop("Status-output directory does not exist: ", dirname(status_path),
         call. = FALSE)
  }
  status_indices <- seq_len(index_value)
}

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
bpath <- file.path(data_root, "ccle")

# ---- Re-derive expected (target, s) grid ----
in_dir <- bpath
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

M_full <- as.matrix(dat[, setdiff(colnames(dat), "Y_real")])
freqs  <- colMeans(M_full)

qs <- quantile(freqs, probs = c(1/3, 2/3), names = FALSE)
pool_low  <- sort(names(freqs)[freqs <= qs[1]])
pool_med  <- sort(names(freqs)[freqs >  qs[1] & freqs <= qs[2]])
pool_high <- sort(names(freqs)[freqs >  qs[2]])

if (any(c(length(pool_low), length(pool_med), length(pool_high)) < N_PER_BIN)) {
  stop(
    "A frequency bin has fewer than ", N_PER_BIN, " mutations: ",
    sprintf(
      "low=%d, med=%d, high=%d",
      length(pool_low), length(pool_med), length(pool_high)
    ),
    call. = FALSE
  )
}

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

representative_match <- match(
  REPRESENTATIVE_EXPECTATIONS$target,
  TARGETS_DF$target
)
if (anyNA(representative_match)) {
  stop(
    "Pre-specified representative target(s) are absent from the canonical ",
    "30-target sample: ",
    paste(
      REPRESENTATIVE_EXPECTATIONS$target[is.na(representative_match)],
      collapse = ", "
    ),
    call. = FALSE
  )
}
REPRESENTATIVE_TARGETS <- data.frame(
  representative_order = REPRESENTATIVE_EXPECTATIONS$representative_order,
  target_idx = TARGETS_DF$target_idx[representative_match],
  target = TARGETS_DF$target[representative_match],
  sparsity_tag = TARGETS_DF$tag[representative_match],
  freq = TARGETS_DF$freq[representative_match],
  stringsAsFactors = FALSE
)
if (!identical(
      REPRESENTATIVE_TARGETS$sparsity_tag,
      REPRESENTATIVE_EXPECTATIONS$sparsity_tag
    )) {
  stop(
    "Representative targets do not have their expected canonical sparsity ",
    "tags.",
    call. = FALSE
  )
}

## The 15-row metadata grid is grid-index-major, matching the ordering of the
## representative records once their already-existing delta0 rows are selected.
REPRESENTATIVE_EFFECT_GRID <- do.call(rbind, lapply(seq_len(5L), function(i) {
  effect_sizes <- vapply(
    REPRESENTATIVE_TARGETS$target,
    function(target) REPRESENTATIVE_EFFECT_SIZES[[target]][[i]],
    numeric(1L)
  )
  transform(
    REPRESENTATIVE_TARGETS,
    grid_idx = i,
    effect_id = effect_id_from_size(effect_sizes),
    effect_size = effect_sizes,
    hypothesis = ifelse(effect_sizes == 0, "null", "alternative")
  )
}))
rownames(REPRESENTATIVE_EFFECT_GRID) <- NULL
REPRESENTATIVE_EFFECT_GRID <- REPRESENTATIVE_EFFECT_GRID[, c(
  "representative_order", "target_idx", "target", "sparsity_tag", "freq",
  "grid_idx", "effect_id", "effect_size", "hypothesis"
)]

null_setups <- transform(
  TARGETS_DF,
  grid_idx = 1L,
  effect_id = "delta0",
  effect_size = 0,
  hypothesis = "null"
)
names(null_setups)[names(null_setups) == "tag"] <- "sparsity_tag"
representative_alternatives <- REPRESENTATIVE_EFFECT_GRID[
  REPRESENTATIVE_EFFECT_GRID$grid_idx > 1L,
  c(
    "target", "sparsity_tag", "freq", "target_idx", "grid_idx",
    "effect_id", "effect_size", "hypothesis"
  )
]
SETUPS <- rbind(
  null_setups[, names(representative_alternatives)],
  representative_alternatives
)
SETUPS$s <- FIXED_S
SETUPS$paired_seed_idx <- SETUPS$target_idx + nrow(TARGETS_DF)
SETUPS$setup_idx <- seq_len(nrow(SETUPS))
SETUPS <- SETUPS[, c(
  "setup_idx", "target_idx", "paired_seed_idx", "target", "sparsity_tag",
  "freq", "s", "grid_idx", "effect_id", "effect_size", "hypothesis"
)]

stopifnot(
  nrow(TARGETS_DF) == 30L,
  nrow(REPRESENTATIVE_TARGETS) == 3L,
  nrow(REPRESENTATIVE_EFFECT_GRID) == 15L,
  nrow(SETUPS) == 42L,
  identical(as.integer(SETUPS$setup_idx), seq_len(42L)),
  identical(as.integer(SETUPS$target_idx[seq_len(30L)]), seq_len(30L)),
  identical(as.integer(SETUPS$grid_idx[seq_len(30L)]), rep(1L, 30L)),
  identical(SETUPS$effect_id[seq_len(30L)], rep("delta0", 30L)),
  identical(as.numeric(SETUPS$effect_size[seq_len(30L)]), rep(0, 30L)),
  identical(
    as.integer(SETUPS$grid_idx[31:42]),
    rep(2:5, each = nrow(REPRESENTATIVE_TARGETS))
  ),
  identical(
    SETUPS$target[31:42],
    rep(REPRESENTATIVE_TARGETS$target, times = 4L)
  ),
  identical(
    as.integer(SETUPS$paired_seed_idx),
    as.integer(SETUPS$target_idx + nrow(TARGETS_DF))
  )
)

# ---- Load each per-setup file ----
expected_files <- file.path(
  bpath,
  sprintf(
    "ccle-typeI-synthY-30targets-matched-v4-results-%s-%s-s%d.rds",
    SETUPS$target,
    SETUPS$effect_id,
    SETUPS$s
  )
)
missing <- expected_files[!file.exists(expected_files)]
if (length(missing) && !status_mode) {
  stop(
    "Missing paired matched-calibration v4 per-setup files (", length(missing),
    " of ", length(expected_files), "):\n  ",
    paste(head(missing, 20), collapse = "\n  "),
    if (length(missing) > 20) "\n  ..." else "",
    call. = FALSE
  )
}

actual_v4_files <- list.files(
  bpath,
  pattern = "^ccle-typeI-synthY-30targets-matched-v4-results-.+[.]rds$",
  full.names = TRUE
)
unexpected <- setdiff(
  basename(actual_v4_files),
  basename(expected_files)
)
if (length(unexpected)) {
  stop(
    "Unexpected matched-v4 per-setup file(s); the directory must contain ",
    "exactly the canonical 42-file null-plus-representative grid:\n  ",
    paste(
      head(file.path(bpath, unexpected), 20L),
      collapse = "\n  "
    ),
    if (length(unexpected) > 20L) "\n  ..." else "",
    call. = FALSE
  )
}

results <- list()

required_record_fields <- c(
  "schema_version", "simulation_design", "calibration_design",
  "matched_calibration", "setup_idx", "target_idx", "paired_seed_idx",
  "target", "sparsity_tag", "freq", "s", "grid_idx", "effect_id",
  "effect_size", "hypothesis", "effect_scale", "R",
  "B_requested", "B_max", "adaptive_h", "adaptive_alpha", "alternative",
  "pvals", "statistics",
  "replication_metadata", "method_definition", "package_versions",
  "rejection_rate", "na_rate"
)
required_metadata_fields <- c(
  "rep", "data_seed", "pcm_split_seed", "gcm_permutation_seed",
  "pcm_permutation_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "pcm_YonXZ_seed", "pcm_YonZ_seed", "pcm_YhatonZ_seed",
  "pcm_RonZ_seed", "outcome_noise_sd", "pcm_D1_n", "pcm_D2_n", "gcm_actual_B",
  "pcm_actual_B", "gcm_stop_reason", "pcm_stop_reason", "gcm_error",
  "pcm_error", "gcm_matched_statistic", "pcm_matched_statistic",
  "gcm_loss_count", "pcm_loss_count", "gcm_t_star", "pcm_t_star",
  "gcm_stopped_early", "pcm_stopped_early", "gcm_adaptive_backend",
  "pcm_adaptive_backend", "pcm_D2_indices"
)
seed_fields <- c(
  "data_seed", "pcm_split_seed", "gcm_permutation_seed",
  "pcm_permutation_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "pcm_YonXZ_seed", "pcm_YonZ_seed", "pcm_YhatonZ_seed",
  "pcm_RonZ_seed"
)
stopifnot(identical(seed_fields, names(EXPECTED_SEED_COMPONENT_OFFSETS)))

is_integerish <- function(x) {
  is.numeric(x) && all(is.finite(x)) &&
    all(x == floor(x)) && all(x >= 1) && all(x <= .Machine$integer.max)
}

validate_record <- function(rec, setup_idx, path) {
  record_fields <- if (is.list(rec)) names(rec) else character()
  missing_fields <- setdiff(required_record_fields, record_fields)
  unexpected_fields <- setdiff(record_fields, required_record_fields)
  if (!is.list(rec) || length(missing_fields) || length(unexpected_fields)) {
    stop(
      "Malformed paired matched-calibration v4 record in ", path,
      if (length(missing_fields)) paste0(
        "; missing field(s): ", paste(missing_fields, collapse = ", ")
      ) else "",
      if (length(unexpected_fields)) paste0(
        "; unexpected field(s): ", paste(unexpected_fields, collapse = ", ")
      ) else "",
      call. = FALSE
    )
  }
  if (!identical(rec$schema_version, SCHEMA_VERSION) ||
      !identical(rec$simulation_design, SIMULATION_DESIGN) ||
      !identical(rec$calibration_design, CALIBRATION_DESIGN) ||
      !isTRUE(rec$matched_calibration)) {
    stop(
      "Wrong-schema, wrong-design, or unmatched calibration record rejected: ",
      path,
      call. = FALSE
    )
  }
  if (!identical(as.integer(rec$setup_idx), as.integer(setup_idx)) ||
      !identical(
        as.integer(rec$target_idx),
        as.integer(SETUPS$target_idx[setup_idx])
      ) ||
      !identical(
        as.integer(rec$paired_seed_idx),
        as.integer(SETUPS$paired_seed_idx[setup_idx])
      ) ||
      !identical(as.character(rec$target), as.character(SETUPS$target[setup_idx])) ||
      !identical(as.integer(rec$s), as.integer(SETUPS$s[setup_idx])) ||
      !identical(
        as.character(rec$sparsity_tag),
        as.character(SETUPS$sparsity_tag[setup_idx])
      ) ||
      !isTRUE(all.equal(
        as.numeric(rec$freq), as.numeric(SETUPS$freq[setup_idx]),
        tolerance = 1e-12
      ))) {
    stop("Setup metadata mismatch in ", path, call. = FALSE)
  }
  expected_grid_idx <- as.integer(SETUPS$grid_idx[setup_idx])
  expected_effect_id <- as.character(SETUPS$effect_id[setup_idx])
  expected_hypothesis <- SETUPS$hypothesis[setup_idx]
  expected_effect <- SETUPS$effect_size[setup_idx]
  if (!identical(as.integer(rec$grid_idx), expected_grid_idx) ||
      !identical(as.character(rec$effect_id), expected_effect_id) ||
      !identical(as.character(rec$hypothesis), expected_hypothesis) ||
      length(rec$effect_size) != 1L ||
      !is.numeric(rec$effect_size) || !is.finite(rec$effect_size) ||
      !isTRUE(all.equal(
        as.numeric(rec$effect_size), expected_effect, tolerance = 0
      )) ||
      !identical(rec$effect_scale, EFFECT_SCALE)) {
    stop("Representative-grid metadata mismatch in ", path, call. = FALSE)
  }
  if (!identical(as.integer(rec$R), EXPECTED_R) ||
      !identical(as.integer(rec$B_requested), EXPECTED_B_MAX) ||
      !identical(as.integer(rec$B_max), EXPECTED_B_MAX) ||
      !identical(as.integer(rec$adaptive_h), EXPECTED_H) ||
      !isTRUE(all.equal(
        as.numeric(rec$adaptive_alpha), EXPECTED_ALPHA, tolerance = 0
      )) ||
      !identical(rec$alternative, "two.sided")) {
    stop("Replication or calibration configuration mismatch in ", path,
         call. = FALSE)
  }
  if (!is.matrix(rec$pvals) || !is.numeric(rec$pvals) ||
      !identical(dim(rec$pvals), c(EXPECTED_R, length(EXPECTED_METHODS))) ||
      !identical(colnames(rec$pvals), EXPECTED_METHODS) ||
      any(!is.finite(rec$pvals)) || any(rec$pvals < 0 | rec$pvals > 1)) {
    stop("Invalid or missing v4 p-value matrix in ", path, call. = FALSE)
  }
  if (!is.matrix(rec$statistics) || !is.numeric(rec$statistics) ||
      !identical(dim(rec$statistics), c(EXPECTED_R, 2L)) ||
      !identical(colnames(rec$statistics), EXPECTED_STATISTICS) ||
      any(!is.finite(rec$statistics))) {
    stop("Invalid or missing shared observed statistics in ", path,
         call. = FALSE)
  }

  expected_asym_gcm <- 2 * stats::pnorm(
    abs(rec$statistics[, "gcm"]), lower.tail = FALSE
  )
  expected_asym_pcm <- 2 * stats::pnorm(
    abs(rec$statistics[, "pcm"]), lower.tail = FALSE
  )
  if (!isTRUE(all.equal(
        unname(rec$pvals[, "asym_gcm"]), unname(expected_asym_gcm),
        tolerance = 1e-14, check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        unname(rec$pvals[, "asym_pcm"]), unname(expected_asym_pcm),
        tolerance = 1e-14, check.attributes = FALSE
      ))) {
    stop(
      "Normal p-values are not calibrated from the stored shared statistic in ",
      path,
      call. = FALSE
    )
  }

  metadata <- rec$replication_metadata
  missing_metadata <- if (is.data.frame(metadata)) {
    setdiff(required_metadata_fields, names(metadata))
  } else {
    required_metadata_fields
  }
  if (!is.data.frame(metadata) || nrow(metadata) != EXPECTED_R ||
      length(missing_metadata) ||
      length(setdiff(names(metadata), required_metadata_fields))) {
    stop(
      "Invalid replication metadata in ", path,
      if (length(missing_metadata)) paste0(
        "; missing column(s): ", paste(missing_metadata, collapse = ", ")
      ) else "",
      call. = FALSE
    )
  }
  if (!identical(as.integer(metadata$rep), seq_len(EXPECTED_R))) {
    stop("Replication IDs are not 1:R in ", path, call. = FALSE)
  }
  if (!is.numeric(metadata$outcome_noise_sd) ||
      any(!is.finite(metadata$outcome_noise_sd)) ||
      any(metadata$outcome_noise_sd <= 0)) {
    stop("Invalid outcome-noise scale metadata in ", path, call. = FALSE)
  }
  bad_seed <- vapply(seed_fields, function(field) {
    !is_integerish(metadata[[field]]) || anyDuplicated(metadata[[field]]) > 0L
  }, logical(1L))
  if (any(bad_seed)) {
    stop(
      "Invalid or repeated seed stream(s) in ", path, ": ",
      paste(seed_fields[bad_seed], collapse = ", "),
      call. = FALSE
    )
  }
  expected_seed_matrix <- vapply(
    seed_fields,
    function(field) {
      as.integer(
        RUN_SEED_BASE +
          as.double(rec$paired_seed_idx) * 100000 +
          as.double(seq_len(EXPECTED_R)) * 100 +
          EXPECTED_SEED_COMPONENT_OFFSETS[[field]]
      )
    },
    integer(EXPECTED_R)
  )
  actual_seed_matrix <- vapply(
    seed_fields,
    function(field) as.integer(metadata[[field]]),
    integer(EXPECTED_R)
  )
  if (!identical(actual_seed_matrix, expected_seed_matrix)) {
    stop("Seed metadata disagrees with the paired v4 schedule in ", path,
         call. = FALSE)
  }
  if (!is_integerish(metadata$pcm_D1_n) ||
      !is_integerish(metadata$pcm_D2_n) ||
      !is.list(metadata$pcm_D2_indices)) {
    stop("Invalid PCM split metadata in ", path, call. = FALSE)
  }
  valid_indices <- vapply(seq_len(EXPECTED_R), function(rep_idx) {
    indices <- metadata$pcm_D2_indices[[rep_idx]]
    is.integer(indices) && !anyNA(indices) &&
      length(indices) == metadata$pcm_D2_n[rep_idx] &&
      !anyDuplicated(indices) && all(indices >= 1L) &&
      all(indices <= metadata$pcm_D1_n[rep_idx] + metadata$pcm_D2_n[rep_idx])
  }, logical(1L))
  if (!all(valid_indices)) {
    stop("Invalid stored PCM D2 indices in ", path, call. = FALSE)
  }
  allowed_stop_reasons <- c("rejection", "futility", "cap")
  allowed_backends <- c("pcmt_native", "ordered_perm_stats_fallback")
  for (prefix in c("gcm", "pcm")) {
    actual_B <- metadata[[paste0(prefix, "_actual_B")]]
    stop_reason <- metadata[[paste0(prefix, "_stop_reason")]]
    loss_count <- metadata[[paste0(prefix, "_loss_count")]]
    t_star <- metadata[[paste0(prefix, "_t_star")]]
    stopped_early <- metadata[[paste0(prefix, "_stopped_early")]]
    backend <- metadata[[paste0(prefix, "_adaptive_backend")]]
    matched <- metadata[[paste0(prefix, "_matched_statistic")]]
    error <- metadata[[paste0(prefix, "_error")]]
    p_value <- rec$pvals[, paste0("perm_", prefix)]

    if (!is_integerish(actual_B) || any(actual_B > EXPECTED_B_MAX) ||
        !is_integerish(loss_count + 1) ||
        any(loss_count < 0L | loss_count > EXPECTED_H) ||
        !is.numeric(t_star) || anyNA(t_star) ||
        any(!(is.finite(t_star) | is.infinite(t_star))) ||
        !is.logical(stopped_early) || anyNA(stopped_early) ||
        !is.character(stop_reason) || anyNA(stop_reason) ||
        any(!stop_reason %in% allowed_stop_reasons) ||
        !is.character(backend) || anyNA(backend) ||
        any(!backend %in% allowed_backends) ||
        !is.logical(matched) || anyNA(matched) || !all(matched) ||
        !is.character(error) ||
        any(!is.na(error) & nzchar(error))) {
      stop("Invalid adaptive ", toupper(prefix), " metadata in ", path,
           call. = FALSE)
    }

    expected_stopped_early <- actual_B < EXPECTED_B_MAX
    valid_t_star <- ifelse(
      loss_count < EXPECTED_H,
      is.infinite(t_star) & t_star > 0,
      is.finite(t_star) & t_star == actual_B
    )
    expected_p <- ifelse(
      loss_count < EXPECTED_H,
      EXPECTED_H / (actual_B + EXPECTED_H - loss_count),
      EXPECTED_H / t_star
    )
    valid_reason <-
      (stop_reason == "rejection" & p_value <= EXPECTED_ALPHA) |
      (stop_reason == "futility" & p_value > EXPECTED_ALPHA &
         loss_count == EXPECTED_H) |
      (stop_reason == "cap" & actual_B == EXPECTED_B_MAX &
         p_value > EXPECTED_ALPHA & loss_count < EXPECTED_H)
    if (any(stopped_early != expected_stopped_early) ||
        !all(valid_t_star) || !all(valid_reason) ||
        !isTRUE(all.equal(
          unname(p_value), unname(expected_p), tolerance = 1e-12,
          check.attributes = FALSE
        ))) {
      stop("Adaptive ", toupper(prefix),
           " stopping metadata disagrees with its anytime-valid p-value in ",
           path, call. = FALSE)
    }
  }

  definition <- rec$method_definition
  required_definition <- c(
    "alternative", "regression", "B_requested", "B_max", "h", "alpha",
    "permutation_mode", "n_cores", "pcm_frac", "rf_num_trees",
    "rf_num_threads", "role_specific_rf_seeds", "gcm", "pcm",
    "seed_schedule"
  )
  if (!is.list(definition) ||
      !setequal(names(definition), required_definition) ||
      !identical(definition$alternative, "two.sided") ||
      !identical(definition$regression, "rf") ||
      !identical(as.integer(definition$B_requested), EXPECTED_B_MAX) ||
      !identical(as.integer(definition$B_max), EXPECTED_B_MAX) ||
      !identical(as.integer(definition$h), EXPECTED_H) ||
      !isTRUE(all.equal(
        as.numeric(definition$alpha), EXPECTED_ALPHA, tolerance = 0
      )) ||
      !identical(definition$permutation_mode, "adaptive_anytime_valid") ||
      !identical(as.integer(definition$n_cores), 1L) ||
      !identical(as.integer(definition$rf_num_threads), 1L) ||
      !isTRUE(definition$role_specific_rf_seeds) ||
      length(definition$rf_num_trees) != 1L ||
      !identical(as.integer(definition$rf_num_trees), EXPECTED_RF_NUM_TREES) ||
      length(definition$pcm_frac) != 1L ||
      !isTRUE(all.equal(
        as.numeric(definition$pcm_frac), EXPECTED_PCM_FRAC, tolerance = 0
      )) ||
      !identical(
        definition$gcm,
        paste(
          "One pcmt::perm_gcm fit; fit$statistic is normal-calibrated and",
          "fit$p.value is permutation-calibrated."
        )
      ) ||
      !identical(
        definition$pcm,
        paste(
          "One pcmt::perm_pcm fit with explicit D2 indices; fit$statistic",
          "is normal-calibrated and fit$p.value is permutation-calibrated."
        )
      )) {
    stop("Invalid v4 method definition in ", path, call. = FALSE)
  }
  seed_schedule <- definition$seed_schedule
  required_schedule_fields <- c(
    "base", "formula", "schedule_index", "setup_order",
    "paired_across_effect_sizes", "component_offsets"
  )
  if (!is.list(seed_schedule) ||
      !setequal(names(seed_schedule), required_schedule_fields) ||
      !identical(as.integer(seed_schedule$base), RUN_SEED_BASE) ||
      !identical(
        seed_schedule$formula,
        paste0(
          "base + paired_seed_idx * 100000 + rep * 100 + ",
          "component_offset"
        )
      ) ||
      !identical(
        seed_schedule$schedule_index,
        "paired_seed_idx = target_idx + 30"
      ) ||
      !identical(seed_schedule$setup_order, SETUP_ORDER) ||
      !isTRUE(seed_schedule$paired_across_effect_sizes) ||
      !identical(
        seed_schedule$component_offsets,
        EXPECTED_SEED_COMPONENT_OFFSETS
      )) {
    stop("Invalid v4 paired seed-schedule definition in ", path,
         call. = FALSE)
  }
  versions <- rec$package_versions
  if (!is.character(versions) || is.null(names(versions)) ||
      !setequal(names(versions), c("R", "pcmt", "comets", "ranger")) ||
      anyNA(versions) || any(!nzchar(versions))) {
    stop("Missing package-version provenance in ", path, call. = FALSE)
  }
  if (!identical(names(rec$rejection_rate), EXPECTED_METHODS) ||
      !identical(names(rec$na_rate), EXPECTED_METHODS) ||
      any(!is.finite(rec$rejection_rate)) ||
      any(!is.finite(rec$na_rate)) || any(rec$na_rate != 0)) {
    stop("Invalid rejection or missingness summary in ", path, call. = FALSE)
  }
  recomputed_rejection <- colMeans(rec$pvals <= 0.05)
  if (!isTRUE(all.equal(
    unname(rec$rejection_rate), unname(recomputed_rejection),
    tolerance = 1e-12
  ))) {
    stop("Stored rejection rates disagree with p-values in ", path,
         call. = FALSE)
  }
  invisible(TRUE)
}

if (status_mode) {
  statuses <- lapply(status_indices, function(i) {
    path <- expected_files[[i]]
    if (!file.exists(path)) {
      return(list(
        state = "missing",
        reason = "expected output is absent",
        package_versions = NULL,
        method_definition = NULL
      ))
    }
    failure <- tryCatch({
      obj <- readRDS(path)
      if (!is.list(obj) || length(obj) != 1L || !is.list(obj[[1L]])) {
        stop("Malformed per-setup result: ", path, call. = FALSE)
      }
      rec <- obj[[1L]]
      validate_record(rec, i, path)
      NULL
    }, error = function(e) conditionMessage(e))
    if (is.null(failure)) {
      list(
        state = "valid",
        reason = "passed the strict per-record v4 validator",
        package_versions = obj[[1L]]$package_versions,
        method_definition = obj[[1L]]$method_definition
      )
    } else {
      list(
        state = "invalid",
        reason = gsub("[\t\r\n]+", " ", failure),
        package_versions = NULL,
        method_definition = NULL
      )
    }
  })

  valid_indices <- which(vapply(
    statuses, function(status) identical(status$state, "valid"), logical(1L)
  ))
  heterogeneous_versions <- FALSE
  heterogeneous_definitions <- FALSE
  if (length(valid_indices) > 1L) {
    reference_versions <- statuses[[valid_indices[[1L]]]]$package_versions
    reference_definition <- statuses[[valid_indices[[1L]]]]$method_definition
    heterogeneous_versions <- !all(vapply(
      valid_indices,
      function(i) identical(statuses[[i]]$package_versions, reference_versions),
      logical(1L)
    ))
    heterogeneous_definitions <- !all(vapply(
      valid_indices,
      function(i) identical(
        statuses[[i]]$method_definition, reference_definition
      ),
      logical(1L)
    ))
  }
  if (heterogeneous_versions || heterogeneous_definitions) {
    heterogeneous_fields <- c(
      if (heterogeneous_versions) "package_versions",
      if (heterogeneous_definitions) "method_definition"
    )
    heterogeneity_reason <- paste0(
      "selected individually valid records are heterogeneous in ",
      paste(heterogeneous_fields, collapse = " and "),
      "; all selected valid records are invalidated for a coherent rerun"
    )
    for (i in valid_indices) {
      statuses[[i]]$state <- "invalid"
      statuses[[i]]$reason <- heterogeneity_reason
    }
  }

  ## A partial resume must not mix newly generated records with otherwise-valid
  ## records from a different runtime. Fully complete, internally homogeneous
  ## results remain eligible for the ordinary no-argument combiner.
  if (any(vapply(
        statuses, function(status) !identical(status$state, "valid"),
        logical(1L)
      ))) {
    current_versions <- c(
      R = R.version.string,
      pcmt = as.character(utils::packageVersion("pcmt")),
      comets = as.character(utils::packageVersion("comets")),
      ranger = as.character(utils::packageVersion("ranger"))
    )
    for (i in seq_along(statuses)) {
      if (identical(statuses[[i]]$state, "valid") &&
          !identical(statuses[[i]]$package_versions, current_versions)) {
        statuses[[i]]$state <- "invalid"
        statuses[[i]]$reason <- paste0(
          "package versions differ from the current runtime while selected ",
          "outputs are incomplete (stored: ",
          paste(
            paste0(
              names(statuses[[i]]$package_versions), "=",
              statuses[[i]]$package_versions
            ),
            collapse = ", "
          ),
          "; current: ",
          paste(
            paste0(names(current_versions), "=", current_versions),
            collapse = ", "
          ),
          ")"
        )
      }
    }
  }

  status_rows <- do.call(rbind, lapply(seq_along(status_indices), function(j) {
    i <- status_indices[[j]]
    data.frame(
      setup_idx = i,
      state = statuses[[j]]$state,
      filename = basename(expected_files[[i]]),
      reason = statuses[[j]]$reason,
      stringsAsFactors = FALSE
    )
  }))
  write.table(
    status_rows,
    file = status_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na = ""
  )
  quit(save = "no", status = 0L)
}

for (i in seq_along(expected_files)) {
  obj <- readRDS(expected_files[i])
  if (!is.list(obj) || length(obj) != 1L || !is.list(obj[[1L]])) {
    stop("Malformed per-setup result: ", expected_files[i], call. = FALSE)
  }

  ## Each per-setup file is a length-1 list containing the result record.
  rec <- obj[[1L]]
  validate_record(rec, i, expected_files[i])

  key <- sprintf("%s_%s_s%d", rec$target, rec$effect_id, rec$s)
  if (key %in% names(results)) {
    stop("Duplicate setup key: ", key, call. = FALSE)
  }
  results[[key]] <- rec
}

if (length(results) != 42L) {
  stop("Expected exactly 42 validated v4 setups; got ", length(results),
       call. = FALSE)
}

setup_ids <- vapply(results, function(rec) as.integer(rec$setup_idx), integer(1L))
target_ids <- vapply(results, function(rec) as.integer(rec$target_idx), integer(1L))
grid_ids <- vapply(results, function(rec) as.integer(rec$grid_idx), integer(1L))
effect_ids <- vapply(results, function(rec) as.character(rec$effect_id), character(1L))
effect_sizes <- vapply(results, function(rec) as.numeric(rec$effect_size), numeric(1L))
target_names <- vapply(results, function(rec) as.character(rec$target), character(1L))
expected_target_counts <- ifelse(
  TARGETS_DF$target %in% REPRESENTATIVE_TARGETS$target, 5L, 1L
)
actual_target_counts <- as.integer(table(factor(
  target_names,
  levels = TARGETS_DF$target
)))
if (!identical(unname(setup_ids), as.integer(SETUPS$setup_idx)) ||
    !identical(unname(target_ids), as.integer(SETUPS$target_idx)) ||
    !identical(unname(grid_ids), as.integer(SETUPS$grid_idx)) ||
    !identical(unname(effect_ids), as.character(SETUPS$effect_id)) ||
    !identical(unname(effect_sizes), as.numeric(SETUPS$effect_size)) ||
    !identical(unname(target_names), as.character(SETUPS$target)) ||
    !identical(actual_target_counts, expected_target_counts)) {
  stop(
    paste0(
      "Validated records do not form the exact 30-null plus ",
      "12-representative-alternative design."
    ),
    call. = FALSE
  )
}

reference_definition <- results[[1L]]$method_definition
reference_versions <- results[[1L]]$package_versions
if (!all(vapply(
      results,
      function(rec) identical(rec$method_definition, reference_definition),
      logical(1L)
    ))) {
  stop("The 42 v4 setups do not share one method definition.", call. = FALSE)
}
if (!all(vapply(
      results,
      function(rec) identical(rec$package_versions, reference_versions),
      logical(1L)
    ))) {
  stop("The 42 v4 setups were produced by different package versions.",
       call. = FALSE)
}

# Every nonzero representative effect must reuse exactly its target's delta0
# stochastic scaffolding. Method outputs and adaptive stopping diagnostics may
# differ with the outcome and are deliberately not compared. Across targets,
# every target/replicate/component stream must remain distinct. Validate
# pairing before checking global uniqueness on the 30 canonical delta0 records.
for (representative_idx in seq_len(nrow(REPRESENTATIVE_TARGETS))) {
  target <- REPRESENTATIVE_TARGETS$target[representative_idx]
  delta0_key <- sprintf("%s_delta0_s%d", target, FIXED_S)
  delta0_rec <- results[[delta0_key]]
  if (is.null(delta0_rec)) {
    stop("Missing delta0 pairing reference for target ", target,
         call. = FALSE)
  }
  delta0_metadata <- delta0_rec$replication_metadata

  target_effect_ids <- REPRESENTATIVE_EFFECT_GRID$effect_id[
    REPRESENTATIVE_EFFECT_GRID$target == target &
      REPRESENTATIVE_EFFECT_GRID$grid_idx > 1L
  ]
  if (length(target_effect_ids) != 4L) {
    stop("Representative target does not have four alternatives: ", target,
         call. = FALSE)
  }
  for (effect_id in target_effect_ids) {
    effect_key <- sprintf("%s_%s_s%d", target, effect_id, FIXED_S)
    effect_rec <- results[[effect_key]]
    if (is.null(effect_rec)) {
      stop("Missing effect ", effect_id, " for target ", target,
           call. = FALSE)
    }
    effect_metadata <- effect_rec$replication_metadata
    paired_seed_fields_match <- all(vapply(seed_fields, function(field) {
      identical(
        as.integer(delta0_metadata[[field]]),
        as.integer(effect_metadata[[field]])
      )
    }, logical(1L)))
    paired_scaffolding_matches <-
      identical(
        as.numeric(delta0_metadata$outcome_noise_sd),
        as.numeric(effect_metadata$outcome_noise_sd)
      ) &&
      identical(
        as.integer(delta0_metadata$pcm_D1_n),
        as.integer(effect_metadata$pcm_D1_n)
      ) &&
      identical(
        as.integer(delta0_metadata$pcm_D2_n),
        as.integer(effect_metadata$pcm_D2_n)
      ) &&
      identical(
        delta0_metadata$pcm_D2_indices,
        effect_metadata$pcm_D2_indices
      )
    if (!paired_seed_fields_match || !paired_scaffolding_matches) {
      stop(
        "Delta0/effect stochastic pairing mismatch for target ", target,
        " at ", effect_id,
        call. = FALSE
      )
    }
  }
}

delta0_results <- results[vapply(
  results,
  function(rec) identical(rec$effect_id, "delta0"),
  logical(1L)
)]
if (length(delta0_results) != nrow(TARGETS_DF)) {
  stop("Expected one canonical delta0 record for each target.", call. = FALSE)
}
canonical_component_seeds <- unlist(lapply(delta0_results, function(rec) {
  unlist(rec$replication_metadata[seed_fields], use.names = FALSE)
}), use.names = FALSE)
expected_seed_count <-
  nrow(TARGETS_DF) * EXPECTED_R * length(seed_fields)
if (length(canonical_component_seeds) != expected_seed_count ||
    anyDuplicated(canonical_component_seeds)) {
  stop(
    paste0(
      "Seed streams overlap across distinct target/replicate/component ",
      "combinations."
    ),
    call. = FALSE
  )
}

# ---- Long-format rejection table ----
rej_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(setup_idx = r$setup_idx,
             target_idx = r$target_idx,
             paired_seed_idx = r$paired_seed_idx,
             target = r$target,
             sparsity_tag = r$sparsity_tag,
             freq = r$freq,
             s = r$s,
             grid_idx = r$grid_idx,
             effect_id = r$effect_id,
             effect_size = r$effect_size,
             hypothesis = r$hypothesis,
             effect_scale = r$effect_scale,
             method = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             na_rate = unname(r$na_rate[names(r$rejection_rate)]),
             stringsAsFactors = FALSE)
}))
rownames(rej_df) <- NULL
rej_df$is_representative <-
  rej_df$target %in% REPRESENTATIVE_TARGETS$target
rej_df$representative_order <- match(
  rej_df$target,
  REPRESENTATIVE_TARGETS$target
)

delta0_rates <- rej_df[
  rej_df$effect_id == "delta0",
  c("target", "method", "rejection_rate")
]
names(delta0_rates)[names(delta0_rates) == "rejection_rate"] <-
  "delta0_rejection_rate"
gain_key <- paste(rej_df$target, rej_df$method, sep = "\r")
delta0_key <- paste(delta0_rates$target, delta0_rates$method, sep = "\r")
if (nrow(delta0_rates) != nrow(TARGETS_DF) * length(EXPECTED_METHODS) ||
    anyDuplicated(delta0_key)) {
  stop("Delta0 reference table is not one row per target and method.",
       call. = FALSE)
}
delta0_match <- match(gain_key, delta0_key)
if (anyNA(delta0_match)) {
  stop("Could not match every target/method row to its delta0 reference.",
       call. = FALSE)
}
rej_df$delta0_rejection_rate <-
  delta0_rates$delta0_rejection_rate[delta0_match]
rej_df$paired_gain_vs_delta0 <-
  rej_df$rejection_rate - rej_df$delta0_rejection_rate
if (any(rej_df$effect_id == "delta0" &
        rej_df$paired_gain_vs_delta0 != 0)) {
  stop("Delta0 rows must have zero paired rejection gain.", call. = FALSE)
}

# The power summary contains one exact cell per representative, grid point, and
# method. No averaging across targets is performed because each sparsity level
# deliberately has one representative power curve.
representative_summary <- rej_df[
  rej_df$is_representative,
  c(
    "representative_order", "target_idx", "target", "sparsity_tag", "freq",
    "s", "grid_idx", "effect_id", "effect_size", "hypothesis", "method",
    "rejection_rate", "delta0_rejection_rate", "paired_gain_vs_delta0",
    "na_rate"
  )
]
representative_summary$n_reps <- EXPECTED_R
representative_summary <- representative_summary[
  order(
    representative_summary$representative_order,
    representative_summary$grid_idx,
    match(representative_summary$method, EXPECTED_METHODS)
  ),
]
rownames(representative_summary) <- NULL
representative_cell_key <- paste(
  representative_summary$target,
  representative_summary$grid_idx,
  representative_summary$method,
  sep = "\r"
)
if (nrow(representative_summary) != 75L ||
    anyDuplicated(representative_cell_key) ||
    any(table(representative_summary$target) != 25L) ||
    any(table(representative_summary$grid_idx) != 15L) ||
    any(table(representative_summary$method) != 15L) ||
    any(representative_summary$n_reps != EXPECTED_R)) {
  stop(
    "Representative summary is not the exact 3-target by 5-grid by ",
    "5-method design.",
    call. = FALSE
  )
}

paired_gain_df <- representative_summary[
  representative_summary$hypothesis == "alternative",
  c(
    "representative_order", "target_idx", "target", "sparsity_tag",
    "grid_idx", "effect_id", "effect_size", "method", "rejection_rate",
    "delta0_rejection_rate", "paired_gain_vs_delta0", "n_reps"
  )
]
rownames(paired_gain_df) <- NULL
if (nrow(paired_gain_df) != 60L) {
  stop("Expected exactly 60 nonnull representative paired-gain cells.",
       call. = FALSE)
}

# Null calibration is still summarized over all 10 targets in each canonical
# group; this is separate from the representative-only power curves.
null_rej_df <- rej_df[rej_df$hypothesis == "null", ]
null_summary_by_sparsity <- aggregate(
  rejection_rate ~ sparsity_tag + method,
  data = null_rej_df,
  FUN = mean
)
names(null_summary_by_sparsity)[
  names(null_summary_by_sparsity) == "rejection_rate"
] <- "mean_rejection_rate"
sparsity_order <- c("low_sparsity", "medium_sparsity", "high_sparsity")
null_summary_by_sparsity$n_targets <- N_PER_BIN
null_summary_by_sparsity$n_reps_per_target <- EXPECTED_R
null_summary_by_sparsity <- null_summary_by_sparsity[
  order(
    match(null_summary_by_sparsity$sparsity_tag, sparsity_order),
    match(null_summary_by_sparsity$method, EXPECTED_METHODS)
  ),
  c(
    "sparsity_tag", "method", "mean_rejection_rate", "n_targets",
    "n_reps_per_target"
  )
]
rownames(null_summary_by_sparsity) <- NULL
null_cell_counts <- aggregate(
  target ~ sparsity_tag + method,
  data = null_rej_df,
  FUN = function(x) length(unique(x))
)
if (nrow(null_summary_by_sparsity) !=
      length(sparsity_order) * length(EXPECTED_METHODS) ||
    nrow(null_cell_counts) != nrow(null_summary_by_sparsity) ||
    any(null_cell_counts$target != N_PER_BIN) ||
    any(null_summary_by_sparsity$n_targets != N_PER_BIN) ||
    any(null_summary_by_sparsity$n_reps_per_target != EXPECTED_R)) {
  stop("Null sparsity-by-method summary cells must contain 10 targets.",
       call. = FALSE)
}

null_summary_by_method <- aggregate(
  rejection_rate ~ method,
  data = null_rej_df,
  FUN = mean
)
names(null_summary_by_method)[
  names(null_summary_by_method) == "rejection_rate"
] <- "mean_rejection_rate"
null_summary_by_method$n_targets <- nrow(TARGETS_DF)
null_summary_by_method$n_reps_per_target <- EXPECTED_R
null_summary_by_method <- null_summary_by_method[
  order(match(null_summary_by_method$method, EXPECTED_METHODS)),
]
rownames(null_summary_by_method) <- NULL

canonical_targets <- TARGETS_DF[, c("target_idx", "target", "tag", "freq")]
names(canonical_targets)[names(canonical_targets) == "tag"] <- "sparsity_tag"
setup_table <- SETUPS

combined <- list(
  schema_version = SCHEMA_VERSION,
  simulation_design = SIMULATION_DESIGN,
  calibration_design = CALIBRATION_DESIGN,
  matched_calibration = TRUE,
  n_targets = nrow(TARGETS_DF),
  n_representatives = nrow(REPRESENTATIVE_TARGETS),
  n_null_setups = nrow(TARGETS_DF),
  n_alternative_setups = nrow(SETUPS) - nrow(TARGETS_DF),
  n_setups = length(results),
  fixed_s = FIXED_S,
  canonical_targets = canonical_targets,
  representative_targets = REPRESENTATIVE_TARGETS,
  representative_effect_grid = REPRESENTATIVE_EFFECT_GRID,
  setup_table = setup_table,
  effect_scale = EFFECT_SCALE,
  setup_order = SETUP_ORDER,
  paired_seed_indices = seq_len(nrow(TARGETS_DF)) + nrow(TARGETS_DF),
  seed_pairing = list(
    formula = "paired_seed_idx = target_idx + 30",
    reference_effect_id = "delta0",
    paired_targets = REPRESENTATIVE_TARGETS$target,
    paired_across_effect_sizes = TRUE,
    component_fields = seed_fields
  ),
  R = EXPECTED_R,
  B_requested = EXPECTED_B_MAX,
  B_max = EXPECTED_B_MAX,
  adaptive_h = EXPECTED_H,
  adaptive_alpha = EXPECTED_ALPHA,
  alternative = "two.sided",
  methods = EXPECTED_METHODS,
  statistics = EXPECTED_STATISTICS,
  method_definition = reference_definition,
  package_versions = reference_versions,
  summary_by_representative_effect_method = representative_summary,
  paired_gain_vs_delta0_by_representative_effect_method = paired_gain_df,
  null_summary_by_sparsity_method = null_summary_by_sparsity,
  null_summary_by_method = null_summary_by_method,
  results = results
)
combined_path <- file.path(
  bpath,
  "ccle-typeI-synthY-30targets-matched-v4-results.rds"
)
rejection_path <- file.path(
  bpath,
  "ccle-typeI-synthY-30targets-matched-v4-rejection-rates.csv"
)
combined_tmp <- tempfile(
  pattern = paste0(basename(combined_path), ".tmp-"),
  tmpdir = dirname(combined_path)
)
rejection_tmp <- tempfile(
  pattern = paste0(basename(rejection_path), ".tmp-"),
  tmpdir = dirname(rejection_path)
)
cleanup_paths <- c(combined_tmp, rejection_tmp)
on.exit(unlink(cleanup_paths, recursive = TRUE, force = TRUE), add = TRUE)
saveRDS(combined, combined_tmp)
write.csv(rej_df, rejection_tmp, row.names = FALSE)

publish_pair <- function(rds_tmp, rds_path, csv_tmp, csv_path) {
  output_paths <- c(rds_path, csv_path)
  temp_paths <- c(rds_tmp, csv_tmp)
  existed <- file.exists(output_paths)
  backups <- rep(NA_character_, length(output_paths))

  for (i in seq_along(output_paths)) {
    if (existed[[i]]) {
      backups[[i]] <- tempfile(
        pattern = paste0(basename(output_paths[[i]]), ".backup-"),
        tmpdir = dirname(output_paths[[i]])
      )
      if (!file.copy(output_paths[[i]], backups[[i]], overwrite = FALSE)) {
        stop("Could not back up existing output before publication: ",
             output_paths[[i]], call. = FALSE)
      }
    }
  }
  backup_paths <- backups[!is.na(backups)]
  on.exit(unlink(backup_paths, recursive = TRUE, force = TRUE), add = TRUE)

  restore_prior_pair <- function() {
    restored <- logical(length(output_paths))
    for (i in seq_along(output_paths)) {
      if (existed[[i]]) {
        restored[[i]] <- file.copy(
          backups[[i]], output_paths[[i]], overwrite = TRUE
        )
      } else {
        unlink(output_paths[[i]], recursive = TRUE, force = TRUE)
        restored[[i]] <- !file.exists(output_paths[[i]])
      }
    }
    all(restored)
  }

  if (!file.rename(temp_paths[[1L]], output_paths[[1L]])) {
    stop("Failed to publish combined v4 RDS: ", output_paths[[1L]],
         call. = FALSE)
  }

  # This failure injection exists only to exercise rollback in isolated tests.
  if (identical(
        Sys.getenv("CCLE_COMBINER_TEST_FAIL_SECOND_PUBLISH", unset = "0"),
        "1"
      )) {
    unlink(temp_paths[[2L]])
  }

  if (!file.rename(temp_paths[[2L]], output_paths[[2L]])) {
    rollback_ok <- restore_prior_pair()
    stop(
      "Failed to publish v4 rejection CSV after publishing the RDS; ",
      if (rollback_ok) {
        "the prior output pair was restored."
      } else {
        "rollback was incomplete; inspect both output paths before reuse."
      },
      call. = FALSE
    )
  }
  invisible(TRUE)
}

publish_pair(combined_tmp, combined_path, rejection_tmp, rejection_path)

cat("=== Combined paired matched-calibration v4 summary (", length(results),
    " setups) ===\n", sep = "")
cat("Representative rejection rates and paired gains by target/grid/method:\n")
print(representative_summary)
cat("\nAll-target delta0 mean rejection rates by sparsity group/method:\n")
print(null_summary_by_sparsity)
cat("\nSaved:\n")
cat("  ", combined_path, "\n")
cat("  ", rejection_path, "\n")
