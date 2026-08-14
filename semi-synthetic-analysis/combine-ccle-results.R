#!/usr/bin/env Rscript
###############################################################################
# Merge the paired matched-calibration v4 effect grid from
# ccle-typeI-synthY-30targets.R (one job per target/effect setting) into a
# single versioned combined file.
#
# Inputs and outputs live in $LOCAL_PDRT_DATA_DIR/ccle:
#   ccle-typeI-synthY-30targets-matched-v4-results-{TARGET}-{EFFECT_ID}-s10.rds
#   ccle-typeI-synthY-30targets-matched-v4-results.rds
#   ccle-typeI-synthY-30targets-matched-v4-rejection-rates.csv
#
# To know which setups are expected, this script re-derives the same 30
# deterministic targets as the runner: it samples 10 targets from each
# frequency tercile of the screened mutations using SAMPLE_SEED = 2026. Each
# target has five paired effects (0, 0.5, 1, 1.5, 2) at s = 10. The setup IDs
# are effect-major: all 30 delta0 jobs, then all 30 delta0p5 jobs, and so on.
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
SIMULATION_DESIGN  <- "paired_effect_grid_s10"
CALIBRATION_DESIGN <- "shared_observed_statistic"
EFFECT_SCALE       <- "noise_sd_equals_sd_mu"
EFFECT_GRID <- data.frame(
  effect_idx = seq_len(5L),
  effect_id = c("delta0", "delta0p5", "delta1", "delta1p5", "delta2"),
  effect_size = c(0, 0.5, 1, 1.5, 2),
  hypothesis = c("null", rep("alternative", 4L)),
  stringsAsFactors = FALSE
)
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

SETUPS <- do.call(rbind, lapply(seq_len(nrow(EFFECT_GRID)), function(effect_idx) {
  effect <- EFFECT_GRID[effect_idx, ]
  transform(
    TARGETS_DF,
    s = FIXED_S,
    effect_idx = effect$effect_idx,
    effect_id = effect$effect_id,
    effect_size = effect$effect_size,
    hypothesis = effect$hypothesis,
    paired_seed_idx = target_idx + nrow(TARGETS_DF)
  )
}))
rownames(SETUPS) <- NULL
SETUPS$setup_idx <- seq_len(nrow(SETUPS))
stopifnot(
  nrow(SETUPS) == 150L,
  identical(as.integer(SETUPS$setup_idx), seq_len(150L)),
  identical(
    as.integer(SETUPS$effect_idx),
    rep(EFFECT_GRID$effect_idx, each = nrow(TARGETS_DF))
  ),
  identical(
    as.character(SETUPS$effect_id),
    rep(EFFECT_GRID$effect_id, each = nrow(TARGETS_DF))
  ),
  identical(
    as.integer(SETUPS$target_idx),
    rep(seq_len(nrow(TARGETS_DF)), times = nrow(EFFECT_GRID))
  ),
  identical(
    as.integer(SETUPS$paired_seed_idx),
    rep(seq_len(nrow(TARGETS_DF)) + nrow(TARGETS_DF),
        times = nrow(EFFECT_GRID))
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
if (length(missing)) {
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
  normalizePath(actual_v4_files, mustWork = TRUE),
  normalizePath(expected_files, mustWork = TRUE)
)
if (length(unexpected)) {
  stop(
    "Unexpected matched-v4 per-setup file(s); the directory must contain ",
    "exactly the canonical 150-file target-by-effect grid:\n  ",
    paste(head(unexpected, 20L), collapse = "\n  "),
    if (length(unexpected) > 20L) "\n  ..." else "",
    call. = FALSE
  )
}

results <- list()

required_record_fields <- c(
  "schema_version", "simulation_design", "calibration_design",
  "matched_calibration", "setup_idx", "target_idx", "paired_seed_idx",
  "target", "sparsity_tag", "freq", "s", "effect_idx", "effect_id",
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
        as.character(SETUPS$tag[setup_idx])
      ) ||
      !isTRUE(all.equal(
        as.numeric(rec$freq), as.numeric(SETUPS$freq[setup_idx]),
        tolerance = 1e-12
      ))) {
    stop("Setup metadata mismatch in ", path, call. = FALSE)
  }
  expected_effect_idx <- as.integer(SETUPS$effect_idx[setup_idx])
  expected_effect_id <- as.character(SETUPS$effect_id[setup_idx])
  expected_hypothesis <- SETUPS$hypothesis[setup_idx]
  expected_effect <- SETUPS$effect_size[setup_idx]
  if (!identical(as.integer(rec$effect_idx), expected_effect_idx) ||
      !identical(as.character(rec$effect_id), expected_effect_id) ||
      !identical(as.character(rec$hypothesis), expected_hypothesis) ||
      length(rec$effect_size) != 1L ||
      !is.numeric(rec$effect_size) || !is.finite(rec$effect_size) ||
      !isTRUE(all.equal(
        as.numeric(rec$effect_size), expected_effect, tolerance = 0
      )) ||
      !identical(rec$effect_scale, EFFECT_SCALE)) {
    stop("Effect-grid metadata mismatch in ", path, call. = FALSE)
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
      !identical(seed_schedule$setup_order, "effect_major_then_target") ||
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

if (length(results) != 150L) {
  stop("Expected exactly 150 validated v4 setups; got ", length(results),
       call. = FALSE)
}

setup_ids <- vapply(results, function(rec) as.integer(rec$setup_idx), integer(1L))
target_ids <- vapply(results, function(rec) as.integer(rec$target_idx), integer(1L))
effect_ids <- vapply(results, function(rec) as.character(rec$effect_id), character(1L))
target_names <- vapply(results, function(rec) as.character(rec$target), character(1L))
if (!identical(unname(sort(setup_ids)), seq_len(150L)) ||
    !identical(unname(sort(unique(target_ids))), seq_len(30L)) ||
    !identical(unname(unique(effect_ids)), EFFECT_GRID$effect_id) ||
    any(table(target_ids, factor(effect_ids, levels = EFFECT_GRID$effect_id)) != 1L) ||
    any(table(target_names, factor(effect_ids, levels = EFFECT_GRID$effect_id)) != 1L)) {
  stop(
    "Validated records do not form the exact 30-target by 5-effect cross-product.",
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
  stop("The 150 v4 setups do not share one method definition.", call. = FALSE)
}
if (!all(vapply(
      results,
      function(rec) identical(rec$package_versions, reference_versions),
      logical(1L)
    ))) {
  stop("The 150 v4 setups were produced by different package versions.",
       call. = FALSE)
}

# Every effect for a target must reuse exactly the delta0 stochastic
# scaffolding. Method outputs and adaptive stopping diagnostics may differ with
# the outcome and are deliberately not compared. Across targets, every
# target/replicate/component stream must remain distinct. Validate pairing
# before checking global uniqueness on the 30 canonical delta0 records only.
for (target_idx in seq_len(nrow(TARGETS_DF))) {
  target <- TARGETS_DF$target[target_idx]
  delta0_key <- sprintf("%s_delta0_s%d", target, FIXED_S)
  delta0_rec <- results[[delta0_key]]
  if (is.null(delta0_rec)) {
    stop("Missing delta0 pairing reference for target ", target,
         call. = FALSE)
  }
  delta0_metadata <- delta0_rec$replication_metadata

  for (effect_id in EFFECT_GRID$effect_id[-1L]) {
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
             effect_idx = r$effect_idx,
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

delta0_rates <- rej_df[
  rej_df$effect_id == "delta0",
  c("target", "method", "rejection_rate")
]
names(delta0_rates)[names(delta0_rates) == "rejection_rate"] <-
  "delta0_rejection_rate"
gain_key <- paste(rej_df$target, rej_df$method, sep = "\r")
delta0_key <- paste(delta0_rates$target, delta0_rates$method, sep = "\r")
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

summary_df <- aggregate(
  cbind(rejection_rate, paired_gain_vs_delta0) ~
    effect_idx + effect_id + effect_size + hypothesis + method,
  data = rej_df,
  FUN = mean
)
names(summary_df)[names(summary_df) == "rejection_rate"] <-
  "mean_rejection_rate"
names(summary_df)[names(summary_df) == "paired_gain_vs_delta0"] <-
  "mean_paired_gain_vs_delta0"
summary_df <- summary_df[
  order(
    summary_df$effect_idx,
    match(summary_df$method, EXPECTED_METHODS)
  ),
  c(
    "effect_idx", "effect_id", "effect_size", "hypothesis", "method",
    "mean_rejection_rate", "mean_paired_gain_vs_delta0"
  )
]
rownames(summary_df) <- NULL

paired_gain_df <- summary_df[
  summary_df$effect_id != "delta0",
  c(
    "effect_idx", "effect_id", "effect_size", "hypothesis", "method",
    "mean_paired_gain_vs_delta0"
  )
]
rownames(paired_gain_df) <- NULL

summary_by_sparsity <- aggregate(
  cbind(rejection_rate, paired_gain_vs_delta0) ~
    sparsity_tag + effect_idx + effect_id + effect_size + hypothesis + method,
  data = rej_df,
  FUN = mean
)
names(summary_by_sparsity)[
  names(summary_by_sparsity) == "rejection_rate"
] <- "mean_rejection_rate"
names(summary_by_sparsity)[
  names(summary_by_sparsity) == "paired_gain_vs_delta0"
] <- "mean_paired_gain_vs_delta0"
summary_by_sparsity$n_targets <- 10L
summary_by_sparsity$n_reps <- EXPECTED_R
sparsity_order <- c("low_sparsity", "medium_sparsity", "high_sparsity")
summary_by_sparsity <- summary_by_sparsity[
  order(
    match(summary_by_sparsity$sparsity_tag, sparsity_order),
    summary_by_sparsity$effect_idx,
    match(summary_by_sparsity$method, EXPECTED_METHODS)
  ),
  c(
    "sparsity_tag", "effect_idx", "effect_id", "effect_size", "hypothesis",
    "method", "mean_rejection_rate", "mean_paired_gain_vs_delta0",
    "n_targets", "n_reps"
  )
]
rownames(summary_by_sparsity) <- NULL

summary_cell_counts <- aggregate(
  target ~ sparsity_tag + effect_id + method,
  data = rej_df,
  FUN = function(x) length(unique(x))
)
if (nrow(summary_by_sparsity) !=
      length(sparsity_order) * nrow(EFFECT_GRID) * length(EXPECTED_METHODS) ||
    nrow(summary_cell_counts) != nrow(summary_by_sparsity) ||
    any(summary_cell_counts$target != 10L) ||
    any(summary_by_sparsity$n_targets != 10L) ||
    any(summary_by_sparsity$n_reps != EXPECTED_R)) {
  stop("Sparsity-by-effect summary cells must contain exactly 10 targets.",
       call. = FALSE)
}

combined <- list(
  schema_version = SCHEMA_VERSION,
  simulation_design = SIMULATION_DESIGN,
  calibration_design = CALIBRATION_DESIGN,
  matched_calibration = TRUE,
  n_targets = nrow(TARGETS_DF),
  n_setups = length(results),
  fixed_s = FIXED_S,
  effect_grid = EFFECT_GRID,
  effect_ids = EFFECT_GRID$effect_id,
  effect_sizes = stats::setNames(EFFECT_GRID$effect_size, EFFECT_GRID$effect_id),
  hypotheses = stats::setNames(EFFECT_GRID$hypothesis, EFFECT_GRID$effect_id),
  effect_scale = EFFECT_SCALE,
  setup_order = "effect_major_then_target",
  paired_seed_indices = seq_len(nrow(TARGETS_DF)) + nrow(TARGETS_DF),
  seed_pairing = list(
    formula = "paired_seed_idx = target_idx + 30",
    reference_effect_id = "delta0",
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
  summary_by_effect_method = summary_df,
  paired_gain_vs_delta0_by_effect_method = paired_gain_df,
  summary_by_sparsity_effect_method = summary_by_sparsity,
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
cat("Mean rejection rate and paired gain per (effect, method):\n")
print(summary_df)
cat("\nSaved:\n")
cat("  ", combined_path, "\n")
cat("  ", rejection_path, "\n")
