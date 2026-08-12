#!/usr/bin/env Rscript
###############################################################################
# Merge matched-calibration v2 results from ccle-typeI-synthY-30targets.R (one
# job per (target, s) pair) into a single versioned combined file.
#
# Inputs and outputs live in $LOCAL_PDRT_DATA_DIR/ccle:
#   ccle-typeI-synthY-30targets-matched-v2-results-{TARGET}-s{S}.rds
#   ccle-typeI-synthY-30targets-matched-v2-results.rds
#   ccle-typeI-synthY-30targets-matched-v2-rejection-rates.csv
#
# To know which setups are expected, this script re-derives the 60-element
# grid the same way the runner does: it samples 10 targets from each frequency
# tercile of the screened mutations using SAMPLE_SEED = 2026.
#
# Usage: Rscript combine-ccle-results.R
###############################################################################

S_VALUES    <- c(5L, 10L)
N_PER_BIN   <- 10L
SAMPLE_SEED <- 2026L
EXPECTED_R  <- 500L
EXPECTED_B_MAX <- 499L
EXPECTED_H <- 15L
EXPECTED_ALPHA <- 0.05
SCHEMA_VERSION     <- "ccle_synthy_matched_v2"
CALIBRATION_DESIGN <- "shared_observed_statistic"
EXPECTED_METHODS <- c(
  "asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive"
)
EXPECTED_STATISTICS <- c("gcm", "pcm")

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

SETUPS <- rbind(
  cbind(TARGETS_DF, s = S_VALUES[1L]),
  cbind(TARGETS_DF, s = S_VALUES[2L])
)
rownames(SETUPS) <- NULL

# ---- Load each per-setup file ----
expected_files <- file.path(
  bpath,
  sprintf(
    "ccle-typeI-synthY-30targets-matched-v2-results-%s-s%d.rds",
    SETUPS$target,
    SETUPS$s
  )
)
missing <- expected_files[!file.exists(expected_files)]
if (length(missing)) {
  legacy_files <- file.path(
    bpath,
    sprintf(
      "ccle-typeI-synthY-30targets-results-%s-s%d.rds",
      SETUPS$target,
      SETUPS$s
    )
  )
  n_legacy <- sum(file.exists(legacy_files[!file.exists(expected_files)]))
  stop(
    "Missing matched-calibration v2 per-setup files (", length(missing),
    " of ", length(expected_files), "):\n  ",
    paste(head(missing, 20), collapse = "\n  "),
    if (length(missing) > 20) "\n  ..." else "",
    if (n_legacy > 0L) paste0(
      "\nFound ", n_legacy, " corresponding legacy file(s); legacy results ",
      "are intentionally rejected. Rerun the v2 analysis."
    ) else "",
    call. = FALSE
  )
}

results <- list()

required_record_fields <- c(
  "schema_version", "calibration_design", "matched_calibration",
  "setup_idx", "target", "sparsity_tag", "freq", "s", "R",
  "B_requested", "B_max", "adaptive_h", "adaptive_alpha", "alternative",
  "pvals", "statistics",
  "replication_metadata", "method_definition", "package_versions",
  "rejection_rate", "na_rate"
)
required_metadata_fields <- c(
  "rep", "data_seed", "pcm_split_seed", "gcm_permutation_seed",
  "pcm_permutation_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "pcm_YonXZ_seed", "pcm_YonZ_seed", "pcm_YhatonZ_seed",
  "pcm_RonZ_seed", "pcm_D1_n", "pcm_D2_n", "gcm_actual_B",
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

is_integerish <- function(x) {
  is.numeric(x) && all(is.finite(x)) &&
    all(x == floor(x)) && all(x >= 1) && all(x <= .Machine$integer.max)
}

validate_record <- function(rec, setup_idx, path) {
  missing_fields <- setdiff(required_record_fields, names(rec))
  if (!is.list(rec) || length(missing_fields)) {
    stop(
      "Malformed matched-calibration v2 record in ", path,
      if (length(missing_fields)) paste0(
        "; missing field(s): ", paste(missing_fields, collapse = ", ")
      ) else "",
      call. = FALSE
    )
  }
  if (!identical(rec$schema_version, SCHEMA_VERSION) ||
      !identical(rec$calibration_design, CALIBRATION_DESIGN) ||
      !isTRUE(rec$matched_calibration)) {
    stop(
      "Legacy or unmatched calibration record rejected: ", path,
      call. = FALSE
    )
  }
  if (!identical(as.integer(rec$setup_idx), as.integer(setup_idx)) ||
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
    stop("Invalid or missing v2 p-value matrix in ", path, call. = FALSE)
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
      length(missing_metadata)) {
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
    "rf_num_threads", "role_specific_rf_seeds"
  )
  if (!is.list(definition) ||
      length(setdiff(required_definition, names(definition))) ||
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
      !is_integerish(definition$rf_num_trees) ||
      length(definition$pcm_frac) != 1L ||
      !is.finite(definition$pcm_frac) || definition$pcm_frac <= 0 ||
      definition$pcm_frac >= 1) {
    stop("Invalid v2 method definition in ", path, call. = FALSE)
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

  key <- sprintf("%s_s%d", rec$target, rec$s)
  if (key %in% names(results)) {
    stop("Duplicate setup key: ", key, call. = FALSE)
  }
  results[[key]] <- rec
}

if (length(results) != 60L) {
  stop("Expected exactly 60 validated v2 setups; got ", length(results),
       call. = FALSE)
}

reference_definition <- results[[1L]]$method_definition
reference_versions <- results[[1L]]$package_versions
if (!all(vapply(
      results,
      function(rec) identical(rec$method_definition, reference_definition),
      logical(1L)
    ))) {
  stop("The 60 v2 setups do not share one method definition.", call. = FALSE)
}
if (!all(vapply(
      results,
      function(rec) identical(rec$package_versions, reference_versions),
      logical(1L)
    ))) {
  stop("The 60 v2 setups were produced by different package versions.",
       call. = FALSE)
}
all_component_seeds <- unlist(lapply(results, function(rec) {
  unlist(rec$replication_metadata[seed_fields], use.names = FALSE)
}), use.names = FALSE)
if (anyDuplicated(all_component_seeds)) {
  stop("Seed streams overlap across the 60 v2 setups.", call. = FALSE)
}

combined <- list(
  schema_version = SCHEMA_VERSION,
  calibration_design = CALIBRATION_DESIGN,
  matched_calibration = TRUE,
  n_setups = length(results),
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
  results = results
)
combined_path <- file.path(
  bpath,
  "ccle-typeI-synthY-30targets-matched-v2-results.rds"
)
combined_tmp <- tempfile(
  pattern = paste0(basename(combined_path), ".tmp-"),
  tmpdir = dirname(combined_path)
)
on.exit(unlink(combined_tmp), add = TRUE)
saveRDS(combined, combined_tmp)
if (!file.rename(combined_tmp, combined_path)) {
  stop("Failed to atomically publish combined v2 result: ", combined_path,
       call. = FALSE)
}

# ---- Long-format rejection table ----
rej_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(target = r$target,
             sparsity_tag = r$sparsity_tag,
             freq = r$freq,
             s = r$s,
             method = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             na_rate = unname(r$na_rate[names(r$rejection_rate)]),
             stringsAsFactors = FALSE)
}))
rownames(rej_df) <- NULL
rejection_path <- file.path(
  bpath,
  "ccle-typeI-synthY-30targets-matched-v2-rejection-rates.csv"
)
write.csv(rej_df, rejection_path, row.names = FALSE)

cat("=== Combined matched-calibration v2 summary (", length(results),
    " setups) ===\n", sep = "")
cat(sprintf("Mean rejection rate per (method, s):\n"))
agg <- aggregate(rejection_rate ~ method + s, data = rej_df, FUN = mean)
print(agg)
cat("\nSaved:\n")
cat("  ", combined_path, "\n")
cat("  ", rejection_path, "\n")
