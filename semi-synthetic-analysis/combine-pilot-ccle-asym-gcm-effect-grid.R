#!/usr/bin/env Rscript
###############################################################################
# Strictly combine all 30 per-target asymptotic-GCM effect-grid pilot files.
# Configuration is read from the same environment variables as the runner:
# PILOT_R, PILOT_RF_TREES, PILOT_RF_THREADS, and PILOT_OUTPUT_DIR.
###############################################################################

SCHEMA_VERSION <- "ccle_asym_gcm_effect_pilot_v3"
COMBINED_SCHEMA_VERSION <- "ccle_asym_gcm_effect_pilot_combined_v3"
SIMULATION_DESIGN <- "paired_effect_grid_s10"
METHOD <- "asym_gcm_rf_cached_x"
FIXED_S <- 10L
EFFECT_SIZES <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)
SOURCE_EFFECT_GRIDS <- stats::setNames(
  list(
    c(0, 0.25, 0.5, 0.75, 1, 1.25),
    c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)
  ),
  c("ccle_asym_gcm_effect_pilot_v1", "ccle_asym_gcm_effect_pilot_v2")
)
RUN_SEED_BASE <- 30260000L
SAMPLE_SEED <- 2026L
N_PER_BIN <- 10L
ALPHA <- 0.05
SPARSITY_DISPLAY_LABELS <- c(
  low_sparsity = "High sparsity",
  medium_sparsity = "Medium sparsity",
  high_sparsity = "Low sparsity"
)
SEED_FORMULA <-
  "base + paired_seed_idx * 100000 + rep * 100 + component_offset"
METHOD_DEFINITION <- list(
  statistic = paste(
    "pcmt signed studentized mean of residual products, using",
    "mean(R^2) - mean(R)^2"
  ),
  normal_p_value = "2 * pnorm(abs(statistic), lower.tail = FALSE)",
  x_fit_scope = "one X|Z ranger fit per target and replicate, reused across deltas",
  y_fit_scope = "one Y|Z ranger fit per target, replicate, and delta",
  regression = "ranger random forest",
  probability_forest_for_binary_x = TRUE
)

read_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L ||
      !identical(as.character(parsed), value)) {
    stop(name, " must be one positive integer; got: ", value, call. = FALSE)
  }
  parsed
}

PILOT_R <- read_positive_integer_env("PILOT_R", 25L)
PILOT_RF_TREES <- read_positive_integer_env("PILOT_RF_TREES", 500L)
PILOT_RF_THREADS <- read_positive_integer_env("PILOT_RF_THREADS", 4L)

data_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(data_root)) {
  stop("LOCAL_PDRT_DATA_DIR must point to the data/results directory.",
       call. = FALSE)
}
data_root <- path.expand(data_root)
if (!dir.exists(data_root)) {
  stop("LOCAL_PDRT_DATA_DIR does not exist: ", data_root, call. = FALSE)
}
data_root <- normalizePath(data_root, mustWork = TRUE)
output_dir <- Sys.getenv(
  "PILOT_OUTPUT_DIR",
  unset = file.path(data_root, "ccle", "pilot-asym-gcm-grid")
)
output_dir <- path.expand(output_dir)
if (!dir.exists(output_dir)) {
  stop("PILOT_OUTPUT_DIR does not exist: ", output_dir, call. = FALSE)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# Re-derive the exact screened data and deterministic 30-target grid used by
# both the production runner and pilot. This prevents a complete-looking set
# of artifacts from a different CCLE preprocessing run from being combined.
ccle_dir <- file.path(data_root, "ccle")
input_files <- file.path(ccle_dir, c("mutation.txt", "response.csv"))
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing CCLE input(s):\n  ", paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}
mutations <- read.delim(file.path(ccle_dir, "mutation.txt"), sep = "\t")
mutations <- mutations[grep("_MUT$", mutations$Name), ]
response <- read.csv(file.path(ccle_dir, "response.csv"))
response <- response[response$Compound == "PLX4720", ]
idx <- match(response$CCLE.Cell.Line.Name, colnames(mutations))
resp <- response[!is.na(idx), ]
mut <- mutations[, c(1, idx[!is.na(idx)])]
rownames(mut) <- mut$Name
mut$Name <- NULL
Y_real <- resp$Amax
cors <- suppressWarnings(apply(mut, 1, function(m) {
  value <- cor(Y_real, m)
  ifelse(!is.na(value), value, 0)
}))
screen <- cors[abs(cors) > 0.05 & !is.na(cors)]
dat <- data.frame(Y_real = Y_real, t(mut[names(screen), ]))
mut_names <- setdiff(colnames(dat), "Y_real")
M_full <- as.matrix(dat[, mut_names])
freqs <- colMeans(M_full)
qs <- quantile(freqs, probs = c(1 / 3, 2 / 3), names = FALSE)
pool_low <- sort(names(freqs)[freqs <= qs[1L]])
pool_med <- sort(names(freqs)[freqs > qs[1L] & freqs <= qs[2L]])
pool_high <- sort(names(freqs)[freqs > qs[2L]])
if (any(c(length(pool_low), length(pool_med), length(pool_high)) < N_PER_BIN)) {
  stop("A frequency bin has fewer than 10 screened mutations.",
       call. = FALSE)
}
set.seed(SAMPLE_SEED)
targets_low <- sort(sample(pool_low, N_PER_BIN))
targets_med <- sort(sample(pool_med, N_PER_BIN))
targets_high <- sort(sample(pool_high, N_PER_BIN))
canonical_targets <- rbind(
  data.frame(
    target = targets_low,
    sparsity_tag = "low_sparsity", stringsAsFactors = FALSE
  ),
  data.frame(
    target = targets_med,
    sparsity_tag = "medium_sparsity", stringsAsFactors = FALSE
  ),
  data.frame(
    target = targets_high,
    sparsity_tag = "high_sparsity", stringsAsFactors = FALSE
  )
)
canonical_targets$target_idx <- seq_len(nrow(canonical_targets))
canonical_targets$freq <- unname(freqs[canonical_targets$target])
canonical_targets$carrier_count <- as.integer(colSums(
  M_full[, canonical_targets$target, drop = FALSE]
))
if (nrow(canonical_targets) != 30L) {
  stop("Canonical pilot target grid must contain exactly 30 targets.",
       call. = FALSE)
}

expected_paths <- vapply(seq_len(30L), function(target_idx) {
  pattern <- sprintf(
    "^ccle-asym-gcm-effect-grid-s10-pilot-target-%02d-.*-R%d-trees%d\\.rds$",
    target_idx, PILOT_R, PILOT_RF_TREES
  )
  matches <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(matches) != 1L) {
    stop(
      "Expected exactly one RDS for target ", target_idx, "; found ",
      length(matches), ". Pattern: ", pattern,
      call. = FALSE
    )
  }
  matches[[1L]]
}, character(1L))
expected_csv_paths <- sub("[.]rds$", ".csv", expected_paths)
if (any(!file.exists(expected_csv_paths))) {
  stop(
    "Missing per-target CSV companion(s):\n  ",
    paste(expected_csv_paths[!file.exists(expected_csv_paths)],
          collapse = "\n  "),
    call. = FALSE
  )
}

required_fields <- c(
  "schema_version", "simulation_design", "calibration", "method",
  "fixed_s", "effect_sizes", "alpha", "R", "rf_num_trees",
  "rf_num_threads", "seed_base", "seed_formula", "target_idx", "paired_seed_idx",
  "target", "sparsity_tag", "freq", "carrier_count", "n", "p",
  "method_definition", "preflight", "package_versions",
  "resume_provenance", "summary", "raw"
)
required_raw_fields <- c(
  "target_idx", "target", "sparsity_tag", "freq", "carrier_count", "rep",
  "paired_seed_idx", "data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "outcome_noise_sd", "effect_size", "statistic", "p_value", "rejected",
  "error", "elapsed_seconds"
)
resume_provenance_fields <- c(
  "mode", "source_schema_version", "retained_effect_sizes",
  "computed_effect_sizes", "source_rds_md5", "source_csv_md5"
)
preflight_fields <- c(
  "passed", "tolerance", "cached_statistic", "pcmt_statistic",
  "statistic_difference", "max_abs_rY_difference",
  "max_abs_rX_difference", "pcmt_B"
)

validate_preflight <- function(preflight, path) {
  if (!is.list(preflight) ||
      !identical(names(preflight), preflight_fields) ||
      !isTRUE(preflight$passed) ||
      !is.finite(preflight$tolerance) || preflight$tolerance <= 0 ||
      !all(is.finite(c(
        preflight$cached_statistic, preflight$pcmt_statistic,
        preflight$statistic_difference, preflight$max_abs_rY_difference,
        preflight$max_abs_rX_difference
      ))) ||
      !isTRUE(all.equal(
        preflight$statistic_difference,
        preflight$pcmt_statistic - preflight$cached_statistic,
        tolerance = 1e-15
      )) ||
      abs(preflight$statistic_difference) > preflight$tolerance ||
      preflight$max_abs_rY_difference < 0 ||
      preflight$max_abs_rY_difference > preflight$tolerance ||
      preflight$max_abs_rX_difference < 0 ||
      preflight$max_abs_rX_difference > preflight$tolerance ||
      !identical(as.integer(preflight$pcmt_B), 1L)) {
    stop("Invalid cached-ranger preflight: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

validate_resume_provenance <- function(provenance, path) {
  if (!is.list(provenance) ||
      !identical(names(provenance), resume_provenance_fields) ||
      length(provenance$mode) != 1L ||
      !provenance$mode %in% c(
        "fresh_eight_effect_grid", "upgraded_effect_grid_prefix"
      )) {
    stop("Invalid pilot resume provenance: ", path, call. = FALSE)
  }
  if (identical(provenance$mode, "fresh_eight_effect_grid")) {
    valid <- length(provenance$source_schema_version) == 1L &&
      is.na(provenance$source_schema_version) &&
      length(provenance$retained_effect_sizes) == 0L &&
      identical(as.numeric(provenance$computed_effect_sizes), EFFECT_SIZES) &&
      length(provenance$source_rds_md5) == 1L &&
      is.na(provenance$source_rds_md5) &&
      length(provenance$source_csv_md5) == 1L &&
      is.na(provenance$source_csv_md5)
  } else {
    md5_pattern <- "^[[:xdigit:]]{32}$"
    source_schema <- provenance$source_schema_version
    source_effects <- if (length(source_schema) == 1L &&
                          !is.na(source_schema)) {
      SOURCE_EFFECT_GRIDS[[source_schema]]
    } else {
      NULL
    }
    missing_effects <- if (is.null(source_effects)) NULL else {
      EFFECT_SIZES[(length(source_effects) + 1L):length(EFFECT_SIZES)]
    }
    valid <- !is.null(source_effects) && identical(
      as.numeric(provenance$retained_effect_sizes), source_effects
    ) && identical(
      as.numeric(provenance$computed_effect_sizes), missing_effects
    ) &&
      length(provenance$source_rds_md5) == 1L &&
      grepl(md5_pattern, provenance$source_rds_md5) &&
      length(provenance$source_csv_md5) == 1L &&
      grepl(md5_pattern, provenance$source_csv_md5)
  }
  if (!isTRUE(valid)) {
    stop("Invalid pilot resume provenance: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

validate_record <- function(record, target_idx, path) {
  canonical <- canonical_targets[target_idx, , drop = FALSE]
  if (!is.list(record) || !identical(names(record), required_fields)) {
    stop("Malformed pilot record: ", path, call. = FALSE)
  }
  if (!identical(record$schema_version, SCHEMA_VERSION) ||
      !identical(record$simulation_design, SIMULATION_DESIGN) ||
      !identical(record$calibration, "normal_reference") ||
      !identical(record$method, METHOD) ||
      !identical(as.integer(record$fixed_s), FIXED_S) ||
      !identical(as.numeric(record$effect_sizes), EFFECT_SIZES) ||
      !identical(as.integer(record$R), PILOT_R) ||
      !identical(as.integer(record$rf_num_trees), PILOT_RF_TREES) ||
      !identical(as.integer(record$rf_num_threads), PILOT_RF_THREADS) ||
      !identical(as.integer(record$seed_base), RUN_SEED_BASE) ||
      !identical(record$seed_formula, SEED_FORMULA) ||
      !identical(as.integer(record$target_idx), target_idx) ||
      !identical(as.integer(record$paired_seed_idx), target_idx + 30L) ||
      !identical(record$target, canonical$target[[1L]]) ||
      !identical(record$sparsity_tag, canonical$sparsity_tag[[1L]]) ||
      !isTRUE(all.equal(
        as.numeric(record$freq), canonical$freq[[1L]], tolerance = 0
      )) ||
      !identical(
        as.integer(record$carrier_count), canonical$carrier_count[[1L]]
      ) ||
      !identical(as.integer(record$n), nrow(M_full)) ||
      !identical(as.integer(record$p), ncol(M_full)) ||
      !isTRUE(all.equal(as.numeric(record$alpha), ALPHA, tolerance = 0)) ||
      !identical(record$method_definition, METHOD_DEFINITION) ||
      !is.character(record$package_versions) ||
      !identical(
        names(record$package_versions), c("R", "pcmt", "comets", "ranger")
      ) || anyNA(record$package_versions) ||
      any(!nzchar(record$package_versions))) {
    stop("Pilot configuration/provenance mismatch: ", path, call. = FALSE)
  }
  validate_preflight(record$preflight, path)
  validate_resume_provenance(record$resume_provenance, path)
  raw <- record$raw
  if (!is.data.frame(raw) ||
      !identical(names(raw), required_raw_fields) ||
      nrow(raw) != PILOT_R * length(EFFECT_SIZES) ||
      any(raw$target_idx != target_idx) ||
      any(raw$target != record$target) ||
      any(raw$sparsity_tag != record$sparsity_tag) ||
      !isTRUE(all.equal(
        as.numeric(raw$freq), rep(as.numeric(record$freq), nrow(raw)),
        tolerance = 0, check.attributes = FALSE
      )) ||
      any(as.integer(raw$carrier_count) != as.integer(record$carrier_count)) ||
      any(raw$paired_seed_idx != target_idx + 30L) ||
      !identical(sort(unique(raw$rep)), seq_len(PILOT_R)) ||
      !identical(sort(unique(raw$effect_size)), EFFECT_SIZES) ||
      any(table(raw$rep) != length(EFFECT_SIZES)) ||
      any(table(raw$effect_size) != PILOT_R) ||
      any(!is.na(raw$error)) || any(!is.finite(raw$statistic)) ||
      any(!is.finite(raw$p_value)) || any(raw$p_value < 0 | raw$p_value > 1) ||
      any(!is.finite(raw$outcome_noise_sd)) ||
      any(raw$outcome_noise_sd <= 0) ||
      any(!is.finite(raw$elapsed_seconds)) || any(raw$elapsed_seconds < 0) ||
      any(raw$rejected != (raw$p_value <= ALPHA)) ||
      !isTRUE(all.equal(
        raw$p_value,
        2 * stats::pnorm(abs(raw$statistic), lower.tail = FALSE),
        tolerance = 1e-15, check.attributes = FALSE
      ))) {
    stop("Invalid raw pilot rows: ", path, call. = FALSE)
  }
  complete_effect_grid <- all(vapply(split(raw, raw$rep), function(block) {
    identical(sort(block$effect_size), EFFECT_SIZES)
  }, logical(1L)))
  expected_data_seed <- as.integer(
    RUN_SEED_BASE + (target_idx + 30L) * 100000 + raw$rep * 100 + 1L
  )
  expected_y_seed <- expected_data_seed + 4L
  expected_x_seed <- expected_data_seed + 5L
  if (!complete_effect_grid ||
      !identical(as.integer(raw$data_seed), expected_data_seed) ||
      !identical(as.integer(raw$gcm_YonZ_seed), expected_y_seed) ||
      !identical(as.integer(raw$gcm_XonZ_seed), expected_x_seed)) {
    stop("Effect grid or exact seed schedule mismatch: ", path,
         call. = FALSE)
  }
  paired_fields <- c(
    "data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed", "outcome_noise_sd"
  )
  paired_ok <- all(vapply(split(raw, raw$rep), function(block) {
    all(vapply(paired_fields, function(field) {
      length(unique(block[[field]])) == 1L
    }, logical(1L)))
  }, logical(1L)))
  if (!paired_ok) {
    stop("Within-replicate effect pairing failed: ", path, call. = FALSE)
  }
  expected_summary <- aggregate(
    rejected ~ effect_size, data = raw, FUN = mean
  )
  names(expected_summary)[2L] <- "rejection_rate"
  expected_summary <- expected_summary[
    match(EFFECT_SIZES, expected_summary$effect_size), , drop = FALSE
  ]
  rownames(expected_summary) <- NULL
  expected_summary$mcse <- sqrt(
    expected_summary$rejection_rate *
      (1 - expected_summary$rejection_rate) / PILOT_R
  )
  expected_summary$n_success <- PILOT_R
  if (!is.data.frame(record$summary) ||
      !identical(names(record$summary), names(expected_summary)) ||
      !isTRUE(all.equal(
        record$summary, expected_summary, tolerance = 1e-15,
        check.attributes = FALSE
      ))) {
    stop("Per-target pilot summary mismatch: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

validate_csv_copy <- function(csv_path, raw, context) {
  csv <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop("Could not read pilot CSV in ", context, ": ",
           conditionMessage(e), call. = FALSE)
    }
  )
  comparable_fields <- setdiff(required_raw_fields, "error")
  if (!identical(names(csv), required_raw_fields) ||
      nrow(csv) != nrow(raw) ||
      any(!is.na(csv$error)) || any(!is.na(raw$error)) ||
      !isTRUE(all.equal(
        csv[comparable_fields], raw[comparable_fields], tolerance = 1e-12,
        check.attributes = FALSE
      ))) {
    stop("Pilot CSV does not reproduce the RDS raw rows in ", context,
         call. = FALSE)
  }
  invisible(TRUE)
}

records <- lapply(seq_along(expected_paths), function(target_idx) {
  record <- readRDS(expected_paths[[target_idx]])
  validate_record(record, target_idx, expected_paths[[target_idx]])
  validate_csv_copy(
    expected_csv_paths[[target_idx]], record$raw,
    paste0("target ", target_idx, " artifact")
  )
  record
})
names(records) <- sprintf("target_%02d", seq_len(30L))

reference_versions <- records[[1L]]$package_versions
if (!all(vapply(records, function(record) {
  identical(record$package_versions, reference_versions)
}, logical(1L)))) {
  stop("Pilot files were produced by different package versions.",
       call. = FALSE)
}
if (length(unique(vapply(records, `[[`, character(1L), "target"))) != 30L ||
    !identical(
      as.integer(vapply(records, `[[`, integer(1L), "target_idx")),
      seq_len(30L)
    ) ||
    any(table(vapply(records, `[[`, character(1L), "sparsity_tag")) != 10L)) {
  stop("The 30 pilot records do not contain the canonical target/bin grid.",
       call. = FALSE)
}

raw <- do.call(rbind, lapply(records, `[[`, "raw"))
rownames(raw) <- NULL
canonical_seeds <- unique(raw[
  raw$effect_size == EFFECT_SIZES[[1L]],
  c("target_idx", "rep", "data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed")
])
seed_values <- unlist(canonical_seeds[
  c("data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed")
], use.names = FALSE)
if (nrow(canonical_seeds) != 30L * PILOT_R || anyDuplicated(seed_values)) {
  stop("Pilot seed streams overlap across target/rep/component.", call. = FALSE)
}

target_summary <- aggregate(
  rejected ~ target_idx + target + sparsity_tag + freq + carrier_count +
    effect_size,
  data = raw,
  FUN = mean
)
names(target_summary)[names(target_summary) == "rejected"] <- "rejection_rate"
target_summary$mcse <- sqrt(
  target_summary$rejection_rate * (1 - target_summary$rejection_rate) / PILOT_R
)
target_summary$n_reps <- PILOT_R
target_summary$sparsity_label <- unname(
  SPARSITY_DISPLAY_LABELS[target_summary$sparsity_tag]
)
if (anyNA(target_summary$sparsity_label)) {
  stop("Unknown pilot sparsity tag while assigning display labels.",
       call. = FALSE)
}

null_rates <- target_summary[
  target_summary$effect_size == 0,
  c("target_idx", "rejection_rate")
]
names(null_rates)[2L] <- "null_rejection_rate"
target_summary <- merge(target_summary, null_rates, by = "target_idx",
                        all.x = TRUE, sort = FALSE)
target_summary$paired_gain <-
  target_summary$rejection_rate - target_summary$null_rejection_rate

paired_difference_rows <- do.call(rbind, lapply(split(raw, raw$target_idx),
                                                function(target_raw) {
  null_rows <- target_raw[target_raw$effect_size == 0,
                          c("rep", "rejected")]
  names(null_rows)[2L] <- "null_rejected"
  merged <- merge(target_raw, null_rows, by = "rep", all.x = TRUE,
                  sort = FALSE)
  merged$paired_difference <-
    as.numeric(merged$rejected) - as.numeric(merged$null_rejected)
  merged
}))
paired_mcse <- aggregate(
  paired_difference ~ target_idx + effect_size,
  data = paired_difference_rows,
  FUN = function(x) {
    if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
  }
)
names(paired_mcse)[3L] <- "paired_gain_mcse"
target_summary <- merge(
  target_summary,
  paired_mcse,
  by = c("target_idx", "effect_size"),
  all.x = TRUE,
  sort = FALSE
)

equal_target_summary <- aggregate(
  cbind(rejection_rate, paired_gain) ~ sparsity_tag + effect_size,
  data = target_summary,
  FUN = mean
)
names(equal_target_summary)[3:4] <- c(
  "equal_target_rejection_rate", "equal_target_paired_gain"
)
equal_target_summary$n_targets <- 10L
equal_target_summary$n_reps_per_target <- PILOT_R
equal_target_rate_mcse <- aggregate(
  rejection_rate ~ sparsity_tag + effect_size,
  data = target_summary,
  FUN = function(target_rates) {
    sqrt(sum(
      target_rates * (1 - target_rates) / PILOT_R
    )) / length(target_rates)
  }
)
names(equal_target_rate_mcse)[3L] <- "mcse"
equal_target_summary <- merge(
  equal_target_summary,
  equal_target_rate_mcse,
  by = c("sparsity_tag", "effect_size"),
  all.x = TRUE,
  sort = FALSE
)
equal_target_paired_mcse <- aggregate(
  paired_gain_mcse ~ sparsity_tag + effect_size,
  data = target_summary,
  FUN = function(x) {
    if (anyNA(x)) NA_real_ else sqrt(sum(x^2)) / length(x)
  },
  na.action = stats::na.pass
)
names(equal_target_paired_mcse)[3L] <- "paired_gain_mcse"
equal_target_summary <- merge(
  equal_target_summary,
  equal_target_paired_mcse,
  by = c("sparsity_tag", "effect_size"),
  all.x = TRUE,
  sort = FALSE
)
equal_target_summary <- equal_target_summary[
  order(
    match(equal_target_summary$sparsity_tag,
          c("low_sparsity", "medium_sparsity", "high_sparsity")),
    equal_target_summary$effect_size
  ),
]
equal_target_summary$sparsity_label <- unname(
  SPARSITY_DISPLAY_LABELS[equal_target_summary$sparsity_tag]
)
equal_target_summary <- equal_target_summary[, c(
  "sparsity_tag", "sparsity_label", "effect_size",
  "equal_target_rejection_rate", "equal_target_paired_gain", "n_targets",
  "n_reps_per_target", "mcse", "paired_gain_mcse"
)]
rownames(equal_target_summary) <- NULL

target_summary <- target_summary[, c(
  "target_idx", "target", "sparsity_tag", "sparsity_label", "freq",
  "carrier_count", "effect_size", "rejection_rate", "mcse", "n_reps",
  "null_rejection_rate", "paired_gain", "paired_gain_mcse"
)]

combined <- list(
  schema_version = COMBINED_SCHEMA_VERSION,
  source_schema_version = SCHEMA_VERSION,
  simulation_design = SIMULATION_DESIGN,
  calibration = "normal_reference",
  method = METHOD,
  fixed_s = FIXED_S,
  effect_sizes = EFFECT_SIZES,
  alpha = ALPHA,
  R = PILOT_R,
  rf_num_trees = PILOT_RF_TREES,
  rf_num_threads = PILOT_RF_THREADS,
  seed_base = RUN_SEED_BASE,
  n_targets = 30L,
  sparsity_display_labels = SPARSITY_DISPLAY_LABELS,
  package_versions = reference_versions,
  equal_target_summary = equal_target_summary,
  target_summary = target_summary,
  raw = raw,
  records = records
)

validate_summary_csv <- function(path, object, context) {
  csv <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop("Could not read staged ", context, ": ", conditionMessage(e),
           call. = FALSE)
    }
  )
  if (!identical(names(csv), names(object)) || nrow(csv) != nrow(object)) {
    stop("Staged ", context, " has the wrong schema.", call. = FALSE)
  }
  for (field in names(object)) {
    if (all(is.na(csv[[field]])) && all(is.na(object[[field]]))) {
      csv[[field]] <- object[[field]]
    }
  }
  if (!isTRUE(all.equal(
    csv, object, tolerance = 1e-12, check.attributes = FALSE
  ))) {
    stop("Staged ", context, " does not reproduce its source data.",
         call. = FALSE)
  }
  invisible(TRUE)
}

artifact_snapshot <- function(paths) {
  exists <- file.exists(paths)
  md5 <- rep(NA_character_, length(paths))
  if (any(exists)) md5[exists] <- unname(tools::md5sum(paths[exists]))
  names(exists) <- names(md5) <- names(paths)
  list(exists = exists, md5 = md5)
}

publish_combined_outputs <- function(combined, equal_summary, target_summary,
                                     paths) {
  snapshot <- artifact_snapshot(paths)
  tmp_paths <- vapply(names(paths), function(name) {
    tempfile(
      pattern = paste0(basename(paths[[name]]), ".staged-"),
      tmpdir = dirname(paths[[name]])
    )
  }, character(1L))
  backup_paths <- vapply(names(paths), function(name) {
    tempfile(
      pattern = paste0(basename(paths[[name]]), ".backup-"),
      tmpdir = dirname(paths[[name]])
    )
  }, character(1L))
  moved_old <- published <- stats::setNames(
    rep(FALSE, length(paths)), names(paths)
  )
  rollback <- TRUE
  on.exit({
    if (rollback) {
      for (name in names(paths)) {
        if (published[[name]] && file.exists(paths[[name]])) {
          unlink(paths[[name]])
        }
      }
      for (name in names(paths)) {
        if (moved_old[[name]] && file.exists(backup_paths[[name]]) &&
            !file.rename(backup_paths[[name]], paths[[name]])) {
          warning("Could not restore original output from ",
                  backup_paths[[name]])
        }
      }
    }
    unlink(tmp_paths)
    if (!rollback) unlink(backup_paths)
  }, add = TRUE)

  saveRDS(combined, tmp_paths[["rds"]])
  write.csv(equal_summary, tmp_paths[["summary"]], row.names = FALSE)
  write.csv(target_summary, tmp_paths[["target_summary"]], row.names = FALSE)
  staged_combined <- readRDS(tmp_paths[["rds"]])
  if (!identical(staged_combined, combined)) {
    stop("Staged combined RDS failed exact round-trip validation.",
         call. = FALSE)
  }
  validate_summary_csv(
    tmp_paths[["summary"]], equal_summary, "equal-target summary CSV"
  )
  validate_summary_csv(
    tmp_paths[["target_summary"]], target_summary, "target summary CSV"
  )
  if (!identical(artifact_snapshot(paths), snapshot)) {
    stop("Combined outputs changed while staging; refusing to overwrite.",
         call. = FALSE)
  }

  for (name in names(paths)) {
    if (snapshot$exists[[name]]) {
      moved_old[[name]] <- file.rename(paths[[name]], backup_paths[[name]])
      if (!moved_old[[name]]) {
        stop("Could not stage existing output for rollback: ", paths[[name]],
             call. = FALSE)
      }
    }
  }
  for (name in names(paths)) {
    published[[name]] <- file.rename(tmp_paths[[name]], paths[[name]])
    if (!published[[name]]) {
      stop("Failed to publish combined output: ", paths[[name]],
           call. = FALSE)
    }
  }
  rollback <- FALSE
  invisible(TRUE)
}

stem <- sprintf(
  "ccle-asym-gcm-effect-grid-s10-pilot-combined-R%d-trees%d",
  PILOT_R, PILOT_RF_TREES
)
rds_path <- file.path(output_dir, paste0(stem, ".rds"))
csv_path <- file.path(output_dir, paste0(stem, "-summary.csv"))
target_csv_path <- file.path(output_dir, paste0(stem, "-target-summary.csv"))
output_paths <- c(
  rds = rds_path,
  summary = csv_path,
  target_summary = target_csv_path
)
publish_combined_outputs(
  combined, equal_target_summary, target_summary, output_paths
)

cat("Combined 30 validated pilot target files.\n")
print(equal_target_summary)
cat("Saved:\n  ", rds_path, "\n  ", csv_path, "\n  ", target_csv_path,
    "\n", sep = "")
