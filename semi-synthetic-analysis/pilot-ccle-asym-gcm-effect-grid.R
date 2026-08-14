#!/usr/bin/env Rscript
###############################################################################
# CCLE asymptotic-GCM effect-size pilot: 30 fixed targets, s = 10.
#
# Usage:
#   Rscript pilot-ccle-asym-gcm-effect-grid.R 1
#   Rscript pilot-ccle-asym-gcm-effect-grid.R 1:15
#   Rscript pilot-ccle-asym-gcm-effect-grid.R 1,3,8:12
#   Rscript pilot-ccle-asym-gcm-effect-grid.R all
#
# Environment:
#   PILOT_R           replicates per target (default 25)
#   PILOT_RF_TREES    trees per nuisance forest (default 500)
#   PILOT_RF_THREADS  ranger threads per forest (default 4)
#   PILOT_OUTPUT_DIR  output directory (default
#                     $LOCAL_PDRT_DATA_DIR/ccle/pilot-asym-gcm-grid)
#
# One atomic RDS and CSV are written per target. This pilot has its own schema
# and namespace and is never consumed by the matched-v4 production combiner.
###############################################################################

suppressPackageStartupMessages({
  library(pcmt)
  library(comets)
  library(ranger)
})

LEGACY_SCHEMA_VERSION <- "ccle_asym_gcm_effect_pilot_v1"
INTERMEDIATE_SCHEMA_VERSION <- "ccle_asym_gcm_effect_pilot_v2"
SCHEMA_VERSION <- "ccle_asym_gcm_effect_pilot_v3"
SIMULATION_DESIGN <- "paired_effect_grid_s10"
CALIBRATION <- "normal_reference"
METHOD <- "asym_gcm_rf_cached_x"
FIXED_S <- 10L
LEGACY_EFFECT_SIZES <- c(0, 0.25, 0.5, 0.75, 1, 1.25)
INTERMEDIATE_EFFECT_SIZES <- c(LEGACY_EFFECT_SIZES, 1.5)
EFFECT_SIZES <- c(INTERMEDIATE_EFFECT_SIZES, 2)
SOURCE_EFFECT_GRIDS <- stats::setNames(
  list(LEGACY_EFFECT_SIZES, INTERMEDIATE_EFFECT_SIZES),
  c(LEGACY_SCHEMA_VERSION, INTERMEDIATE_SCHEMA_VERSION)
)
SAMPLE_SEED <- 2026L
RUN_SEED_BASE <- 30260000L
N_PER_BIN <- 10L
ALPHA <- 0.05

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not determine this script's path", call. = FALSE)
}
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "matched-calibration-utils.R"))

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
ccle_dir <- file.path(data_root, "ccle")

output_dir <- Sys.getenv(
  "PILOT_OUTPUT_DIR",
  unset = file.path(ccle_dir, "pilot-asym-gcm-grid")
)
output_dir <- path.expand(output_dir)
if (!dir.exists(output_dir) &&
    !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create PILOT_OUTPUT_DIR: ", output_dir, call. = FALSE)
}
output_dir <- normalizePath(output_dir, mustWork = TRUE)

parse_target_selection <- function(value, n_targets = 30L) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("Provide a target selection: 1..30, comma/range syntax, or all.",
         call. = FALSE)
  }
  value <- gsub("[[:space:]]", "", value)
  if (tolower(value) == "all") {
    return(seq_len(n_targets))
  }
  pieces <- strsplit(value, ",", fixed = TRUE)[[1L]]
  if (!length(pieces) || any(!nzchar(pieces))) {
    stop("Invalid target selection: ", value, call. = FALSE)
  }
  parsed <- unlist(lapply(pieces, function(piece) {
    if (grepl("^[0-9]+$", piece)) {
      return(as.integer(piece))
    }
    if (grepl("^[0-9]+[:-][0-9]+$", piece)) {
      bounds <- as.integer(strsplit(piece, "[:-]")[[1L]])
      if (bounds[1L] > bounds[2L]) {
        stop("Descending target ranges are not allowed: ", piece,
             call. = FALSE)
      }
      return(seq.int(bounds[1L], bounds[2L]))
    }
    stop("Invalid target-selection token: ", piece, call. = FALSE)
  }), use.names = FALSE)
  if (anyNA(parsed) || any(parsed < 1L | parsed > n_targets)) {
    stop("Target indices must lie in 1..", n_targets, call. = FALSE)
  }
  unique(as.integer(parsed))
}

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) != 1L) {
  stop(
    paste0(
      "Usage: Rscript pilot-ccle-asym-gcm-effect-grid.R ",
      "<1..30|comma/range selection|all>"
    ),
    call. = FALSE
  )
}
target_selection <- parse_target_selection(cli_args[[1L]])

# ---- Exact production preprocessing and deterministic target selection ----
input_files <- file.path(ccle_dir, c("mutation.txt", "response.csv"))
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Missing CCLE input(s):\n  ", paste(missing_inputs, collapse = "\n  "),
       call. = FALSE)
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
n <- nrow(M_full)
freqs <- colMeans(M_full)

qs <- quantile(freqs, probs = c(1 / 3, 2 / 3), names = FALSE)
pool_low <- sort(names(freqs)[freqs <= qs[1L]])
pool_med <- sort(names(freqs)[freqs > qs[1L] & freqs <= qs[2L]])
pool_high <- sort(names(freqs)[freqs > qs[2L]])
if (any(c(length(pool_low), length(pool_med), length(pool_high)) < N_PER_BIN)) {
  stop("A frequency bin has fewer than 10 screened mutations.", call. = FALSE)
}
set.seed(SAMPLE_SEED)
targets_low <- sort(sample(pool_low, N_PER_BIN))
targets_med <- sort(sample(pool_med, N_PER_BIN))
targets_high <- sort(sample(pool_high, N_PER_BIN))
TARGETS_DF <- rbind(
  data.frame(target = targets_low, sparsity_tag = "low_sparsity",
             stringsAsFactors = FALSE),
  data.frame(target = targets_med, sparsity_tag = "medium_sparsity",
             stringsAsFactors = FALSE),
  data.frame(target = targets_high, sparsity_tag = "high_sparsity",
             stringsAsFactors = FALSE)
)
TARGETS_DF$target_idx <- seq_len(nrow(TARGETS_DF))
TARGETS_DF$freq <- unname(freqs[TARGETS_DF$target])
TARGETS_DF$carrier_count <- as.integer(colSums(M_full[, TARGETS_DF$target,
                                                     drop = FALSE]))
stopifnot(nrow(TARGETS_DF) == 30L)

# ---- Exact production generator ----
generate_Y_components <- function(Z, s) {
  p <- ncol(Z)
  s <- max(1L, min(s, p))
  S <- sort(sample.int(p, size = s))
  beta <- rnorm(s, mean = 0, sd = 1)
  mu <- as.numeric(Z[, S, drop = FALSE] %*% beta)
  if (s >= 2L) {
    n_pairs_max <- as.integer(choose(s, 2L))
    n_inter <- min(2L * s, n_pairs_max)
    pairs <- utils::combn(S, 2)
    if (ncol(pairs) > n_inter) {
      pairs <- pairs[, sample.int(ncol(pairs), size = n_inter), drop = FALSE]
    }
    gamma <- rnorm(ncol(pairs), mean = 0, sd = 1)
    for (k in seq_len(ncol(pairs))) {
      mu <- mu + gamma[k] * Z[, pairs[1L, k]] * Z[, pairs[2L, k]]
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

studentized_statistic <- function(rY, rX) {
  product <- as.numeric(rY) * as.numeric(rX)
  product_mean <- mean(product)
  product_variance <- mean(product^2) - product_mean^2
  if (!is.finite(product_variance) || product_variance <= 0) return(0)
  sqrt(length(product)) * product_mean / sqrt(product_variance)
}

fit_x_residuals <- function(X, Z, seed) {
  fit <- ranger::ranger(
    y = factor(X),
    x = Z,
    probability = TRUE,
    seed = as.integer(seed),
    num.threads = PILOT_RF_THREADS,
    num.trees = PILOT_RF_TREES
  )
  predictions <- predict(fit, data = Z)$predictions
  if (!is.matrix(predictions) || ncol(predictions) != 2L) {
    stop("Binary X|Z ranger fit did not return two probability columns.",
         call. = FALSE)
  }
  as.numeric(X) - as.numeric(predictions[, 2L])
}

fit_y_residuals <- function(Y, Z, seed) {
  fit <- ranger::ranger(
    y = as.numeric(Y),
    x = Z,
    seed = as.integer(seed),
    num.threads = PILOT_RF_THREADS,
    num.trees = PILOT_RF_TREES
  )
  as.numeric(Y) - as.numeric(predict(fit, data = Z)$predictions)
}

package_versions <- c(
  R = R.version.string,
  pcmt = as.character(utils::packageVersion("pcmt")),
  comets = as.character(utils::packageVersion("comets")),
  ranger = as.character(utils::packageVersion("ranger"))
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
RAW_FIELDS <- c(
  "target_idx", "target", "sparsity_tag", "freq", "carrier_count", "rep",
  "paired_seed_idx", "data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "outcome_noise_sd", "effect_size", "statistic", "p_value", "rejected",
  "error", "elapsed_seconds"
)
RESULT_FIELDS_V1 <- c(
  "schema_version", "simulation_design", "calibration", "method",
  "fixed_s", "effect_sizes", "alpha", "R", "rf_num_trees",
  "rf_num_threads", "seed_base", "seed_formula", "target_idx",
  "paired_seed_idx", "target", "sparsity_tag", "freq", "carrier_count",
  "n", "p", "method_definition", "preflight", "package_versions",
  "summary", "raw"
)
RESULT_FIELDS_VERSIONED <- append(
  RESULT_FIELDS_V1,
  "resume_provenance",
  after = match("package_versions", RESULT_FIELDS_V1)
)
PREFLIGHT_FIELDS <- c(
  "passed", "tolerance", "cached_statistic", "pcmt_statistic",
  "statistic_difference", "max_abs_rY_difference",
  "max_abs_rX_difference", "pcmt_B"
)
SUMMARY_FIELDS <- c(
  "effect_size", "rejection_rate", "mcse", "n_success"
)
RESUME_PROVENANCE_FIELDS <- c(
  "mode", "source_schema_version", "retained_effect_sizes",
  "computed_effect_sizes", "source_rds_md5", "source_csv_md5"
)

same_numeric <- function(x, y, tolerance = 0) {
  isTRUE(all.equal(
    as.numeric(x), as.numeric(y), tolerance = tolerance,
    check.attributes = FALSE
  ))
}

make_rejection_summary <- function(raw, effect_sizes) {
  result <- aggregate(rejected ~ effect_size, data = raw, FUN = mean)
  names(result)[2L] <- "rejection_rate"
  result <- result[match(effect_sizes, result$effect_size), , drop = FALSE]
  rownames(result) <- NULL
  result$mcse <- sqrt(
    result$rejection_rate * (1 - result$rejection_rate) / PILOT_R
  )
  result$n_success <- PILOT_R
  result
}

validate_preflight <- function(preflight, context) {
  if (!is.list(preflight) || !identical(names(preflight), PREFLIGHT_FIELDS) ||
      !isTRUE(preflight$passed) ||
      !is.finite(preflight$tolerance) || preflight$tolerance <= 0 ||
      !all(is.finite(c(
        preflight$cached_statistic, preflight$pcmt_statistic,
        preflight$statistic_difference, preflight$max_abs_rY_difference,
        preflight$max_abs_rX_difference
      ))) ||
      !same_numeric(
        preflight$statistic_difference,
        preflight$pcmt_statistic - preflight$cached_statistic,
        tolerance = 1e-15
      ) ||
      abs(preflight$statistic_difference) > preflight$tolerance ||
      preflight$max_abs_rY_difference < 0 ||
      preflight$max_abs_rY_difference > preflight$tolerance ||
      preflight$max_abs_rX_difference < 0 ||
      preflight$max_abs_rX_difference > preflight$tolerance ||
      !identical(as.integer(preflight$pcmt_B), 1L)) {
    stop("Invalid cached-ranger preflight in ", context, call. = FALSE)
  }
  invisible(TRUE)
}

validate_resume_provenance <- function(provenance, schema_version, context) {
  if (!is.list(provenance) ||
      !identical(names(provenance), RESUME_PROVENANCE_FIELDS) ||
      length(provenance$mode) != 1L) {
    stop("Invalid resume provenance in ", context, call. = FALSE)
  }
  if (identical(schema_version, INTERMEDIATE_SCHEMA_VERSION)) {
    allowed_modes <- c(
      "fresh_seven_effect_grid", "upgraded_legacy_six_effect_grid"
    )
    if (!provenance$mode %in% allowed_modes) {
      stop("Invalid resume provenance in ", context, call. = FALSE)
    }
    if (identical(provenance$mode, "fresh_seven_effect_grid")) {
      valid <- length(provenance$source_schema_version) == 1L &&
        is.na(provenance$source_schema_version) &&
        length(provenance$retained_effect_sizes) == 0L &&
        identical(
          as.numeric(provenance$computed_effect_sizes),
          INTERMEDIATE_EFFECT_SIZES
        ) && length(provenance$source_rds_md5) == 1L &&
        is.na(provenance$source_rds_md5) &&
        length(provenance$source_csv_md5) == 1L &&
        is.na(provenance$source_csv_md5)
    } else {
      md5_pattern <- "^[[:xdigit:]]{32}$"
      valid <- identical(
        provenance$source_schema_version, LEGACY_SCHEMA_VERSION
      ) && identical(
        as.numeric(provenance$retained_effect_sizes), LEGACY_EFFECT_SIZES
      ) && identical(
        as.numeric(provenance$computed_effect_sizes), 1.5
      ) && length(provenance$source_rds_md5) == 1L &&
        grepl(md5_pattern, provenance$source_rds_md5) &&
        length(provenance$source_csv_md5) == 1L &&
        grepl(md5_pattern, provenance$source_csv_md5)
    }
  } else if (identical(schema_version, SCHEMA_VERSION) &&
             identical(provenance$mode, "fresh_eight_effect_grid")) {
    valid <- length(provenance$source_schema_version) == 1L &&
      is.na(provenance$source_schema_version) &&
      length(provenance$retained_effect_sizes) == 0L &&
      identical(as.numeric(provenance$computed_effect_sizes), EFFECT_SIZES) &&
      length(provenance$source_rds_md5) == 1L &&
      is.na(provenance$source_rds_md5) &&
      length(provenance$source_csv_md5) == 1L &&
      is.na(provenance$source_csv_md5)
  } else if (identical(schema_version, SCHEMA_VERSION) &&
             identical(provenance$mode, "upgraded_effect_grid_prefix")) {
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
    valid <- !is.null(source_effects) &&
      identical(
        as.numeric(provenance$retained_effect_sizes), source_effects
      ) && identical(
        as.numeric(provenance$computed_effect_sizes), missing_effects
      ) && length(provenance$source_rds_md5) == 1L &&
      grepl(md5_pattern, provenance$source_rds_md5) &&
      length(provenance$source_csv_md5) == 1L &&
      grepl(md5_pattern, provenance$source_csv_md5)
  } else {
    valid <- FALSE
  }
  if (!isTRUE(valid)) {
    stop("Invalid resume provenance in ", context, call. = FALSE)
  }
  invisible(TRUE)
}

validate_raw_grid <- function(raw, target_idx, target_row, effect_sizes,
                              context) {
  n_effects <- length(effect_sizes)
  expected_rep <- rep(seq_len(PILOT_R), each = n_effects)
  expected_effect <- rep(effect_sizes, times = PILOT_R)
  paired_seed_idx <- target_idx + 30L
  if (!is.data.frame(raw) || !identical(names(raw), RAW_FIELDS) ||
      nrow(raw) != PILOT_R * n_effects ||
      !identical(as.integer(raw$rep), expected_rep) ||
      !identical(as.numeric(raw$effect_size), expected_effect) ||
      any(raw$target_idx != target_idx) ||
      any(raw$target != target_row$target) ||
      any(raw$sparsity_tag != target_row$sparsity_tag) ||
      !same_numeric(raw$freq, rep(target_row$freq, nrow(raw))) ||
      any(raw$carrier_count != target_row$carrier_count) ||
      any(raw$paired_seed_idx != paired_seed_idx) ||
      any(!is.na(raw$error)) || any(!is.finite(raw$statistic)) ||
      any(!is.finite(raw$p_value)) || any(raw$p_value < 0 | raw$p_value > 1) ||
      anyNA(raw$rejected) || any(raw$rejected != (raw$p_value <= ALPHA)) ||
      any(!is.finite(raw$elapsed_seconds)) || any(raw$elapsed_seconds < 0) ||
      any(!is.finite(raw$outcome_noise_sd)) ||
      any(raw$outcome_noise_sd <= 0) ||
      !same_numeric(
        raw$p_value,
        2 * stats::pnorm(abs(raw$statistic), lower.tail = FALSE),
        tolerance = 1e-15
      )) {
    stop("Invalid or incomplete effect-grid rows in ", context,
         call. = FALSE)
  }

  expected_seeds <- lapply(seq_len(PILOT_R), function(rep_idx) {
    ccle_seed_schedule(paired_seed_idx, rep_idx, RUN_SEED_BASE)
  })
  expected_data_seed <- rep(vapply(
    expected_seeds, `[[`, integer(1L), "data_seed"
  ), each = n_effects)
  expected_y_seed <- rep(vapply(
    expected_seeds, `[[`, integer(1L), "gcm_YonZ_seed"
  ), each = n_effects)
  expected_x_seed <- rep(vapply(
    expected_seeds, `[[`, integer(1L), "gcm_XonZ_seed"
  ), each = n_effects)
  if (!identical(as.integer(raw$data_seed), expected_data_seed) ||
      !identical(as.integer(raw$gcm_YonZ_seed), expected_y_seed) ||
      !identical(as.integer(raw$gcm_XonZ_seed), expected_x_seed)) {
    stop("Exact seed schedule mismatch in ", context, call. = FALSE)
  }

  paired_fields <- c(
    "data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed", "outcome_noise_sd"
  )
  paired_ok <- all(vapply(split(raw, raw$rep), function(block) {
    all(vapply(paired_fields, function(field) {
      length(unique(block[[field]])) == 1L
    }, logical(1L)))
  }, logical(1L)))
  canonical <- raw[raw$effect_size == effect_sizes[[1L]], ]
  seed_values <- unlist(
    canonical[c("data_seed", "gcm_YonZ_seed", "gcm_XonZ_seed")],
    use.names = FALSE
  )
  if (!paired_ok || nrow(canonical) != PILOT_R || anyDuplicated(seed_values)) {
    stop("Within-target effect pairing or seed uniqueness failed in ",
         context, call. = FALSE)
  }
  invisible(TRUE)
}

validate_result_record <- function(record, target_idx, target_row,
                                   effect_sizes, schema_version, context) {
  expected_fields <- if (identical(schema_version, LEGACY_SCHEMA_VERSION)) {
    RESULT_FIELDS_V1
  } else {
    RESULT_FIELDS_VERSIONED
  }
  if (!is.list(record) || !identical(names(record), expected_fields) ||
      !identical(record$schema_version, schema_version) ||
      !identical(record$simulation_design, SIMULATION_DESIGN) ||
      !identical(record$calibration, CALIBRATION) ||
      !identical(record$method, METHOD) ||
      !identical(as.integer(record$fixed_s), FIXED_S) ||
      !identical(as.numeric(record$effect_sizes), effect_sizes) ||
      !same_numeric(record$alpha, ALPHA) ||
      !identical(as.integer(record$R), PILOT_R) ||
      !identical(as.integer(record$rf_num_trees), PILOT_RF_TREES) ||
      !identical(as.integer(record$rf_num_threads), PILOT_RF_THREADS) ||
      !identical(as.integer(record$seed_base), RUN_SEED_BASE) ||
      !identical(record$seed_formula, SEED_FORMULA) ||
      !identical(as.integer(record$target_idx), target_idx) ||
      !identical(as.integer(record$paired_seed_idx), target_idx + 30L) ||
      !identical(record$target, target_row$target) ||
      !identical(record$sparsity_tag, target_row$sparsity_tag) ||
      !same_numeric(record$freq, target_row$freq) ||
      !identical(as.integer(record$carrier_count), target_row$carrier_count) ||
      !identical(as.integer(record$n), n) ||
      !identical(as.integer(record$p), ncol(M_full)) ||
      !identical(record$method_definition, METHOD_DEFINITION) ||
      !identical(record$package_versions, package_versions)) {
    stop("Pilot configuration or provenance mismatch in ", context,
         call. = FALSE)
  }
  validate_preflight(record$preflight, context)
  validate_raw_grid(record$raw, target_idx, target_row, effect_sizes, context)
  expected_summary <- make_rejection_summary(record$raw, effect_sizes)
  if (!is.data.frame(record$summary) ||
      !identical(names(record$summary), SUMMARY_FIELDS) ||
      !isTRUE(all.equal(
        record$summary, expected_summary, tolerance = 1e-15,
        check.attributes = FALSE
      ))) {
    stop("Pilot rejection summary mismatch in ", context, call. = FALSE)
  }
  if (!identical(schema_version, LEGACY_SCHEMA_VERSION)) {
    validate_resume_provenance(
      record$resume_provenance, schema_version, context
    )
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
  comparable_fields <- setdiff(RAW_FIELDS, "error")
  if (!identical(names(csv), RAW_FIELDS) || nrow(csv) != nrow(raw) ||
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

artifact_snapshot <- function(paths) {
  exists <- file.exists(paths)
  md5 <- rep(NA_character_, length(paths))
  if (any(exists)) {
    md5[exists] <- unname(tools::md5sum(paths[exists]))
  }
  names(exists) <- names(md5) <- names(paths)
  list(exists = exists, md5 = md5)
}

same_artifact_snapshot <- function(x, y) {
  identical(x$exists, y$exists) && identical(x$md5, y$md5)
}

load_existing_artifact <- function(paths, target_idx, target_row) {
  snapshot <- artifact_snapshot(paths)
  if (xor(snapshot$exists[["rds"]], snapshot$exists[["csv"]])) {
    stop(
      "Refusing to resume target ", target_idx,
      ": exactly one of the expected RDS/CSV artifacts exists.",
      call. = FALSE
    )
  }
  if (!snapshot$exists[["rds"]]) {
    return(list(mode = "fresh", snapshot = snapshot, record = NULL))
  }
  record <- tryCatch(
    readRDS(paths[["rds"]]),
    error = function(e) {
      stop("Could not read existing pilot RDS: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  context <- paste0("existing target ", target_idx, " artifact")
  if (identical(record$schema_version, SCHEMA_VERSION)) {
    validate_result_record(
      record, target_idx, target_row, EFFECT_SIZES, SCHEMA_VERSION, context
    )
    validate_csv_copy(paths[["csv"]], record$raw, context)
    return(list(mode = "complete", snapshot = snapshot, record = record))
  }
  source_effects <- if (length(record$schema_version) == 1L &&
                        !is.na(record$schema_version)) {
    SOURCE_EFFECT_GRIDS[[record$schema_version]]
  } else {
    NULL
  }
  if (!is.null(source_effects)) {
    validate_result_record(
      record, target_idx, target_row, source_effects,
      record$schema_version, context
    )
    validate_csv_copy(paths[["csv"]], record$raw, context)
    return(list(mode = "resume", snapshot = snapshot, record = record))
  }
  stop(
    "Existing target ", target_idx,
    " artifact has an unsupported schema and was left untouched.",
    call. = FALSE
  )
}

validate_outcome_components <- function(outcome, expected_noise_sd = NULL,
                                        context) {
  if (!is.list(outcome) ||
      !identical(names(outcome), c("mu", "eps", "outcome_noise_sd")) ||
      length(outcome$mu) != n || length(outcome$eps) != n ||
      any(!is.finite(outcome$mu)) || any(!is.finite(outcome$eps)) ||
      !is.finite(outcome$outcome_noise_sd) || outcome$outcome_noise_sd <= 0 ||
      !same_numeric(
        outcome$outcome_noise_sd, stats::sd(outcome$mu), tolerance = 1e-14
      ) ||
      (!is.null(expected_noise_sd) && !same_numeric(
        outcome$outcome_noise_sd, expected_noise_sd, tolerance = 1e-14
      ))) {
    stop("Regenerated outcome components failed validation for ", context,
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_x_residuals <- function(rX, X, context) {
  if (!is.numeric(rX) || length(rX) != length(X) || any(!is.finite(rX)) ||
      any(rX < -1 - 1e-12 | rX > 1 + 1e-12) ||
      any(X == 0 & rX > 1e-12) || any(X == 1 & rX < -1e-12)) {
    stop("Regenerated cached X|Z residuals failed validation for ", context,
         call. = FALSE)
  }
  invisible(TRUE)
}

publish_artifact_pair <- function(result, paths, target_idx, target_row,
                                  expected_snapshot) {
  if (!same_artifact_snapshot(artifact_snapshot(paths), expected_snapshot)) {
    stop(
      "Target ", target_idx,
      " artifacts changed during computation; refusing to overwrite them.",
      call. = FALSE
    )
  }
  tmp_rds <- tempfile(
    pattern = paste0(basename(paths[["rds"]]), ".staged-"),
    tmpdir = dirname(paths[["rds"]])
  )
  tmp_csv <- tempfile(
    pattern = paste0(basename(paths[["csv"]]), ".staged-"),
    tmpdir = dirname(paths[["csv"]])
  )
  backup_rds <- tempfile(
    pattern = paste0(basename(paths[["rds"]]), ".backup-"),
    tmpdir = dirname(paths[["rds"]])
  )
  backup_csv <- tempfile(
    pattern = paste0(basename(paths[["csv"]]), ".backup-"),
    tmpdir = dirname(paths[["csv"]])
  )
  rollback <- TRUE
  moved_old_rds <- moved_old_csv <- FALSE
  published_rds <- published_csv <- FALSE
  on.exit({
    if (rollback) {
      if (published_rds && file.exists(paths[["rds"]])) unlink(paths[["rds"]])
      if (published_csv && file.exists(paths[["csv"]])) unlink(paths[["csv"]])
      if (moved_old_rds && file.exists(backup_rds) &&
          !file.rename(backup_rds, paths[["rds"]])) {
        warning("Could not restore the original RDS from ", backup_rds)
      }
      if (moved_old_csv && file.exists(backup_csv) &&
          !file.rename(backup_csv, paths[["csv"]])) {
        warning("Could not restore the original CSV from ", backup_csv)
      }
    }
    unlink(c(tmp_rds, tmp_csv))
    if (!rollback) unlink(c(backup_rds, backup_csv))
  }, add = TRUE)

  saveRDS(result, tmp_rds)
  write.csv(result$raw, tmp_csv, row.names = FALSE)
  staged <- readRDS(tmp_rds)
  validate_result_record(
    staged, target_idx, target_row, EFFECT_SIZES, SCHEMA_VERSION,
    paste0("staged target ", target_idx, " artifact")
  )
  validate_csv_copy(
    tmp_csv, staged$raw, paste0("staged target ", target_idx, " artifact")
  )
  if (!same_artifact_snapshot(artifact_snapshot(paths), expected_snapshot)) {
    stop(
      "Target ", target_idx,
      " artifacts changed while staging; refusing to overwrite them.",
      call. = FALSE
    )
  }

  if (expected_snapshot$exists[["rds"]]) {
    moved_old_rds <- file.rename(paths[["rds"]], backup_rds)
    if (!moved_old_rds) stop("Could not stage the original RDS for rollback.",
                             call. = FALSE)
    moved_old_csv <- file.rename(paths[["csv"]], backup_csv)
    if (!moved_old_csv) stop("Could not stage the original CSV for rollback.",
                             call. = FALSE)
  }
  published_rds <- file.rename(tmp_rds, paths[["rds"]])
  if (!published_rds) stop("Failed to publish the upgraded RDS.",
                           call. = FALSE)
  published_csv <- file.rename(tmp_csv, paths[["csv"]])
  if (!published_csv) stop("Failed to publish the upgraded CSV.",
                           call. = FALSE)
  rollback <- FALSE
  invisible(TRUE)
}

run_target <- function(target_idx) {
  target_row <- TARGETS_DF[target_idx, ]
  target <- target_row$target
  X <- as.numeric(M_full[, target])
  Z <- M_full[, setdiff(mut_names, target), drop = FALSE]
  paired_seed_idx <- target_idx + 30L
  stem <- sprintf(
    "ccle-asym-gcm-effect-grid-s10-pilot-target-%02d-%s-R%d-trees%d",
    target_idx, target, PILOT_R, PILOT_RF_TREES
  )
  paths <- c(
    rds = file.path(output_dir, paste0(stem, ".rds")),
    csv = file.path(output_dir, paste0(stem, ".csv"))
  )
  existing <- load_existing_artifact(paths, target_idx, target_row)
  if (identical(existing$mode, "complete")) {
    cat(sprintf(
      "\n[Target %d/30] %s already has a validated eight-effect artifact; skipping.\n",
      target_idx, target
    ))
    return(invisible(existing$record))
  }

  effects_to_compute <- if (identical(existing$mode, "resume")) {
    setdiff(EFFECT_SIZES, existing$record$effect_sizes)
  } else {
    EFFECT_SIZES
  }
  retained_raw <- if (identical(existing$mode, "resume")) {
    existing$record$raw
  } else {
    NULL
  }
  n_rows <- PILOT_R * length(effects_to_compute)
  rows <- vector("list", n_rows)
  row_idx <- 0L
  preflight <- if (identical(existing$mode, "resume")) {
    existing$record$preflight
  } else {
    NULL
  }
  target_t0 <- proc.time()[["elapsed"]]

  cat(sprintf(
    "\n[Target %d/30] %s; carriers=%d; R=%d; computing effects=%s%s\n",
    target_idx, target, target_row$carrier_count, PILOT_R,
    paste(effects_to_compute, collapse = ","),
    if (identical(existing$mode, "resume")) {
      paste0(
        "; retaining validated effects=",
        paste(existing$record$effect_sizes, collapse = ",")
      )
    } else {
      ""
    }
  ))

  for (rep_idx in seq_len(PILOT_R)) {
    seeds <- ccle_seed_schedule(paired_seed_idx, rep_idx, RUN_SEED_BASE)
    outcome <- ccle_with_seed(
      seeds[["data_seed"]],
      generate_Y_components(Z, FIXED_S)
    )
    expected_noise_sd <- NULL
    if (identical(existing$mode, "resume")) {
      legacy_block <- retained_raw[retained_raw$rep == rep_idx, ]
      expected_noise_sd <- legacy_block$outcome_noise_sd[[1L]]
    }
    validate_outcome_components(
      outcome, expected_noise_sd,
      paste0("target ", target_idx, ", replicate ", rep_idx)
    )
    rX <- fit_x_residuals(X, Z, seeds[["gcm_XonZ_seed"]])
    validate_x_residuals(
      rX, X, paste0("target ", target_idx, ", replicate ", rep_idx)
    )

    for (effect_idx in seq_along(effects_to_compute)) {
      effect_size <- effects_to_compute[[effect_idx]]
      Y <- outcome$mu + outcome$eps +
        effect_size * outcome$outcome_noise_sd * X
      fit_t0 <- proc.time()[["elapsed"]]
      error <- NA_character_
      statistic <- p_value <- NA_real_
      attempt <- tryCatch({
        rY <- fit_y_residuals(Y, Z, seeds[["gcm_YonZ_seed"]])
        statistic <- studentized_statistic(rY, rX)
        p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
        list(rY = rY, statistic = statistic, p_value = p_value)
      }, error = function(e) {
        error <<- conditionMessage(e)
        NULL
      })
      elapsed <- proc.time()[["elapsed"]] - fit_t0

      if (identical(existing$mode, "fresh") && rep_idx == 1L &&
          effect_idx == 1L) {
        if (is.null(attempt)) {
          stop("Cannot run preflight after cached fit failure: ", error,
               call. = FALSE)
        }
        preflight_args_y <- ccle_rf_args(
          seeds[["gcm_YonZ_seed"]], PILOT_RF_TREES
        )
        preflight_args_y$num.threads <- PILOT_RF_THREADS
        preflight_args_x <- ccle_rf_args(
          seeds[["gcm_XonZ_seed"]], PILOT_RF_TREES
        )
        preflight_args_x$num.threads <- PILOT_RF_THREADS
        preflight_fit <- ccle_with_seed(
          seeds[["gcm_permutation_seed"]],
          pcmt::perm_gcm(
            Y = Y,
            X = X,
            Z = Z,
            reg_YonZ = "rf",
            reg_XonZ = "rf",
            args_YonZ = preflight_args_y,
            args_XonZ = preflight_args_x,
            B = 1L,
            n_cores = 1L,
            alternative = "two.sided"
          )
        )
        statistic_difference <-
          as.numeric(preflight_fit$statistic) - attempt$statistic
        residual_y_difference <- max(abs(
          as.numeric(preflight_fit$rY) - attempt$rY
        ))
        residual_x_difference <- max(abs(
          as.numeric(preflight_fit$rX) - rX
        ))
        tolerance <- 1e-12
        preflight <- list(
          passed = isTRUE(abs(statistic_difference) <= tolerance) &&
            isTRUE(residual_y_difference <= tolerance) &&
            isTRUE(residual_x_difference <= tolerance),
          tolerance = tolerance,
          cached_statistic = attempt$statistic,
          pcmt_statistic = as.numeric(preflight_fit$statistic),
          statistic_difference = statistic_difference,
          max_abs_rY_difference = residual_y_difference,
          max_abs_rX_difference = residual_x_difference,
          pcmt_B = as.integer(preflight_fit$parameter[["B"]])
        )
        if (!isTRUE(preflight$passed) || preflight$pcmt_B != 1L) {
          stop("Cached public-ranger implementation failed pcmt preflight.",
               call. = FALSE)
        }
      }

      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        target_idx = target_idx,
        target = target,
        sparsity_tag = target_row$sparsity_tag,
        freq = target_row$freq,
        carrier_count = target_row$carrier_count,
        rep = rep_idx,
        paired_seed_idx = paired_seed_idx,
        data_seed = seeds[["data_seed"]],
        gcm_YonZ_seed = seeds[["gcm_YonZ_seed"]],
        gcm_XonZ_seed = seeds[["gcm_XonZ_seed"]],
        outcome_noise_sd = outcome$outcome_noise_sd,
        effect_size = effect_size,
        statistic = statistic,
        p_value = p_value,
        rejected = is.finite(p_value) && p_value <= ALPHA,
        error = error,
        elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
    }
    if (rep_idx %% 5L == 0L || rep_idx == PILOT_R) {
      cat(sprintf("  rep %d/%d; elapsed %.1f min\n", rep_idx, PILOT_R,
                  (proc.time()[["elapsed"]] - target_t0) / 60))
    }
  }

  computed_raw <- do.call(rbind, rows)
  raw <- if (is.null(retained_raw)) {
    computed_raw
  } else {
    rbind(retained_raw, computed_raw)
  }
  raw <- raw[order(raw$rep, match(raw$effect_size, EFFECT_SIZES)), ]
  rownames(raw) <- NULL
  validate_raw_grid(
    raw, target_idx, target_row, EFFECT_SIZES,
    paste0("computed target ", target_idx, " artifact")
  )

  rejection_summary <- make_rejection_summary(raw, EFFECT_SIZES)
  resume_provenance <- if (identical(existing$mode, "resume")) {
    list(
      mode = "upgraded_effect_grid_prefix",
      source_schema_version = existing$record$schema_version,
      retained_effect_sizes = existing$record$effect_sizes,
      computed_effect_sizes = effects_to_compute,
      source_rds_md5 = existing$snapshot$md5[["rds"]],
      source_csv_md5 = existing$snapshot$md5[["csv"]]
    )
  } else {
    list(
      mode = "fresh_eight_effect_grid",
      source_schema_version = NA_character_,
      retained_effect_sizes = numeric(),
      computed_effect_sizes = EFFECT_SIZES,
      source_rds_md5 = NA_character_,
      source_csv_md5 = NA_character_
    )
  }

  result <- list(
    schema_version = SCHEMA_VERSION,
    simulation_design = SIMULATION_DESIGN,
    calibration = CALIBRATION,
    method = METHOD,
    fixed_s = FIXED_S,
    effect_sizes = EFFECT_SIZES,
    alpha = ALPHA,
    R = PILOT_R,
    rf_num_trees = PILOT_RF_TREES,
    rf_num_threads = PILOT_RF_THREADS,
    seed_base = RUN_SEED_BASE,
    seed_formula = SEED_FORMULA,
    target_idx = target_idx,
    paired_seed_idx = paired_seed_idx,
    target = target,
    sparsity_tag = target_row$sparsity_tag,
    freq = target_row$freq,
    carrier_count = target_row$carrier_count,
    n = n,
    p = ncol(M_full),
    method_definition = METHOD_DEFINITION,
    preflight = preflight,
    package_versions = package_versions,
    resume_provenance = resume_provenance,
    summary = rejection_summary,
    raw = raw
  )

  validate_result_record(
    result, target_idx, target_row, EFFECT_SIZES, SCHEMA_VERSION,
    paste0("in-memory target ", target_idx, " artifact")
  )
  publish_artifact_pair(
    result, paths, target_idx, target_row, existing$snapshot
  )
  cat("  Saved:\n    ", paths[["rds"]], "\n    ", paths[["csv"]], "\n",
      sep = "")
  invisible(result)
}

cat(sprintf(
  paste0(
    "Loaded CCLE: n=%d, p=%d; selected target indices: %s\n",
    "Pilot: R=%d, trees=%d, RF threads=%d; output=%s\n"
  ),
  n, ncol(M_full), paste(target_selection, collapse = ","), PILOT_R,
  PILOT_RF_TREES, PILOT_RF_THREADS, output_dir
))
for (target_idx in target_selection) run_target(target_idx)
cat("\nPilot target selection completed successfully.\n")
