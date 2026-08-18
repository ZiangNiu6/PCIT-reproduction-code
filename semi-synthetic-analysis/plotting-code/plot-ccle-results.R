#!/usr/bin/env Rscript
###############################################################################
# Plot the paired CCLE synthetic-outcome analysis at s = 10.
#
# The preferred input is the matched-v4 result containing 30 null targets and
# 12 additional nonnull settings for three fixed representatives at s = 10.
# The plotter can also validate and use a transitional matched-v3 or matched-v2
# artifact when selected through
# CCLE_RESULTS_SCHEMA=auto|v4|v3|v2 (default: auto, priority v4 > v3 > v2).
# Unmatched legacy results are never accepted. Existing combined-null and QQ
# filenames are preserved. V4 additionally writes three representative-target
# GCM nominal-rejection curves on target-specific effect axes, a two-panel
# introduction preview, and a three-row main-text composite. Pass --intro-only
# to regenerate the standalone EIF5A preview and, for v4, its introduction
# composite; pass --validate-only to validate without writing figures.
###############################################################################

script_args <- commandArgs(trailingOnly = TRUE)
known_args <- c("--intro-only", "--validate-only")
unknown_args <- setdiff(script_args, known_args)
if (length(unknown_args) > 0L) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "),
       call. = FALSE)
}
intro_only <- "--intro-only" %in% script_args
validate_only <- "--validate-only" %in% script_args
if (intro_only && validate_only) {
  stop("Use at most one of --intro-only and --validate-only.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
})

script_file_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_file_arg) != 1L) {
  stop("Could not determine the plotting script path.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_file_arg[[1L]]),
  mustWork = TRUE
)
script_dir <- dirname(script_path)
project_root <- normalizePath(
  file.path(script_dir, "..", "..", ".."),
  mustWork = TRUE
)
output_dir <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables"
)
if (!identical(
      normalizePath(dirname(dirname(dirname(output_dir))), mustWork = FALSE),
      project_root
    )) {
  stop("CCLE output directory is not rooted directly under project_root.",
       call. = FALSE)
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

style_path <- file.path(project_root, "reproduction-code", "plotting-style.R")
if (!file.exists(style_path)) {
  stop("Shared plotting style not found: ", style_path, call. = FALSE)
}
source(style_path)

SCHEMA_V2          <- "ccle_synthy_matched_v2"
SCHEMA_V3          <- "ccle_synthy_matched_v3"
SCHEMA_V4          <- "ccle_synthy_matched_v4"
SIMULATION_DESIGN_V3 <- "paired_null_alternative_s10"
SIMULATION_DESIGN_V4 <- "representative_effect_grid_s10"
CALIBRATION_DESIGN <- "shared_observed_statistic"
EXPECTED_R <- 500L
EXPECTED_B_MAX <- 499L
EXPECTED_H <- 15L
EXPECTED_ALPHA <- 0.05
EXPECTED_S <- 10L
EXPECTED_ALT_EFFECT <- 0.75
EXPECTED_EFFECT_SCALE <- "noise_sd_equals_sd_mu"
EXPECTED_SETUP_ORDER_V4 <-
  "null_all_targets_then_grid_index_major_representatives"
EXPECTED_REPRESENTATIVE_TARGETS_V4 <- data.frame(
  representative_order = 1:3,
  target_idx = c(4L, 15L, 22L),
  target = c("EIF5A_MUT", "PLK1_MUT", "ATP11A_MUT"),
  sparsity_tag = c("low_sparsity", "medium_sparsity", "high_sparsity"),
  stringsAsFactors = FALSE
)
EXPECTED_REPRESENTATIVE_EFFECTS_V4 <- list(
  EIF5A_MUT = c(0, 1, 2, 3, 4),
  PLK1_MUT = c(0, 0.5, 1, 1.5, 2),
  ATP11A_MUT = c(0, 0.5, 1, 1.5, 2)
)
effect_id_from_size_v4 <- function(x) {
  lookup <- c(
    `0` = "delta0", `0.5` = "delta0p5", `1` = "delta1",
    `1.5` = "delta1p5", `2` = "delta2", `3` = "delta3",
    `4` = "delta4"
  )
  ids <- unname(lookup[as.character(x)])
  if (anyNA(ids)) {
    stop("No canonical v4 effect ID for effect size: ",
         paste(x[is.na(ids)], collapse = ", "), call. = FALSE)
  }
  ids
}

results_path_v4 <- file.path(
  data_root,
  "ccle",
  "ccle-typeI-synthY-30targets-matched-v4-results.rds"
)
results_path_v3 <- file.path(
  data_root,
  "ccle",
  "ccle-typeI-synthY-30targets-matched-v3-results.rds"
)
results_path_v2 <- file.path(
  data_root,
  "ccle",
  "ccle-typeI-synthY-30targets-matched-v2-results.rds"
)
schema_request <- tolower(Sys.getenv("CCLE_RESULTS_SCHEMA", unset = "auto"))
valid_schema_requests <- c("auto", "v4", "v3", "v2")
if (!schema_request %in% valid_schema_requests) {
  stop(
    "CCLE_RESULTS_SCHEMA must be one of: ",
    paste(valid_schema_requests, collapse = ", "),
    call. = FALSE
  )
}
schema_paths <- c(
  v4 = results_path_v4,
  v3 = results_path_v3,
  v2 = results_path_v2
)
schema_versions <- c(v4 = SCHEMA_V4, v3 = SCHEMA_V3, v2 = SCHEMA_V2)
if (identical(schema_request, "auto")) {
  available <- names(schema_paths)[file.exists(schema_paths)]
  selected_schema <- if (length(available)) available[[1L]] else NA_character_
} else {
  selected_schema <- schema_request
}
if (is.na(selected_schema) || !file.exists(schema_paths[[selected_schema]])) {
  stop(
    if (!identical(schema_request, "auto")) paste0(
      "The explicitly selected ", schema_request,
      " CCLE result is unavailable: ", schema_paths[[schema_request]], "\n"
    ) else "No matched CCLE result found. Expected one of:\n",
    if (identical(schema_request, "auto")) {
      paste0("  ", paste(unname(schema_paths), collapse = "\n  "), "\n")
    } else "",
    "Run combine-ccle-results.R after the required setup jobs finish.",
    " Unmatched legacy results are intentionally not accepted.",
    call. = FALSE
  )
}
results_path <- unname(schema_paths[[selected_schema]])
input_schema <- unname(schema_versions[[selected_schema]])
v4_manuscript_composites <- file.path(
  output_dir,
  c(
    "ccle-typeI-synthY-intro-calibration-power.pdf",
    "ccle-typeI-synthY-30targets-typeI-power-combined.pdf"
  )
)
if (!identical(selected_schema, "v4") &&
    any(file.exists(v4_manuscript_composites))) {
  stop(
    "Matched-v4 manuscript composite(s) already exist, but schema ",
    selected_schema, " was selected. Refusing to mix stale power panels with ",
    "legacy null-only results. Move the v4 composites out of the figure ",
    "directory or plot the matching v4 artifact.",
    call. = FALSE
  )
}
if (!identical(selected_schema, "v4")) {
  warning(
    "Using explicitly validated transitional matched-", selected_schema,
    " CCLE results; the v4 representative-grid diagnostic is unavailable.",
    call. = FALSE
  )
}
combined <- readRDS(results_path)

storage_methods <- c("asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive")
storage_statistics <- c("gcm", "pcm")
seed_fields <- c(
  "data_seed", "pcm_split_seed", "gcm_permutation_seed",
  "pcm_permutation_seed", "gcm_YonZ_seed", "gcm_XonZ_seed",
  "pcm_YonXZ_seed", "pcm_YonZ_seed", "pcm_YhatonZ_seed",
  "pcm_RonZ_seed"
)
required_combined_fields <- c(
  "schema_version", "calibration_design", "matched_calibration", "n_setups",
  "R", "B_requested", "B_max", "adaptive_h", "adaptive_alpha",
  "alternative", "methods", "statistics", "method_definition",
  "package_versions", "results"
)
expected_n_setups <- switch(
  selected_schema,
  v4 = 42L,
  v3 = 60L,
  v2 = 60L
)
if (!is.list(combined) ||
    length(setdiff(required_combined_fields, names(combined))) ||
    !identical(combined$schema_version, input_schema) ||
    !identical(combined$calibration_design, CALIBRATION_DESIGN) ||
    !isTRUE(combined$matched_calibration) ||
    !identical(as.integer(combined$n_setups), expected_n_setups) ||
    !identical(as.integer(combined$R), EXPECTED_R) ||
    !identical(as.integer(combined$B_requested), EXPECTED_B_MAX) ||
    !identical(as.integer(combined$B_max), EXPECTED_B_MAX) ||
    !identical(as.integer(combined$adaptive_h), EXPECTED_H) ||
    !isTRUE(all.equal(
      as.numeric(combined$adaptive_alpha), EXPECTED_ALPHA, tolerance = 0
    )) ||
    !identical(combined$alternative, "two.sided") ||
    !identical(combined$methods, storage_methods) ||
    !identical(combined$statistics, storage_statistics)) {
  stop(
    "Legacy, malformed, or unmatched combined CCLE result rejected: ",
    results_path,
    call. = FALSE
  )
}
if (identical(input_schema, SCHEMA_V3)) {
  required_v3_combined <- c(
    "simulation_design", "fixed_s", "hypotheses", "effect_sizes",
    "effect_scale"
  )
  if (length(setdiff(required_v3_combined, names(combined))) ||
      !identical(combined$simulation_design, SIMULATION_DESIGN_V3) ||
      !identical(as.integer(combined$fixed_s), EXPECTED_S) ||
      !identical(as.character(combined$hypotheses), c("null", "alternative")) ||
      !isTRUE(all.equal(
        as.numeric(combined$effect_sizes), c(0, EXPECTED_ALT_EFFECT),
        tolerance = 0
      )) ||
      !identical(combined$effect_scale, EXPECTED_EFFECT_SCALE)) {
    stop("Malformed paired matched-v3 CCLE design: ", results_path,
         call. = FALSE)
  }
}
if (identical(input_schema, SCHEMA_V4)) {
  required_v4_combined <- c(
    "simulation_design", "n_targets", "n_representatives",
    "n_null_setups", "n_alternative_setups", "fixed_s",
    "canonical_targets", "representative_targets",
    "representative_effect_grid", "setup_table", "setup_order",
    "paired_seed_indices", "seed_pairing", "effect_scale",
    "summary_by_representative_effect_method",
    "paired_gain_vs_delta0_by_representative_effect_method",
    "null_summary_by_sparsity_method", "null_summary_by_method"
  )
  if (length(setdiff(required_v4_combined, names(combined))) ||
      !identical(combined$simulation_design, SIMULATION_DESIGN_V4) ||
      !identical(as.integer(combined$n_targets), 30L) ||
      !identical(as.integer(combined$n_representatives), 3L) ||
      !identical(as.integer(combined$n_null_setups), 30L) ||
      !identical(as.integer(combined$n_alternative_setups), 12L) ||
      !identical(as.integer(combined$fixed_s), EXPECTED_S) ||
      !identical(combined$setup_order, EXPECTED_SETUP_ORDER_V4) ||
      !identical(as.integer(combined$paired_seed_indices), 31:60) ||
      !is.list(combined$seed_pairing) ||
      !identical(
        names(combined$seed_pairing),
        c(
          "formula", "reference_effect_id", "paired_targets",
          "paired_across_effect_sizes", "component_fields"
        )
      ) ||
      !identical(
        combined$seed_pairing$formula,
        "paired_seed_idx = target_idx + 30"
      ) ||
      !identical(combined$seed_pairing$reference_effect_id, "delta0") ||
      !identical(
        as.character(combined$seed_pairing$paired_targets),
        EXPECTED_REPRESENTATIVE_TARGETS_V4$target
      ) ||
      !isTRUE(combined$seed_pairing$paired_across_effect_sizes) ||
      !identical(combined$seed_pairing$component_fields, seed_fields) ||
      !identical(combined$effect_scale, EXPECTED_EFFECT_SCALE)) {
    stop("Malformed paired matched-v4 CCLE representative-grid design: ",
         results_path,
         call. = FALSE)
  }
}
if (identical(input_schema, SCHEMA_V4)) {
  canonical_targets_v4 <- combined$canonical_targets
  representative_targets_v4 <- combined$representative_targets
  representative_grid_v4 <- combined$representative_effect_grid
  setup_table_v4 <- combined$setup_table
  same_ordered_frame <- function(actual, expected) {
    is.data.frame(actual) && is.data.frame(expected) &&
      identical(names(actual), names(expected)) &&
      isTRUE(all.equal(
        actual, expected, tolerance = 0, check.attributes = FALSE
      ))
  }

  if (!is.data.frame(canonical_targets_v4) ||
      !identical(
        names(canonical_targets_v4),
        c("target_idx", "target", "sparsity_tag", "freq")
      ) ||
      nrow(canonical_targets_v4) != 30L ||
      !identical(as.integer(canonical_targets_v4$target_idx), 1:30) ||
      anyNA(canonical_targets_v4) ||
      anyDuplicated(canonical_targets_v4$target) ||
      any(!is.finite(canonical_targets_v4$freq)) ||
      any(canonical_targets_v4$freq <= 0 | canonical_targets_v4$freq >= 1) ||
      !identical(
        as.integer(table(factor(
          canonical_targets_v4$sparsity_tag,
          levels = c("low_sparsity", "medium_sparsity", "high_sparsity")
        ))),
        rep(10L, 3L)
      )) {
    stop("Malformed canonical-target table in matched-v4 CCLE result.",
         call. = FALSE)
  }

  expected_representative_core <-
    EXPECTED_REPRESENTATIVE_TARGETS_V4[, c(
      "representative_order", "target_idx", "target", "sparsity_tag"
    )]
  if (!is.data.frame(representative_targets_v4) ||
      !identical(
        names(representative_targets_v4),
        c(
          "representative_order", "target_idx", "target", "sparsity_tag",
          "freq"
        )
      ) ||
      nrow(representative_targets_v4) != 3L ||
      !same_ordered_frame(
        representative_targets_v4[, names(expected_representative_core)],
        expected_representative_core
      ) ||
      any(!is.finite(representative_targets_v4$freq))) {
    stop("Malformed representative-target table in matched-v4 CCLE result.",
         call. = FALSE)
  }
  representative_canonical_match <- match(
    representative_targets_v4$target,
    canonical_targets_v4$target
  )
  if (anyNA(representative_canonical_match) ||
      !isTRUE(all.equal(
        representative_targets_v4$freq,
        canonical_targets_v4$freq[representative_canonical_match],
        tolerance = 0, check.attributes = FALSE
      ))) {
    stop("Representative targets disagree with the canonical v4 target table.",
         call. = FALSE)
  }

  expected_representative_grid_v4 <- do.call(
    rbind,
    lapply(1:5, function(grid_idx) {
      effect_sizes <- vapply(
        representative_targets_v4$target,
        function(target) {
          EXPECTED_REPRESENTATIVE_EFFECTS_V4[[target]][[grid_idx]]
        },
        numeric(1L)
      )
      transform(
        representative_targets_v4,
        grid_idx = grid_idx,
        effect_id = effect_id_from_size_v4(effect_sizes),
        effect_size = effect_sizes,
        hypothesis = ifelse(effect_sizes == 0, "null", "alternative")
      )
    })
  )
  rownames(expected_representative_grid_v4) <- NULL
  expected_representative_grid_v4 <-
    expected_representative_grid_v4[, c(
      "representative_order", "target_idx", "target", "sparsity_tag",
      "freq", "grid_idx", "effect_id", "effect_size", "hypothesis"
    )]
  if (!same_ordered_frame(
        representative_grid_v4, expected_representative_grid_v4
      )) {
    stop("Malformed target-specific v4 representative effect grids.",
         call. = FALSE)
  }

  expected_null_setups_v4 <- transform(
    canonical_targets_v4,
    paired_seed_idx = target_idx + 30L,
    s = EXPECTED_S,
    grid_idx = 1L,
    effect_id = "delta0",
    effect_size = 0,
    hypothesis = "null"
  )
  expected_alternative_setups_v4 <- expected_representative_grid_v4[
    expected_representative_grid_v4$grid_idx > 1L,
    c(
      "target_idx", "target", "sparsity_tag", "freq", "grid_idx",
      "effect_id", "effect_size", "hypothesis"
    )
  ]
  expected_alternative_setups_v4$paired_seed_idx <-
    expected_alternative_setups_v4$target_idx + 30L
  expected_alternative_setups_v4$s <- EXPECTED_S
  setup_columns_v4 <- c(
    "target_idx", "paired_seed_idx", "target", "sparsity_tag", "freq", "s",
    "grid_idx", "effect_id", "effect_size", "hypothesis"
  )
  expected_setup_table_v4 <- rbind(
    expected_null_setups_v4[, setup_columns_v4],
    expected_alternative_setups_v4[, setup_columns_v4]
  )
  expected_setup_table_v4$setup_idx <- seq_len(nrow(expected_setup_table_v4))
  expected_setup_table_v4 <- expected_setup_table_v4[, c(
    "setup_idx", setup_columns_v4
  )]
  rownames(expected_setup_table_v4) <- NULL
  if (!same_ordered_frame(setup_table_v4, expected_setup_table_v4)) {
    stop("Malformed exact 42-row v4 setup table.", call. = FALSE)
  }
}
results <- combined$results

required_record_fields <- c(
  "schema_version", "calibration_design", "matched_calibration",
  "setup_idx", "target", "sparsity_tag", "freq", "s", "R",
  "B_requested", "B_max", "adaptive_h", "adaptive_alpha", "alternative",
  "pvals", "statistics",
  "replication_metadata", "method_definition", "package_versions",
  "rejection_rate", "na_rate"
)
required_record_fields_v3 <- c(
  "simulation_design", "target_idx", "hypothesis", "effect_size",
  "effect_scale", "paired_seed_idx"
)
required_record_fields_v4 <- c(
  required_record_fields_v3, "grid_idx", "effect_id"
)
validate_record <- function(record, key) {
  if (!is.list(record) ||
      !all(required_record_fields %in% names(record)) ||
      (identical(input_schema, SCHEMA_V3) &&
       !all(required_record_fields_v3 %in% names(record))) ||
      (identical(input_schema, SCHEMA_V4) &&
       !all(required_record_fields_v4 %in% names(record)))) {
    stop("Malformed matched-calibration CCLE record: ", key,
         call. = FALSE)
  }
  if (!identical(record$schema_version, input_schema) ||
      !identical(record$calibration_design, CALIBRATION_DESIGN) ||
      !isTRUE(record$matched_calibration)) {
    stop("Legacy or unmatched CCLE record rejected: ", key, call. = FALSE)
  }
  if (!identical(as.integer(record$R), EXPECTED_R) ||
      !identical(as.integer(record$B_requested), EXPECTED_B_MAX) ||
      !identical(as.integer(record$B_max), EXPECTED_B_MAX) ||
      !identical(as.integer(record$adaptive_h), EXPECTED_H) ||
      !isTRUE(all.equal(
        as.numeric(record$adaptive_alpha), EXPECTED_ALPHA, tolerance = 0
      )) ||
      !identical(record$alternative, "two.sided") ||
      !is.matrix(record$pvals) || !is.numeric(record$pvals) ||
      !identical(dim(record$pvals), c(EXPECTED_R, length(storage_methods))) ||
      !identical(colnames(record$pvals), storage_methods) ||
      !identical(names(record$rejection_rate), storage_methods) ||
      !identical(names(record$na_rate), storage_methods)) {
    stop("Noncanonical CCLE method/replication schema in: ", key,
         call. = FALSE)
  }
  if (!is.matrix(record$statistics) || !is.numeric(record$statistics) ||
      !identical(dim(record$statistics), c(EXPECTED_R, 2L)) ||
      !identical(colnames(record$statistics), storage_statistics) ||
      any(!is.finite(record$statistics))) {
    stop("Invalid shared observed statistics in: ", key, call. = FALSE)
  }
  expected_key <- if (identical(input_schema, SCHEMA_V4)) {
    sprintf("%s_%s_s%d", record$target, record$effect_id, record$s)
  } else if (identical(input_schema, SCHEMA_V3)) {
    sprintf("%s_%s_s%d", record$target, record$hypothesis, record$s)
  } else {
    sprintf("%s_s%d", record$target, record$s)
  }
  if (length(record$setup_idx) != 1L || is.na(record$setup_idx) ||
      length(record$target) != 1L || is.na(record$target) ||
      length(record$sparsity_tag) != 1L || is.na(record$sparsity_tag) ||
      length(record$freq) != 1L || !is.finite(record$freq) ||
      length(record$s) != 1L || is.na(record$s) ||
      !identical(key, expected_key)) {
    stop("Invalid CCLE setup metadata in: ", key, call. = FALSE)
  }
  if (identical(input_schema, SCHEMA_V3)) {
    expected_effect <- if (identical(record$hypothesis, "null")) {
      0
    } else if (identical(record$hypothesis, "alternative")) {
      EXPECTED_ALT_EFFECT
    } else {
      NA_real_
    }
    if (!identical(record$simulation_design, SIMULATION_DESIGN_V3) ||
        !identical(as.integer(record$s), EXPECTED_S) ||
        length(record$target_idx) != 1L || is.na(record$target_idx) ||
        record$target_idx < 1L || record$target_idx > 30L ||
        length(record$paired_seed_idx) != 1L ||
        !identical(
          as.integer(record$paired_seed_idx),
          as.integer(record$target_idx) + 30L
        ) ||
        !is.finite(expected_effect) ||
        !isTRUE(all.equal(
          as.numeric(record$effect_size), expected_effect, tolerance = 0
        )) ||
        !identical(record$effect_scale, EXPECTED_EFFECT_SCALE)) {
      stop("Invalid paired simulation metadata in: ", key, call. = FALSE)
    }
  }
  if (identical(input_schema, SCHEMA_V4)) {
    record_setup_idx <- as.integer(record$setup_idx)
    expected_setup <- if (
      length(record_setup_idx) == 1L && !is.na(record_setup_idx) &&
      record_setup_idx >= 1L && record_setup_idx <= 42L
    ) expected_setup_table_v4[record_setup_idx, , drop = FALSE] else NULL
    if (!identical(record$simulation_design, SIMULATION_DESIGN_V4) ||
        !identical(as.integer(record$s), EXPECTED_S) ||
        is.null(expected_setup) ||
        !identical(as.integer(record$target_idx),
                   as.integer(expected_setup$target_idx)) ||
        length(record$paired_seed_idx) != 1L ||
        !identical(
          as.integer(record$paired_seed_idx),
          as.integer(expected_setup$paired_seed_idx)
        ) ||
        !identical(as.character(record$target),
                   as.character(expected_setup$target)) ||
        !identical(as.character(record$sparsity_tag),
                   as.character(expected_setup$sparsity_tag)) ||
        !isTRUE(all.equal(
          as.numeric(record$freq), as.numeric(expected_setup$freq),
          tolerance = 0
        )) ||
        !identical(as.integer(record$grid_idx),
                   as.integer(expected_setup$grid_idx)) ||
        !identical(as.character(record$effect_id),
                   as.character(expected_setup$effect_id)) ||
        !isTRUE(all.equal(
          as.numeric(record$effect_size),
          as.numeric(expected_setup$effect_size), tolerance = 0
        )) ||
        !identical(as.character(record$hypothesis),
                   as.character(expected_setup$hypothesis)) ||
        !identical(record$effect_scale, EXPECTED_EFFECT_SCALE)) {
      stop("Invalid paired v4 representative-grid metadata in: ", key,
           call. = FALSE)
    }
  }
  if (any(!is.finite(record$pvals)) ||
      any(record$pvals < 0 | record$pvals > 1) ||
      any(!is.finite(record$rejection_rate)) ||
      any(!is.finite(record$na_rate)) || any(record$na_rate != 0)) {
    stop("Missing or invalid CCLE p-values in: ", key, call. = FALSE)
  }
  expected_asym_gcm <- 2 * stats::pnorm(
    abs(record$statistics[, "gcm"]), lower.tail = FALSE
  )
  expected_asym_pcm <- 2 * stats::pnorm(
    abs(record$statistics[, "pcm"]), lower.tail = FALSE
  )
  if (!isTRUE(all.equal(
        unname(record$pvals[, "asym_gcm"]), unname(expected_asym_gcm),
        tolerance = 1e-14, check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        unname(record$pvals[, "asym_pcm"]), unname(expected_asym_pcm),
        tolerance = 1e-14, check.attributes = FALSE
      ))) {
    stop(
      "Normal p-values do not use the stored shared statistic in: ", key,
      call. = FALSE
    )
  }
  metadata <- record$replication_metadata
  required_metadata <- c(
    "rep", seed_fields, "pcm_D1_n", "pcm_D2_n", "gcm_actual_B",
    "pcm_actual_B",
    "gcm_stop_reason", "pcm_stop_reason", "gcm_error", "pcm_error",
    "gcm_matched_statistic", "pcm_matched_statistic", "gcm_loss_count",
    "pcm_loss_count", "gcm_t_star", "pcm_t_star", "gcm_stopped_early",
    "pcm_stopped_early", "gcm_adaptive_backend", "pcm_adaptive_backend",
    "pcm_D2_indices"
  )
  if (input_schema %in% c(SCHEMA_V3, SCHEMA_V4)) {
    required_metadata <- c(required_metadata, "outcome_noise_sd")
  }
  if (!is.data.frame(metadata) || nrow(metadata) != EXPECTED_R ||
      length(setdiff(required_metadata, names(metadata))) ||
      !identical(as.integer(metadata$rep), seq_len(EXPECTED_R)) ||
      !is.list(metadata$pcm_D2_indices)) {
    stop("Incomplete adaptive metadata in: ", key, call. = FALSE)
  }
  if (input_schema %in% c(SCHEMA_V3, SCHEMA_V4) &&
      (!is.numeric(metadata$outcome_noise_sd) ||
       any(!is.finite(metadata$outcome_noise_sd)) ||
       any(metadata$outcome_noise_sd <= 0))) {
    stop("Invalid outcome-noise scales in: ", key, call. = FALSE)
  }
  invalid_seed <- vapply(seed_fields, function(field) {
    values <- metadata[[field]]
    !is.numeric(values) || anyNA(values) ||
      any(values < 1 | values > .Machine$integer.max |
            values != floor(values)) ||
      anyDuplicated(values) > 0L
  }, logical(1L))
  valid_indices <- vapply(seq_len(EXPECTED_R), function(rep_idx) {
    indices <- metadata$pcm_D2_indices[[rep_idx]]
    D1_n <- metadata$pcm_D1_n[[rep_idx]]
    D2_n <- metadata$pcm_D2_n[[rep_idx]]
    is.numeric(D1_n) && is.finite(D1_n) && D1_n >= 1 &&
      D1_n == floor(D1_n) && is.numeric(D2_n) && is.finite(D2_n) &&
      D2_n >= 1 && D2_n == floor(D2_n) && is.integer(indices) &&
      !anyNA(indices) && length(indices) == D2_n && !anyDuplicated(indices) &&
      all(indices >= 1L) && all(indices <= D1_n + D2_n)
  }, logical(1L))
  if (any(invalid_seed) || !all(valid_indices)) {
    stop("Invalid seed or PCM split metadata in: ", key, call. = FALSE)
  }
  for (prefix in c("gcm", "pcm")) {
    actual_B <- metadata[[paste0(prefix, "_actual_B")]]
    stop_reason <- metadata[[paste0(prefix, "_stop_reason")]]
    loss_count <- metadata[[paste0(prefix, "_loss_count")]]
    t_star <- metadata[[paste0(prefix, "_t_star")]]
    stopped_early <- metadata[[paste0(prefix, "_stopped_early")]]
    backend <- metadata[[paste0(prefix, "_adaptive_backend")]]
    matched <- metadata[[paste0(prefix, "_matched_statistic")]]
    error <- metadata[[paste0(prefix, "_error")]]
    p_value <- record$pvals[, paste0("perm_", prefix)]
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
    if (!is.numeric(actual_B) || anyNA(actual_B) ||
        any(actual_B < 1 | actual_B > EXPECTED_B_MAX |
              actual_B != floor(actual_B)) ||
        !is.numeric(loss_count) || anyNA(loss_count) ||
        any(loss_count < 0 | loss_count > EXPECTED_H |
              loss_count != floor(loss_count)) ||
        !is.numeric(t_star) || anyNA(t_star) ||
        any(!(is.finite(t_star) | is.infinite(t_star))) ||
        !is.logical(stopped_early) || anyNA(stopped_early) ||
        any(stopped_early != (actual_B < EXPECTED_B_MAX)) ||
        !is.character(stop_reason) || anyNA(stop_reason) ||
        any(!stop_reason %in% c("rejection", "futility", "cap")) ||
        !is.character(backend) || anyNA(backend) ||
        any(!backend %in% c("pcmt_native", "ordered_perm_stats_fallback")) ||
        !is.logical(matched) || anyNA(matched) || !all(matched) ||
        !is.character(error) ||
        any(!is.na(error) & nzchar(error)) ||
        any(loss_count < EXPECTED_H &
              !(is.infinite(t_star) & t_star > 0)) ||
        any(loss_count == EXPECTED_H &
              !(is.finite(t_star) & t_star == actual_B)) ||
        !all(valid_reason) ||
        !isTRUE(all.equal(
          unname(p_value), unname(expected_p), tolerance = 1e-12,
          check.attributes = FALSE
        ))) {
      stop("Adaptive ", toupper(prefix),
           " stopping metadata is inconsistent in: ", key, call. = FALSE)
    }
  }
  definition <- record$method_definition
  if (!is.list(definition) ||
      !identical(definition$permutation_mode, "adaptive_anytime_valid") ||
      !identical(definition$alternative, "two.sided") ||
      !identical(as.integer(definition$B_requested), EXPECTED_B_MAX) ||
      !identical(as.integer(definition$B_max), EXPECTED_B_MAX) ||
      !identical(as.integer(definition$h), EXPECTED_H) ||
      !isTRUE(all.equal(
        as.numeric(definition$alpha), EXPECTED_ALPHA, tolerance = 0
      )) ||
      !identical(as.integer(definition$n_cores), 1L) ||
      !identical(as.integer(definition$rf_num_threads), 1L) ||
      !isTRUE(definition$role_specific_rf_seeds) ||
      length(definition$rf_num_trees) != 1L ||
      !is.numeric(definition$rf_num_trees) ||
      !is.finite(definition$rf_num_trees) ||
      definition$rf_num_trees < 1 ||
      definition$rf_num_trees != floor(definition$rf_num_trees)) {
    stop("Invalid matched-calibration method definition in: ", key,
         call. = FALSE)
  }
  versions <- record$package_versions
  if (!is.character(versions) || is.null(names(versions)) ||
      !setequal(names(versions), c("R", "pcmt", "comets", "ranger")) ||
      anyNA(versions) || any(!nzchar(versions))) {
    stop("Missing package-version provenance in: ", key, call. = FALSE)
  }
  expected_rejection <- colMeans(record$pvals <= 0.05)
  if (!isTRUE(all.equal(
    unname(record$rejection_rate), unname(expected_rejection),
    tolerance = 1e-12
  ))) {
    stop("Stored CCLE rejection rates disagree with p-values in: ", key,
         call. = FALSE)
  }
  invisible(TRUE)
}

if (!is.list(results) || length(results) != expected_n_setups ||
    is.null(names(results)) ||
    anyNA(names(results)) || any(names(results) == "") ||
    anyDuplicated(names(results))) {
  stop("The canonical CCLE result must contain ", expected_n_setups,
       " uniquely named setups.",
       call. = FALSE)
}
invisible(Map(validate_record, results, names(results)))
if (!all(vapply(
      results,
      function(record) identical(
        record$method_definition, combined$method_definition
      ),
      logical(1L)
    )) ||
    !all(vapply(
      results,
      function(record) identical(record$package_versions, combined$package_versions),
      logical(1L)
    ))) {
  stop("Combined provenance disagrees with one or more CCLE setup records.",
       call. = FALSE)
}
setup_ids <- vapply(results, function(x) as.integer(x$setup_idx), integer(1L))
targets <- vapply(results, function(x) as.character(x$target), character(1L))
s_values <- vapply(results, function(x) as.integer(x$s), integer(1L))
sparsity_tags <- vapply(
  results, function(x) as.character(x$sparsity_tag), character(1L)
)
valid_tags <- c("low_sparsity", "medium_sparsity", "high_sparsity")
if (!identical(unname(sort(setup_ids)), seq_len(expected_n_setups)) ||
    length(unique(targets)) != 30L ||
    !setequal(unique(sparsity_tags), valid_tags)) {
  stop("The CCLE result does not contain the canonical setup IDs, targets, and frequency bins.",
       call. = FALSE)
}

if (identical(input_schema, SCHEMA_V2)) {
  if (!identical(sort(unique(s_values)), c(5L, 10L)) ||
      any(table(targets) != 2L) ||
      any(table(s_values) != 30L) ||
      any(table(targets, s_values) != 1L) ||
      any(table(sparsity_tags) != 20L)) {
    stop("The matched-v2 result does not match the canonical 30-target x 2-s grid.",
         call. = FALSE)
  }
  all_component_seeds <- unlist(lapply(results, function(record) {
    unlist(record$replication_metadata[seed_fields], use.names = FALSE)
  }), use.names = FALSE)
  if (anyDuplicated(all_component_seeds)) {
    stop("Seed streams overlap across the matched-v2 CCLE setups.",
         call. = FALSE)
  }
  null_results <- results[s_values == EXPECTED_S]
  alternative_results <- list()
} else if (identical(input_schema, SCHEMA_V3)) {
  hypotheses <- vapply(
    results, function(x) as.character(x$hypothesis), character(1L)
  )
  target_indices <- vapply(
    results, function(x) as.integer(x$target_idx), integer(1L)
  )
  if (!identical(unique(s_values), EXPECTED_S) ||
      !setequal(unique(hypotheses), c("null", "alternative")) ||
      any(table(hypotheses) != 30L) ||
      any(table(targets) != 2L) ||
      any(table(targets, hypotheses) != 1L) ||
      !identical(unname(sort(unique(target_indices))), seq_len(30L)) ||
      any(table(target_indices) != 2L) ||
      any(table(sparsity_tags, hypotheses) != 10L)) {
    stop("The matched-v3 result does not match the canonical paired 30-target s=10 design.",
         call. = FALSE)
  }
  for (target_name in unique(targets)) {
    pair <- results[targets == target_name]
    null_record <- pair[[which(vapply(
      pair, function(x) identical(x$hypothesis, "null"), logical(1L)
    ))]]
    alternative_record <- pair[[which(vapply(
      pair, function(x) identical(x$hypothesis, "alternative"), logical(1L)
    ))]]
    paired_fields_match <-
      identical(null_record$target_idx, alternative_record$target_idx) &&
      identical(null_record$paired_seed_idx,
                alternative_record$paired_seed_idx) &&
      identical(null_record$sparsity_tag, alternative_record$sparsity_tag) &&
      isTRUE(all.equal(null_record$freq, alternative_record$freq,
                       tolerance = 0)) &&
      identical(
        null_record$replication_metadata[seed_fields],
        alternative_record$replication_metadata[seed_fields]
      ) &&
      identical(
        null_record$replication_metadata$pcm_D2_indices,
        alternative_record$replication_metadata$pcm_D2_indices
      ) &&
      isTRUE(all.equal(
        null_record$replication_metadata$outcome_noise_sd,
        alternative_record$replication_metadata$outcome_noise_sd,
        tolerance = 0
      ))
    if (!paired_fields_match) {
      stop("Null/alternative pairing metadata disagree for target: ",
           target_name, call. = FALSE)
    }
  }
  null_results <- results[hypotheses == "null"]
  alternative_results <- results[hypotheses == "alternative"]
  unique_pair_seeds <- unlist(lapply(null_results, function(record) {
    unlist(record$replication_metadata[seed_fields], use.names = FALSE)
  }), use.names = FALSE)
  if (anyDuplicated(unique_pair_seeds)) {
    stop("Seed streams overlap across distinct matched-v3 target pairs.",
         call. = FALSE)
  }
} else {
  hypotheses <- vapply(
    results, function(x) as.character(x$hypothesis), character(1L)
  )
  target_indices <- vapply(
    results, function(x) as.integer(x$target_idx), integer(1L)
  )
  grid_indices <- vapply(
    results, function(x) as.integer(x$grid_idx), integer(1L)
  )
  effect_ids <- vapply(
    results, function(x) as.character(x$effect_id), character(1L)
  )
  effect_sizes <- vapply(
    results, function(x) as.numeric(x$effect_size), numeric(1L)
  )
  if (!identical(unique(s_values), EXPECTED_S) ||
      !setequal(unique(hypotheses), c("null", "alternative")) ||
      sum(hypotheses == "null") != 30L ||
      sum(hypotheses == "alternative") != 12L ||
      !identical(unname(setup_ids), seq_len(42L)) ||
      !identical(unname(target_indices),
                 as.integer(expected_setup_table_v4$target_idx)) ||
      !identical(unname(grid_indices),
                 as.integer(expected_setup_table_v4$grid_idx)) ||
      !identical(unname(effect_ids),
                 as.character(expected_setup_table_v4$effect_id)) ||
      !isTRUE(all.equal(
        unname(effect_sizes), as.numeric(expected_setup_table_v4$effect_size),
        tolerance = 0
      )) ||
      !identical(unname(targets),
                 as.character(expected_setup_table_v4$target)) ||
      !identical(unname(sort(unique(target_indices))), seq_len(30L)) ||
      !identical(
        as.integer(table(factor(grid_indices, levels = 1:5))),
        c(30L, 3L, 3L, 3L, 3L)
      )) {
    stop(
      paste(
        "The matched-v4 result does not match the canonical 30-null plus",
        "12-representative-alternative s=10 design."
      ),
      call. = FALSE
    )
  }
  target_counts <- as.integer(table(factor(
    targets, levels = canonical_targets_v4$target
  )))
  expected_target_counts <- ifelse(
    canonical_targets_v4$target %in%
      EXPECTED_REPRESENTATIVE_TARGETS_V4$target,
    5L, 1L
  )
  if (!identical(target_counts, expected_target_counts)) {
    stop("V4 target counts do not match 30 nulls and three five-point grids.",
         call. = FALSE)
  }
  for (target_name in EXPECTED_REPRESENTATIVE_TARGETS_V4$target) {
    target_grid <- results[targets == target_name]
    target_grid_indices <- vapply(
      target_grid, function(x) as.integer(x$grid_idx), integer(1L)
    )
    target_grid <- target_grid[order(target_grid_indices)]
    target_effect_sizes <- vapply(
      target_grid, function(x) as.numeric(x$effect_size), numeric(1L)
    )
    if (!identical(
          unname(vapply(
            target_grid, function(x) as.integer(x$grid_idx), integer(1L)
          )),
          1:5
        ) ||
        !isTRUE(all.equal(
          unname(target_effect_sizes),
          EXPECTED_REPRESENTATIVE_EFFECTS_V4[[target_name]],
          tolerance = 0
        )) ||
        !identical(
          unname(vapply(
            target_grid, function(x) as.character(x$effect_id), character(1L)
          )),
          effect_id_from_size_v4(unname(target_effect_sizes))
        )) {
      stop("Target-specific v4 grid is malformed for: ", target_name,
           call. = FALSE)
    }
    base_record <- target_grid[[1L]]
    paired_fields_match <- vapply(target_grid[-1L], function(record) {
      identical(base_record$target_idx, record$target_idx) &&
        identical(base_record$paired_seed_idx, record$paired_seed_idx) &&
        identical(base_record$sparsity_tag, record$sparsity_tag) &&
        isTRUE(all.equal(base_record$freq, record$freq, tolerance = 0)) &&
        identical(
          base_record$replication_metadata[seed_fields],
          record$replication_metadata[seed_fields]
        ) &&
        identical(
          base_record$replication_metadata$pcm_D2_indices,
          record$replication_metadata$pcm_D2_indices
        ) &&
        identical(
          as.integer(base_record$replication_metadata$pcm_D1_n),
          as.integer(record$replication_metadata$pcm_D1_n)
        ) &&
        identical(
          as.integer(base_record$replication_metadata$pcm_D2_n),
          as.integer(record$replication_metadata$pcm_D2_n)
        ) &&
        isTRUE(all.equal(
          base_record$replication_metadata$outcome_noise_sd,
          record$replication_metadata$outcome_noise_sd,
          tolerance = 0
        ))
    }, logical(1L))
    if (!all(paired_fields_match)) {
      stop("Effect-grid pairing metadata disagree for target: ",
           target_name, call. = FALSE)
    }
  }
  null_results <- results[effect_sizes == 0]
  alternative_results <- results[effect_sizes > 0]
  effect_grid_results <- results[
    targets %in% EXPECTED_REPRESENTATIVE_TARGETS_V4$target
  ]
  unique_pair_seeds <- unlist(lapply(null_results, function(record) {
    unlist(record$replication_metadata[seed_fields], use.names = FALSE)
  }), use.names = FALSE)
  if (anyDuplicated(unique_pair_seeds)) {
    stop("Seed streams overlap across distinct matched-v4 null targets.",
         call. = FALSE)
  }

  record_rejection_v4 <- do.call(rbind, lapply(results, function(record) {
    data.frame(
      target_idx = record$target_idx,
      target = record$target,
      sparsity_tag = record$sparsity_tag,
      freq = record$freq,
      s = record$s,
      grid_idx = record$grid_idx,
      effect_id = record$effect_id,
      effect_size = record$effect_size,
      hypothesis = record$hypothesis,
      method = names(record$rejection_rate),
      rejection_rate = unname(record$rejection_rate),
      na_rate = unname(record$na_rate[names(record$rejection_rate)]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(record_rejection_v4) <- NULL
  record_rejection_v4$representative_order <- match(
    record_rejection_v4$target,
    EXPECTED_REPRESENTATIVE_TARGETS_V4$target
  )
  delta0_rejection_v4 <- record_rejection_v4[
    record_rejection_v4$effect_id == "delta0",
    c("target", "method", "rejection_rate")
  ]
  names(delta0_rejection_v4)[3L] <- "delta0_rejection_rate"
  delta0_match_v4 <- match(
    paste(record_rejection_v4$target, record_rejection_v4$method, sep = "\r"),
    paste(delta0_rejection_v4$target, delta0_rejection_v4$method, sep = "\r")
  )
  if (anyNA(delta0_match_v4)) {
    stop("Every v4 result must match one target/method delta0 reference.",
         call. = FALSE)
  }
  record_rejection_v4$delta0_rejection_rate <-
    delta0_rejection_v4$delta0_rejection_rate[delta0_match_v4]
  record_rejection_v4$paired_gain_vs_delta0 <-
    record_rejection_v4$rejection_rate -
    record_rejection_v4$delta0_rejection_rate

  expected_representative_summary_v4 <- record_rejection_v4[
    !is.na(record_rejection_v4$representative_order),
    c(
      "representative_order", "target_idx", "target", "sparsity_tag",
      "freq", "s", "grid_idx", "effect_id", "effect_size", "hypothesis",
      "method", "rejection_rate", "delta0_rejection_rate",
      "paired_gain_vs_delta0", "na_rate"
    )
  ]
  expected_representative_summary_v4$n_reps <- EXPECTED_R
  expected_representative_summary_v4 <- expected_representative_summary_v4[
    order(
      expected_representative_summary_v4$representative_order,
      expected_representative_summary_v4$grid_idx,
      match(expected_representative_summary_v4$method, storage_methods)
    ),
  ]
  rownames(expected_representative_summary_v4) <- NULL
  if (!same_ordered_frame(
        combined$summary_by_representative_effect_method,
        expected_representative_summary_v4
      )) {
    stop("Stored representative v4 summary disagrees with setup records.",
         call. = FALSE)
  }

  expected_paired_gain_v4 <- expected_representative_summary_v4[
    expected_representative_summary_v4$hypothesis == "alternative",
    c(
      "representative_order", "target_idx", "target", "sparsity_tag",
      "grid_idx", "effect_id", "effect_size", "method", "rejection_rate",
      "delta0_rejection_rate", "paired_gain_vs_delta0", "n_reps"
    )
  ]
  rownames(expected_paired_gain_v4) <- NULL
  if (!same_ordered_frame(
        combined$paired_gain_vs_delta0_by_representative_effect_method,
        expected_paired_gain_v4
      )) {
    stop("Stored v4 paired-gain summary disagrees with setup records.",
         call. = FALSE)
  }

  null_rejection_v4 <- record_rejection_v4[
    record_rejection_v4$hypothesis == "null", , drop = FALSE
  ]
  expected_null_sparsity_v4 <- aggregate(
    rejection_rate ~ sparsity_tag + method,
    data = null_rejection_v4,
    FUN = mean
  )
  names(expected_null_sparsity_v4)[3L] <- "mean_rejection_rate"
  expected_null_sparsity_v4$n_targets <- 10L
  expected_null_sparsity_v4$n_reps_per_target <- EXPECTED_R
  expected_null_sparsity_v4 <- expected_null_sparsity_v4[
    order(
      match(expected_null_sparsity_v4$sparsity_tag,
            c("low_sparsity", "medium_sparsity", "high_sparsity")),
      match(expected_null_sparsity_v4$method, storage_methods)
    ),
  ]
  rownames(expected_null_sparsity_v4) <- NULL
  if (!same_ordered_frame(
        combined$null_summary_by_sparsity_method,
        expected_null_sparsity_v4
      )) {
    stop("Stored v4 null sparsity summary disagrees with setup records.",
         call. = FALSE)
  }

  expected_null_method_v4 <- aggregate(
    rejection_rate ~ method,
    data = null_rejection_v4,
    FUN = mean
  )
  names(expected_null_method_v4)[2L] <- "mean_rejection_rate"
  expected_null_method_v4$n_targets <- 30L
  expected_null_method_v4$n_reps_per_target <- EXPECTED_R
  expected_null_method_v4 <- expected_null_method_v4[
    order(match(expected_null_method_v4$method, storage_methods)),
  ]
  rownames(expected_null_method_v4) <- NULL
  if (!same_ordered_frame(
        combined$null_summary_by_method,
        expected_null_method_v4
      )) {
    stop("Stored v4 null method summary disagrees with setup records.",
         call. = FALSE)
  }
}

if (!exists("effect_grid_results", inherits = FALSE)) {
  effect_grid_results <- list()
}

# Normalize keys for the null-only plotting code below.  This also makes the
# v2 fallback and preferred v3 artifact produce byte-for-byte comparable plot
# selections whenever their null records agree.
names(null_results) <- vapply(null_results, function(record) {
  sprintf("%s_s%d", record$target, record$s)
}, character(1L))
if (length(alternative_results)) {
  names(alternative_results) <- vapply(alternative_results, function(record) {
    if (identical(input_schema, SCHEMA_V4)) {
      sprintf("%s_%s_s%d", record$target, record$effect_id, record$s)
    } else {
      sprintf("%s_s%d", record$target, record$s)
    }
  }, character(1L))
}
results <- null_results

if (length(results) != 30L ||
    any(table(vapply(
      results, function(x) as.character(x$sparsity_tag), character(1L)
    )) != 10L) ||
    !all(c("EIF5A_MUT_s10", "PLK1_MUT_s10", "ATP11A_MUT_s10") %in%
         names(results))) {
  stop("The canonical 30 s=10 null records or representative targets are missing.",
       call. = FALSE)
}
if (validate_only) {
  cat("Validated ", input_schema, " CCLE result: ", results_path, "\n",
      sep = "")
  cat("Validated 30 s=10 null records",
      if (identical(input_schema, SCHEMA_V4)) {
        " and 12 paired alternatives across three representative grids"
      } else if (length(alternative_results)) {
        " and 30 paired alternative records"
      } else "",
      ".\n", sep = "")
  cat("Validated figure output directory: ", output_dir, "\n", sep = "")
  quit(save = "no", status = 0L)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Focus the CCLE figures on the rare-carrier GCM mechanism. PCM remains in the
# validated result schema and reproduction archive, but is not displayed here
# because this null-only benchmark does not cleanly isolate its calibration
# behavior.
method_levels <- c("perm_gcm", "asym_gcm", "naive")
method_labels <- c(
  "asym_gcm" = "GCM (asym)",
  "perm_gcm" = "GCM (adaptive perm)",
  "naive" = "Pearson"
)
stopifnot(all(method_levels %in% storage_methods))

## Internal tags refer to mutation FREQUENCY tercile; display labels invert that
## to X-vector SPARSITY: low-frequency mutations correspond to high-sparsity X.
sparsity_levels <- c("low_sparsity", "medium_sparsity", "high_sparsity")
sparsity_labels <- c("low_sparsity"    = "High sparsity",
                     "medium_sparsity" = "Medium sparsity",
                     "high_sparsity"   = "Low sparsity")

# ---------------------------------------------------------------------------
# 1. Long-format rejection table for the 30 s=10 null setups
# ---------------------------------------------------------------------------
rej_long <- bind_rows(lapply(results, function(r) {
  data.frame(target         = r$target,
             sparsity_tag   = r$sparsity_tag,
             freq           = r$freq,
             s              = r$s,
             method         = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             stringsAsFactors = FALSE)
})) |>
  filter(method %in% method_levels) |>
  mutate(
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    ),
    sparsity_tag = factor(
      sparsity_tag,
      levels = sparsity_levels,
      labels = sparsity_labels[sparsity_levels]
    )
  ) |>
  add_calibration_test_aesthetics()

if (nrow(rej_long) != length(results) * length(method_levels) ||
    anyNA(rej_long$method) ||
    anyNA(rej_long$method_label)) {
  stop("The displayed CCLE result must contain exactly three GCM-focused procedures per setup",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# 2. Pick 3 null-calibration representatives (one per tercile) at s = 10.
#    EIF5A is also the v4 high-sparsity conditional-effect representative, so
#    the null and effect panels use the same target in that sparsity group.
# ---------------------------------------------------------------------------
chosen <- list(
  list(key = "EIF5A_MUT_s10",  tag = "low_sparsity"),
  list(key = "PLK1_MUT_s10",   tag = "medium_sparsity"),
  list(key = "ATP11A_MUT_s10", tag = "high_sparsity")
)
chosen_keys <- vapply(chosen, function(x) x$key, character(1))
stopifnot(all(chosen_keys %in% names(results)))

setup_label <- function(rec) {
  sprintf("%s   (freq = %.3f, s = %d)",
          gsub("_MUT$", "", rec$target), rec$freq, rec$s)
}

setup_order <- vapply(chosen_keys, function(k) setup_label(results[[k]]), character(1))

# ---------------------------------------------------------------------------
# 3. Stopped anytime-value diagnostic for the three chosen setups.
#    The uniform line is descriptive only: adaptive permutation values stopped
#    for rejection/futility are not expected to be Uniform(0,1).
# ---------------------------------------------------------------------------
pv_long <- bind_rows(lapply(chosen_keys, function(k) {
  r  <- results[[k]]
  pv <- as.data.frame(r$pvals)
  pv$rep <- seq_len(nrow(pv))
  pv |>
    pivot_longer(-rep, names_to = "method", values_to = "p_value") |>
    mutate(setup = setup_label(r))
})) |>
  filter(method %in% method_levels, !is.na(p_value)) |>
  mutate(
    setup = factor(setup, levels = setup_order),
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

pv_qq <- pv_long |>
  group_by(setup, method) |>
  arrange(p_value, .by_group = TRUE) |>
  mutate(theoretical = ppoints(n())) |>
  ungroup()

if (nrow(pv_qq) != length(chosen_keys) * EXPECTED_R * length(method_levels) ||
    anyNA(pv_qq$method) ||
    anyNA(pv_qq$method_label)) {
  stop("The CCLE diagnostic must contain exactly three displayed procedures per replicate",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# 3a. V4 representative conditional-effect curves.  Build this once so the
#     same validated target-level rates feed the introduction preview, the
#     standalone diagnostic, and row (c) of the main-text composite.
# ---------------------------------------------------------------------------
effect_grid_summary <- NULL
p_power_core <- NULL
if (identical(input_schema, SCHEMA_V4)) {
  alternative_method_levels <- c("perm_gcm", "asym_gcm")
  effect_grid_long <- bind_rows(lapply(effect_grid_results, function(record) {
    data.frame(
      target = record$target,
      sparsity_tag = record$sparsity_tag,
      freq = record$freq,
      s = record$s,
      grid_idx = record$grid_idx,
      effect_id = record$effect_id,
      effect_size = record$effect_size,
      method = names(record$rejection_rate),
      rejection_rate = unname(record$rejection_rate),
      stringsAsFactors = FALSE
    )
  })) |>
    filter(method %in% alternative_method_levels) |>
    mutate(
      method = factor(method, levels = alternative_method_levels),
      method_label = factor(
        unname(method_labels[as.character(method)]),
        levels = unname(method_labels[alternative_method_levels])
      ),
      target_label = factor(
        target,
        levels = EXPECTED_REPRESENTATIVE_TARGETS_V4$target,
        labels = c(
          "EIF5A (high sparsity)",
          "PLK1 (medium sparsity)",
          "ATP11A (low sparsity)"
        )
      )
    ) |>
    add_calibration_test_aesthetics()

  if (nrow(effect_grid_long) !=
        3L * 5L * length(alternative_method_levels) ||
      anyNA(effect_grid_long$method) ||
      anyNA(effect_grid_long$method_label) ||
      anyNA(effect_grid_long$target_label) ||
      any(table(
        effect_grid_long$target,
        effect_grid_long$grid_idx,
        effect_grid_long$method
      ) != 1L)) {
    stop("The paired v4 representative GCM diagnostic is incomplete.",
         call. = FALSE)
  }

  effect_grid_summary <- effect_grid_long |>
    arrange(target_label, method, grid_idx)

  p_power_core <- ggplot(
    effect_grid_summary,
    aes(
      x = effect_size,
      y = rejection_rate,
      colour = test,
      linetype = calibration,
      shape = calibration,
      group = method_label
    )
  ) +
    geom_hline(
      yintercept = EXPECTED_ALPHA,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.55
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.1, stroke = 0.7) +
    facet_wrap(~ target_label, nrow = 1, scales = "free_x") +
    scale_test_colour(guide = "none") +
    scale_calibration_linetype() +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      x = "Signal strength",
      y = "Rejection rate"
    ) +
    guides(
      linetype = guide_legend(order = 1),
      shape = "none"
    ) +
    theme_figure3(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

# ---------------------------------------------------------------------------
# 3b. Standalone EIF5A calibration preview for the introduction. This uses
#     test colour and calibration border type without repeating EIF5A metadata
#     in a facet strip.
# ---------------------------------------------------------------------------
wilson_interval <- function(successes, repetitions, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  rate <- successes / repetitions
  denominator <- 1 + z^2 / repetitions
  center <- (rate + z^2 / (2 * repetitions)) / denominator
  half_width <- z * sqrt(
    rate * (1 - rate) / repetitions + z^2 / (4 * repetitions^2)
  ) / denominator
  data.frame(
    ci_lower = pmax(0, center - half_width),
    ci_upper = pmin(1, center + half_width)
  )
}

eif5a_null_record <- results[["EIF5A_MUT_s10"]]
if (!is.list(eif5a_null_record) ||
    !all(c("target", "s", "R", "pvals", "rejection_rate") %in%
         names(eif5a_null_record))) {
  stop("Malformed EIF5A Figure 1 result", call. = FALSE)
}
if (!identical(eif5a_null_record$target, "EIF5A_MUT") ||
    !identical(as.integer(eif5a_null_record$s), 10L) ||
    !identical(as.integer(eif5a_null_record$R), 500L)) {
  stop("The Figure 1 result must be EIF5A_MUT, s=10, R=500", call. = FALSE)
}
eif5a_null_pvals <- as.data.frame(eif5a_null_record$pvals)[method_levels]
if (nrow(eif5a_null_pvals) != eif5a_null_record$R ||
    any(!is.finite(as.matrix(eif5a_null_pvals))) ||
    any(as.matrix(eif5a_null_pvals) < 0 |
        as.matrix(eif5a_null_pvals) > 1)) {
  stop("Figure 1 contains missing or invalid p-values", call. = FALSE)
}
eif5a_null_repetitions <- colSums(!is.na(eif5a_null_pvals))
eif5a_null_rejections <- colSums(eif5a_null_pvals <= 0.05, na.rm = TRUE)
eif5a_null_summary <- data.frame(
  method = method_levels,
  repetitions = unname(eif5a_null_repetitions),
  rejections = unname(eif5a_null_rejections),
  rejection_rate = unname(eif5a_null_rejections / eif5a_null_repetitions),
  stringsAsFactors = FALSE
) |>
  bind_cols(wilson_interval(
    eif5a_null_rejections,
    eif5a_null_repetitions
  )) |>
  mutate(
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

stopifnot(isTRUE(all.equal(
  unname(eif5a_null_summary$rejection_rate),
  unname(eif5a_null_record$rejection_rate[method_levels])
)))

theme_eif5a_intro <- theme_bw(base_size = 9) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, colour = "black"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(0.28, "in"),
    legend.spacing.x = grid::unit(0.01, "in"),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 7.2, margin = margin(r = 2, l = 4)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0.06, "in"),
    plot.margin = margin(3, 2, 8, 2)
  )

eif5a_null_bar_upper <- min(
  1,
  1.20 * max(eif5a_null_summary$ci_upper, 0.05)
)
p_eif5a_null_bar <- ggplot(
  eif5a_null_summary,
  aes(x = method_label, y = rejection_rate)
) +
  geom_col(
    aes(fill = test, linetype = calibration),
    colour = "black",
    linewidth = 0.55,
    width = 0.64
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.13,
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.55
  ) +
  geom_text(
    aes(y = ci_upper, label = sprintf("%.3f", rejection_rate)),
    vjust = -0.5,
    size = 2.5
  ) +
  scale_test_fill() +
  scale_calibration_linetype() +
  scale_y_continuous(
    limits = c(0, eif5a_null_bar_upper),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Rejection rate") +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(
        colour = "black",
        linetype = "solid",
        linewidth = 0.55
      )
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(
        fill = "white",
        colour = "black",
        linewidth = 0.8
      )
    )
  ) +
  theme_eif5a_intro +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
eif5a_calibration_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-EIF5A-calibration.pdf"
)
ggsave(
  eif5a_calibration_path,
  p_eif5a_null_bar,
  width = 5.5,
  height = 3.4
)
cat("Saved: ", eif5a_calibration_path, "\n", sep = "")

if (identical(input_schema, SCHEMA_V4)) {
  eif5a_effect_summary <- effect_grid_summary |>
    filter(target == "EIF5A_MUT") |>
    arrange(method, grid_idx)
  if (nrow(eif5a_effect_summary) !=
        5L * length(alternative_method_levels) ||
      any(table(
        eif5a_effect_summary$grid_idx,
        eif5a_effect_summary$method
      ) != 1L)) {
    stop("The EIF5A introduction power panel is incomplete.", call. = FALSE)
  }

  p_eif5a_intro <- ggplot(
    eif5a_effect_summary,
    aes(
      x = effect_size,
      y = rejection_rate,
      colour = test,
      linetype = calibration,
      shape = calibration,
      group = method_label
    )
  ) +
    geom_hline(
      yintercept = EXPECTED_ALPHA,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.55
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.9, stroke = 0.65) +
    scale_test_colour(guide = "none") +
    scale_calibration_linetype(guide = "none") +
    scale_x_continuous(
      limits = c(0, 4),
      breaks = 0:4,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "EIF5A conditional-effect response",
      x = "Signal strength",
      y = "Rejection rate"
    ) +
    theme_eif5a_intro +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
      plot.margin = margin(3, 4, 8, 2)
    )

  intro_legend <- cowplot::get_plot_component(
    p_eif5a_null_bar + theme(legend.position = "bottom"),
    "guide-box-bottom"
  )
  p_intro_panels <- cowplot::plot_grid(
    p_eif5a_null_bar +
      labs(title = "EIF5A null calibration") +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5)
      ),
    p_eif5a_intro,
    ncol = 2,
    align = "hv",
    axis = "tb",
    labels = c("(a)", "(b)"),
    label_size = 10,
    label_fontface = "bold",
    rel_widths = c(0.46, 0.54)
  )
  p_intro_combined <- cowplot::plot_grid(
    p_intro_panels,
    intro_legend,
    ncol = 1,
    rel_heights = c(1, 0.14)
  )
  intro_combined_path <- file.path(
    output_dir,
    "ccle-typeI-synthY-intro-calibration-power.pdf"
  )
  ggsave(
    intro_combined_path,
    p_intro_combined,
    width = 8.4,
    height = 3.6
  )
  cat("Saved: ", intro_combined_path, "\n", sep = "")
}

if (intro_only) {
  quit(save = "no", status = 0L)
}

p_qq <- ggplot(
  pv_qq,
  aes(
    x = theoretical,
    y = p_value,
    colour = test,
    linetype = calibration,
    group = method_label
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_grid(method_label ~ setup) +
  scale_test_colour(guide = "none") +
  scale_calibration_linetype(guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "CCLE GCM null-calibration diagnostic",
    subtitle = paste0(
      "Adaptive GCM permutation values use h = ", EXPECTED_H,
      " and Bmax = ", EXPECTED_B_MAX,
      "; the 45-degree line is descriptive, not a uniform-null reference."
    ),
    x = "Uniform(0,1) descriptive quantile",
    y = "Observed p-value quantile"
  ) +
  theme_figure3(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

qq_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-30targets-pvalqq.pdf"
)
ggsave(qq_path, p_qq, width = 12, height = 7.5)
cat("Saved: ", qq_path, "\n", sep = "")
cat(
  "  Adaptive GCM permutation values are stopped anytime-valid values; the ",
  "uniform line is descriptive and must not be interpreted as their null ",
  "distribution.\n",
  sep = ""
)

# ---------------------------------------------------------------------------
# 4. Settings shared by the combined figure's all-setup boxplot.
# ---------------------------------------------------------------------------
n_setups_per_facet <- rej_long |>
  filter(method == method_levels[1]) |>
  count(sparsity_tag)

## Pseudo-log y-axis (asinh-style; handles the few asym-GCM rejection rates
## that are exactly zero) with breaks at log decades, which are evenly spaced
## visually. The 0.05 nominal level is flagged by the dashed black line.
typeI_y_scale <- scale_y_continuous(
  trans  = scales::pseudo_log_trans(sigma = 0.001, base = 10),
  breaks = c(0.001, 0.01, 0.1, 1),
  labels = c("0.001", "0.01", "0.1", "1"),
  limits = c(0, 1),
  expand = expansion(mult = c(0.02, 0.02))
)

# ---------------------------------------------------------------------------
# 5. Rejection rates for the combined figure's representative-setting bars.
# ---------------------------------------------------------------------------
rej_chosen <- bind_rows(lapply(chosen_keys, function(k) {
  r <- results[[k]]
  data.frame(setup          = setup_label(r),
             method         = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             stringsAsFactors = FALSE)
})) |>
  filter(method %in% method_levels) |>
  mutate(
    setup = factor(setup, levels = setup_order),
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

if (nrow(rej_chosen) != length(chosen_keys) * length(method_levels) ||
    anyNA(rej_chosen$method) ||
    anyNA(rej_chosen$method_label)) {
  stop("Representative CCLE panels must contain exactly three displayed procedures per setup",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# 6. Combined two-row figure: 30 null setups on top, 3 setups below.
#    Shared test-colour and calibration-linetype legends appear at the bottom;
#    x-axis labels are omitted because the legends identify all procedures.
# ---------------------------------------------------------------------------
strip_x <- theme(axis.text.x  = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank())

p_box_c <- ggplot(
  rej_long,
  aes(x = method_label, y = rejection_rate, fill = test)
) +
  geom_boxplot(
    aes(linetype = calibration),
    colour = "black",
    width = 0.65,
    linewidth = 0.55,
    outlier.size = 0.8,
    alpha = 0.85
  ) +
  geom_jitter(aes(colour = test), width = 0.12, height = 0,
              size = 0.55, alpha = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  facet_wrap(~ sparsity_tag, nrow = 1) +
  scale_test_fill() +
  scale_test_colour(guide = "none") +
  scale_calibration_linetype() +
  typeI_y_scale +
  labs(y = expression("Empirical Type-I rate at " * alpha * " = 0.05")) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(alpha = 1, colour = "black", linetype = "solid")
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(fill = "white", colour = "black", alpha = 1)
    )
  ) +
  theme_figure3(base_size = 13) +
  theme(strip.text       = element_text(face = "bold", size = 11),
        legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.text      = element_text(size = 11),
        legend.key       = element_rect(colour = NA, fill = NA)) +
  strip_x

p_bar_c <- ggplot(
  rej_chosen,
  aes(x = method_label, y = rejection_rate, fill = test)
) +
  geom_col(
    aes(linetype = calibration),
    colour = "black",
    linewidth = 0.55,
    width = 0.7
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  facet_wrap(~ setup, nrow = 1) +
  scale_test_fill() +
  scale_calibration_linetype() +
  scale_y_continuous(limits = c(0, min(
                       1, 1.15 * max(rej_chosen$rejection_rate)
                     )),
                     breaks = scales::breaks_pretty(n = 5),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(y = expression("Rejection rate at " * alpha * " = 0.05")) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(alpha = 1, colour = "black", linetype = "solid")
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(fill = "white", colour = "black", alpha = 1)
    )
  ) +
  theme_figure3(base_size = 13) +
  theme(strip.text       = element_text(size = 11, face = "bold"),
        legend.box       = "horizontal") +
  strip_x

# Extract the explicitly bottom-positioned guide box.  ggplot2 4 creates
# placeholder guide boxes for several positions, so matching the generic
# "guide-box" can select an empty placeholder instead of the visible legend.
legend_plot <- p_bar_c +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.width = grid::unit(0.40, "in"),
    legend.key = element_rect(colour = NA, fill = NA)
  )
shared_legend <- cowplot::get_plot_component(
  legend_plot,
  "guide-box-bottom"
)

p_box_nl <- p_box_c + theme(legend.position = "none")
p_bar_nl <- p_bar_c + theme(legend.position = "none")

p_panels <- cowplot::plot_grid(
  p_box_nl, p_bar_nl,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("(a)", "(b)"),
  label_size = 12,
  label_fontface = "bold",
  rel_heights = c(1, 1)
)

p_combined <- cowplot::plot_grid(
  p_panels, shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.10)
)

combined_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-30targets-typeI-combined.pdf"
)
ggsave(combined_path, p_combined, width = 12, height = 7.8)
cat("Saved: ", combined_path, "\n", sep = "")

# ---------------------------------------------------------------------------
# 6b. Main-text composite (v4): a 2 x 3 grid whose ROWS are the two extreme
#     mutation-sparsity groups and whose COLUMNS are the three panel types.
#     The medium-sparsity group is deliberately omitted here; it remains in the
#     standalone diagnostics and in the appendix figures.
# ---------------------------------------------------------------------------
if (identical(input_schema, SCHEMA_V4)) {
  p_power_row <- p_power_core +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank()
    )

  ## Row keys: display sparsity label, its representative target, and the
  ## representative's setup label used by the null-calibration bar panel.
  composite_rows <- list(
    list(
      sparsity = sparsity_labels[["low_sparsity"]],
      target   = "EIF5A_MUT",
      setup    = setup_label(results[["EIF5A_MUT_s10"]])
    ),
    list(
      sparsity = sparsity_labels[["high_sparsity"]],
      target   = "ATP11A_MUT",
      setup    = setup_label(results[["ATP11A_MUT_s10"]])
    )
  )

  ## The sparsity group names each row from the left margin, so the group is
  ## stated once per row instead of on the column-(a) facet strip.
  row_side_label <- function(label) {
    cowplot::ggdraw() +
      cowplot::draw_label(
        label,
        angle = 90,
        fontface = "bold",
        size = 13
      )
  }

  composite_row_grid <- function(row, labels) {
    box_data <- rej_long[
      as.character(rej_long$sparsity_tag) == row$sparsity, , drop = FALSE
    ]
    bar_data <- rej_chosen[
      as.character(rej_chosen$setup) == row$setup, , drop = FALSE
    ]
    power_data <- effect_grid_summary[
      effect_grid_summary$target == row$target, , drop = FALSE
    ]
    if (!nrow(box_data) || !nrow(bar_data) || !nrow(power_data)) {
      stop("Empty composite row for sparsity group: ", row$sparsity,
           call. = FALSE)
    }
    ## The row label already names the sparsity group and its representative,
    ## so each strip states what its own column measures.  Keeping a strip on
    ## every column also keeps the three panel heights identical within a row.
    target_short <- sub("_MUT$", "", row$target)
    box_data$panel_label <- "All 10 null targets"
    bar_data$panel_label <- paste0(target_short, ": Type-I error")
    power_data$panel_label <- paste0(target_short, ": Power")
    panels <- cowplot::plot_grid(
      p_box_nl + box_data + facet_wrap(~ panel_label, nrow = 1),
      p_bar_nl + bar_data + facet_wrap(~ panel_label, nrow = 1),
      p_power_row + power_data + facet_wrap(~ panel_label, nrow = 1),
      nrow = 1,
      align = "h",
      axis = "tb",
      labels = labels,
      label_size = 12,
      label_fontface = "bold",
      rel_widths = c(1, 1, 1.15)
    )
    cowplot::plot_grid(
      row_side_label(row$sparsity),
      panels,
      nrow = 1,
      rel_widths = c(0.045, 1)
    )
  }

  p_two_row_panels <- cowplot::plot_grid(
    composite_row_grid(composite_rows[[1L]], c("(a)", "(b)", "(c)")),
    composite_row_grid(composite_rows[[2L]], NULL),
    ncol = 1,
    rel_heights = c(1, 1)
  )
  p_two_row_combined <- cowplot::plot_grid(
    p_two_row_panels,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.09)
  )
  combined_power_path <- file.path(
    output_dir,
    "ccle-typeI-synthY-30targets-typeI-power-combined.pdf"
  )
  ggsave(
    combined_power_path,
    p_two_row_combined,
    width = 13.5,
    height = 7.6
  )
  cat("Saved: ", combined_power_path, "\n", sep = "")
}

# ---------------------------------------------------------------------------
# 7. Paired effect diagnostic (v4 representative grids; v3 transitional arm).
# ---------------------------------------------------------------------------
if (identical(input_schema, SCHEMA_V4)) {
  p_alternative <- p_power_core +
    labs(
      title = "CCLE paired conditional-effect experiment",
      subtitle = paste(
        "s = 10; one target per sparsity group;",
        "target-specific effect axes; nominal cutoff = 0.05"
      ),
      caption = paste(
        paste(
          "GCM only. Points are target-level raw nominal-cutoff rejection,",
          "including the paired null at delta = 0."
        ),
        paste(
          "Curves are not size-adjusted. EIF5A was selected in screening-only",
          "R = 25 pilot checks with wide Monte Carlo uncertainty; pilot runs are not pooled."
        ),
        sep = "\n"
      )
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )

  alternative_path <- file.path(
    output_dir,
    "ccle-alternative-synthY-30targets-rejection.pdf"
  )
  ggsave(alternative_path, p_alternative, width = 12.5, height = 5.0)
  cat("Saved: ", alternative_path, "\n", sep = "")
} else if (length(alternative_results)) {
  alternative_method_levels <- c("perm_gcm", "asym_gcm")
  alternative_long <- bind_rows(lapply(alternative_results, function(record) {
    data.frame(
      target = record$target,
      sparsity_tag = record$sparsity_tag,
      freq = record$freq,
      s = record$s,
      effect_size = record$effect_size,
      method = names(record$rejection_rate),
      rejection_rate = unname(record$rejection_rate),
      stringsAsFactors = FALSE
    )
  })) |>
    filter(method %in% alternative_method_levels) |>
    mutate(
      method = factor(method, levels = alternative_method_levels),
      method_label = factor(
        unname(method_labels[as.character(method)]),
        levels = unname(method_labels[alternative_method_levels])
      ),
      sparsity_tag = factor(
        sparsity_tag,
        levels = sparsity_levels,
        labels = sparsity_labels[sparsity_levels]
      )
    ) |>
    add_calibration_test_aesthetics()

  if (nrow(alternative_long) != 30L * length(alternative_method_levels) ||
      anyNA(alternative_long$method) || anyNA(alternative_long$method_label) ||
      any(alternative_long$effect_size != EXPECTED_ALT_EFFECT)) {
    stop("The paired conditional-alternative diagnostic is incomplete.",
         call. = FALSE)
  }

  p_alternative <- ggplot(
    alternative_long,
    aes(x = method_label, y = rejection_rate, fill = test)
  ) +
    geom_boxplot(
      aes(linetype = calibration),
      colour = "black",
      width = 0.65,
      linewidth = 0.55,
      outlier.size = 0.8,
      alpha = 0.85
    ) +
    geom_jitter(
      aes(colour = test),
      width = 0.12,
      height = 0,
      size = 0.7,
      alpha = 0.65,
      show.legend = FALSE
    ) +
    facet_wrap(~ sparsity_tag, nrow = 1) +
    scale_test_fill() +
    scale_test_colour(guide = "none") +
    scale_calibration_linetype() +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "CCLE paired conditional-alternative experiment",
      subtitle = paste0(
        "s = 10; per-carrier shift = ", EXPECTED_ALT_EFFECT,
        " outcome-noise SD; nominal cutoff = 0.05"
      ),
      x = NULL,
      y = "Rejection rate under conditional alternative",
      caption = paste(
        "GCM only. Null rejection rates are reported separately;",
        "these rates are not size-adjusted."
      )
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        order = 1,
        override.aes = list(alpha = 1, colour = "black", linetype = "solid")
      ),
      linetype = guide_legend(
        nrow = 1,
        byrow = TRUE,
        order = 2,
        override.aes = list(fill = "white", colour = "black", alpha = 1)
      )
    ) +
    theme_figure3(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 10, colour = "black"),
      legend.position = "none",
      legend.box = "horizontal",
      plot.title = element_text(face = "bold")
    )

  alternative_path <- file.path(
    output_dir,
    "ccle-alternative-synthY-30targets-rejection.pdf"
  )
  ggsave(alternative_path, p_alternative, width = 12, height = 4.8)
  cat("Saved: ", alternative_path, "\n", sep = "")
} else {
  cat(
    "Skipped paired effect-rejection figure: the matched-v2 fallback contains null records only.\n",
    sep = ""
  )
}

# ---------------------------------------------------------------------------
# 8. Console summary
# ---------------------------------------------------------------------------
cat("\n=== Setups per sparsity facet (one row per method) ===\n")
print(n_setups_per_facet)

cat("\n=== GCM-focused mean / median Type-I rate by method x sparsity ===\n")
rej_long |>
  group_by(sparsity_tag, method_label) |>
  summarise(mean   = mean(rejection_rate),
            median = median(rejection_rate),
            min    = min(rejection_rate),
            max    = max(rejection_rate),
            .groups = "drop") |>
  print(n = Inf)

if (identical(input_schema, SCHEMA_V4)) {
  cat("\n=== Representative-target GCM nominal rejection by effect ===\n")
  effect_grid_summary |>
    arrange(target_label, method_label, grid_idx) |>
    as_tibble() |>
    print(n = Inf)
} else if (length(alternative_results)) {
  cat("\n=== Transitional v3 GCM rejection under the alternative ===\n")
  alternative_long |>
    group_by(sparsity_tag, method_label) |>
    summarise(mean = mean(rejection_rate), .groups = "drop") |>
    print(n = Inf)
}

cat("\n=== Three chosen GCM-focused setups for the diagnostic ===\n")
for (k in chosen_keys) {
  r <- results[[k]]
  cat(sprintf("  %-18s  tag=%-15s  freq=%.4f  s=%2d\n",
              r$target, r$sparsity_tag, r$freq, r$s))
  displayed_rates <- r$rejection_rate[method_levels]
  cat("    displayed rejection rates: ",
      paste(sprintf("%s=%.3f", names(displayed_rates), displayed_rates),
            collapse = "  "), "\n", sep = "")
}
