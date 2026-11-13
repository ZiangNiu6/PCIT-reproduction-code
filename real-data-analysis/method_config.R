## Shared method and output-schema configuration for the Frangieh score study.
##
## Keep these values synchronized with the Frangieh section of the manuscript. Both
## the per-replicate runner and the downstream scripts source this file so a
## result cannot silently change method order or plotting colours.

PERTURB_ANALYSIS_SCHEMA <- "gb2024"
PERTURB_SCHEMA_VERSION <- "perturb_score5_v4_gb2024_q10"
PERTURB_RESULTS_SUBDIR <- "score-calibration-gb2024-q10"
PERTURB_INTERMEDIATE_SUBDIR <- "intermediate-data-gb2024-single-guide"
PERTURB_PLOTS_SUBDIR <- "plots"
PERTURB_FIGURES_SUBDIR <- "figures-and-tables"
PERTURB_COHORT_POLICY <- "gb2024_single_guide_historical_stored_counts_v1"
PERTURB_INTERMEDIATE_SCHEMA <- "frangieh_intermediate_bundle_v2"
PERTURB_PROVENANCE_SCHEMA <- "frangieh_io_provenance_v3"
PERTURB_Q_LEVELS <- 0.10
PERTURB_HEADLINE_Q <- 0.10
PERTURB_RNG_LAYOUT_ID <- "frangieh_paired_reference_slot2_v1"
PERTURB_RNG_KIND <- c("Mersenne-Twister", "Inversion", "Rejection")
PERTURB_NULL_GENE_POLICY <- "exclude_positive_control_targets"

PERTURB_METHOD_LEVELS <- c(
  "parametric_score",
  "information_studentized_permutation",
  "studentized_permutation",
  "hgf_signflip_mc",
  "standardized_signflip"
)

PERTURB_RESAMPLING_METHODS <- PERTURB_METHOD_LEVELS[-1L]

PERTURB_METHOD_DISPLAY <- c(
  "parametric_score" = "Parametric score",
  "information_studentized_permutation" =
    "Barry et al.",
  "studentized_permutation" = "Empirical-studentized permutation",
  "hgf_signflip_mc" = "Hemerik et al.",
  "standardized_signflip" = "De Santis et al."
)

PERTURB_METHOD_COLORS <- c(
  "parametric_score" = "#F8766D",
  "information_studentized_permutation" = "#C77CFF",
  "studentized_permutation" = "#619CFF",
  "hgf_signflip_mc" = "#00BA38",
  "standardized_signflip" = "#00BFC4"
)

if (!identical(names(PERTURB_METHOD_DISPLAY), PERTURB_METHOD_LEVELS) ||
    !identical(names(PERTURB_METHOD_COLORS), PERTURB_METHOD_LEVELS) ||
    !identical(PERTURB_RESAMPLING_METHODS, PERTURB_METHOD_LEVELS[-1L])) {
  stop("Perturb-seq method configuration is internally inconsistent")
}

perturb_q_key <- function(q) {
  vapply(as.numeric(q), function(value) format(value, scientific = FALSE, trim = TRUE), character(1))
}

## Construct and validate the response-gene pool used for synthetic nulls.
## The preprocessing object historically contained the full analysis panel,
## including the positive controls, so the set difference is deliberately
## repeated here as a defensive check.
perturb_null_gene_pool <- function(
    gene_list, negative_gene_ids, positive_targets, expected_size = NULL) {
  gene_list <- as.character(gene_list)
  negative_gene_ids <- as.character(negative_gene_ids)
  positive_targets <- as.character(positive_targets)

  if (anyNA(gene_list) || anyNA(negative_gene_ids) ||
      anyNA(positive_targets) || any(gene_list == "") ||
      any(negative_gene_ids == "") || any(positive_targets == "")) {
    stop("gene identifiers must be non-missing, non-empty strings")
  }
  if (anyDuplicated(gene_list)) {
    stop("gene_list must contain unique gene identifiers")
  }

  pool <- unique(setdiff(negative_gene_ids, positive_targets))
  if (length(pool) == 0L) {
    stop("the non-positive null-gene pool is empty")
  }
  if (!all(pool %in% gene_list)) {
    stop("the null-gene pool contains genes outside gene_list")
  }
  if (any(pool %in% positive_targets)) {
    stop("the null-gene pool overlaps the positive-control targets")
  }
  if (!is.null(expected_size) && length(pool) != as.integer(expected_size)) {
    stop(
      "expected ", as.integer(expected_size),
      " non-positive null genes but found ", length(pool)
    )
  }
  pool
}

## Fixed reference layout: four blocks per replicate, using coordinate slot 2.
## This preserves the selected experiment's exact base-R draws independently
## of the number of reported q values. Stream 1 pairs the permutation methods;
## stream 2 pairs the sign-flip methods. The maximum full-run seed is 2e9.
perturb_resample_seed <- function(
    rep_index, max_reps, q, q_levels, iteration, pair_idx, n_pairs,
    B_max, stream) {
  if (!identical(as.numeric(q), PERTURB_Q_LEVELS) ||
      !identical(as.numeric(q_levels), PERTURB_Q_LEVELS)) stop("Only q = 0.10 is supported")
  if (length(rep_index) != 1L || rep_index != as.integer(rep_index) ||
      rep_index < 1L || rep_index > max_reps) {
    stop("replicate index is outside 1:max_reps")
  }
  if (length(iteration) != 1L || iteration < 1L || iteration > B_max) {
    stop("resampling iteration is outside 1:B_max")
  }
  if (length(pair_idx) != 1L || pair_idx < 1L || pair_idx > n_pairs) {
    stop("pair index is outside 1:n_pairs")
  }
  if (length(stream) != 1L || !stream %in% c(1L, 2L)) {
    stop("resampling stream must be 1 (permutation) or 2 (sign flip)")
  }

  block_size <- as.double(B_max) * n_pairs
  block_index <- 2L + (stream - 1L)
  coordinate <- (as.double(iteration) - 1) * n_pairs + (pair_idx - 1L)
  blocks_per_rep <- 4L
  total_span <- as.double(max_reps) * blocks_per_rep * block_size
  if (total_span > 2147483646) {
    stop("configured resampling grid exceeds the collision-free seed space")
  }
  seed <- (as.double(rep_index) - 1) * blocks_per_rep * block_size +
    block_index * block_size + coordinate + 1
  as.integer(seed)
}

## One canonical configuration shared by native runs and verified adoption.
perturb_run_config <- function(rep_idx, preprocessing_contract, n_nt_exposure_guides,
                               median_size, null_gene_pool) {
  list(schema_version = PERTURB_SCHEMA_VERSION, method_levels = PERTURB_METHOD_LEVELS,
       n_reps_max = 20L, rep_idx = as.integer(rep_idx), seed = 1000L + as.integer(rep_idx),
       n_pos = 50L, n_null = 950L, n_total = 1000L, median_size = as.integer(median_size),
       n_cells_full = preprocessing_contract$counts$retained_cells,
       n_pure_nt = preprocessing_contract$counts$pure_nt_cells,
       preprocessing_policy = PERTURB_COHORT_POLICY, preprocessing_contract = preprocessing_contract,
       covariate_names = preprocessing_contract$covariates,
       n_nt_exposure_guides = as.integer(n_nt_exposure_guides),
       null_gene_policy = PERTURB_NULL_GENE_POLICY, null_gene_pool_size = length(null_gene_pool),
       null_gene_pool = as.character(null_gene_pool), q_levels = PERTURB_Q_LEVELS,
       h = 15L, B_max = 25000L, alternative = "two.sided",
       rng_layout = PERTURB_RNG_LAYOUT_ID, RNGkind = PERTURB_RNG_KIND)
}
