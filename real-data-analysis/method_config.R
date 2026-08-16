## Shared method and output-schema configuration for the Frangieh score study.
##
## Keep these values synchronized with the Frangieh section of the manuscript. Both
## the per-replicate runner and the downstream scripts source this file so a
## result cannot silently change method order or plotting colours.

PERTURB_SCHEMA_VERSION <- "perturb_score5_v2_nonpositive_nulls"
PERTURB_RESULTS_SUBDIR <- "score-calibration-nonpositive-nulls"
PERTURB_PLOTS_SUBDIR <- "score-calibration-nonpositive-nulls"
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
  format(as.numeric(q), scientific = FALSE, trim = TRUE)
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

## Map a resampling coordinate to a reproducible base-R seed.  Across the full
## 20-replicate experiment, every (replicate, q, stream, iteration, pair) tuple
## occupies a distinct position below R's seed modulus.  Stream 1 is shared by
## the two permutation methods and stream 2 by the two flip methods.
perturb_resample_seed <- function(
    rep_index, max_reps, q, q_levels, iteration, pair_idx, n_pairs,
    B_max, stream) {
  q_index <- match(perturb_q_key(q), perturb_q_key(q_levels))
  if (is.na(q_index)) stop("q is not in the configured q grid")
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
  block_index <- (q_index - 1L) * 2L + (stream - 1L)
  coordinate <- (as.double(iteration) - 1) * n_pairs + (pair_idx - 1L)
  blocks_per_rep <- length(q_levels) * 2L
  total_span <- as.double(max_reps) * blocks_per_rep * block_size
  if (total_span > 2147483646) {
    stop("configured resampling grid exceeds the collision-free seed space")
  }
  seed <- (as.double(rep_index) - 1) * blocks_per_rep * block_size +
    block_index * block_size + coordinate + 1
  as.integer(seed)
}
