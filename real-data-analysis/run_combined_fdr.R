#!/usr/bin/env Rscript
## run_combined_fdr.R
##
## Combined null + positive-control FDR experiment for Frangieh 2021.
##
## Each replicate contains 50 positive gene--gRNA pairs and 950 synthetic
## null pairs.  All five methods use one fitted negative-binomial efficient
## score and match the comparison in the manuscript Frangieh section:
##   * parametric_score
##   * information_studentized_permutation
##   * studentized_permutation (empirical studentization)
##   * hgf_signflip_mc (Monte Carlo, no SPA)
##   * standardized_signflip
##
## The four resampling methods use pcmt::adaptive_fdr_bh separately at each
## target q.  Information/empirical permutation draws are paired by
## (replicate, q, round, pair), as are HGF/standardized sign draws.
##
## Usage: Rscript run_combined_fdr.R <rep_idx>
## Output: results/score-calibration-nonpositive-nulls/
##         score_fdr_rep<rep_idx>.rds

suppressPackageStartupMessages({
  library(pcmt)
  library(ondisc)
  library(dplyr)
  library(MASS)
})

if (!"adaptive_fdr_bh" %in% getNamespaceExports("pcmt")) {
  stop(
    "This analysis requires a pcmt version exporting adaptive_fdr_bh. ",
    "Install the same development version used for the original Frangieh run."
  )
}

## --- Paths -----------------------------------------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
results_root <- Sys.getenv(
  "LOCAL_PDRT_DATA_DIR",
  unset = file.path(project_root, "results")
)
intermediate_dir <- Sys.getenv(
  "FRANGIEH_INTERMEDIATE_DIR",
  unset = file.path(results_root, "frangieh-2021", "intermediate-data")
)
source(file.path(script_dir, "method_config.R"))
source(file.path(script_dir, "score_helpers.R"))
results_dir <- Sys.getenv(
  "FRANGIEH_SCORE_RESULTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_RESULTS_SUBDIR)
)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

## --- Parse args ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript run_combined_fdr.R <rep_idx>")
rep_idx <- as.integer(args[[1L]])
N_REPS_MAX <- 20L
stopifnot(!is.na(rep_idx), rep_idx >= 1L, rep_idx <= N_REPS_MAX)
rep_seed <- 1000L + rep_idx

## --- Load preprocessed data ------------------------------------------------
cat("Loading intermediate data...\n")
subsetted_ODM <- readRDS(file.path(intermediate_dir, "subsetted_ODM.rds"))
covariate_matrix <- readRDS(file.path(intermediate_dir, "covariate_matrix.rds"))
gene_list <- readRDS(file.path(intermediate_dir, "gene_list.rds"))
positive_control_pairs <- readRDS(
  file.path(intermediate_dir, "positive_control_pairs.rds")
)
negative_control_pairs <- readRDS(
  file.path(intermediate_dir, "negative_control_pairs.rds")
)
pure_nt <- readRDS(file.path(intermediate_dir, "pure_nt.rds"))
pos_cell_masks <- readRDS(file.path(intermediate_dir, "pos_cell_masks.rds"))

n_cells_full <- nrow(covariate_matrix)
Z_full <- as.matrix(covariate_matrix)
Z_design_full <- cbind(intercept = 1, Z_full)
n_pure_nt <- sum(pure_nt)
pure_nt_idx <- which(pure_nt)
Z_design_pure_nt <- Z_design_full[pure_nt_idx, , drop = FALSE]
cat(sprintf("Pure-NT cells (universal control): %d / %d\n",
            n_pure_nt, n_cells_full))

## --- Experiment constants --------------------------------------------------
Q_LEVELS <- c(0.01, 0.05)
H_STOP <- 15L
B_MAX <- 25000L
pos_targets <- unique(positive_control_pairs$grna_target)
n_pos <- length(pos_targets)
n_total <- as.integer(round(n_pos / 0.05))
n_null <- n_total - n_pos
median_size <- as.integer(median(table(positive_control_pairs$grna_target)))
nt_grna_ids <- negative_control_pairs$grna$grna_id
null_gene_pool <- perturb_null_gene_pool(
  gene_list = gene_list,
  negative_gene_ids = negative_control_pairs$gene_id,
  positive_targets = pos_targets,
  expected_size = n_null
)

cat(sprintf("Pair counts: %d alternatives + %d nulls = %d total\n",
            n_pos, n_null, n_total))
cat(sprintf(
  "Eligible non-positive response genes for null pairs: %d\n",
  length(null_gene_pool)
))
cat(sprintf("Median targeting gRNAs per positive (= null NT subset): %d\n",
            median_size))

## --- Exposure vectors ------------------------------------------------------
cat("Building positive-control exposure vectors...\n")
pos_X_full <- vector("list", n_pos)
names(pos_X_full) <- pos_targets
for (g in pos_targets) {
  grna_ids <- positive_control_pairs |>
    filter(grna_target == !!g) |>
    pull(grna_id)
  presences <- subsetted_ODM$grna[[grna_ids, ]] |>
    Matrix::t() |>
    Matrix::as.matrix()
  pos_X_full[[g]] <- as.integer(rowSums(presences) >= 1L)
}

cat("Loading non-targeting gRNA assignments...\n")
nt_mat_full <- subsetted_ODM$grna[[nt_grna_ids, ]] |>
  Matrix::t() |>
  Matrix::as.matrix()
colnames(nt_mat_full) <- nt_grna_ids

## --- Gene-level NB caches --------------------------------------------------
score_cache_env <- new.env(hash = TRUE)
y_cache_env <- new.env(hash = TRUE)

get_Y <- function(gene_id) {
  key <- as.character(gene_id)
  if (!is.null(y_cache_env[[key]])) return(y_cache_env[[key]])
  Y <- as.numeric(subsetted_ODM$gene[[gene_id, ]])
  y_cache_env[[key]] <- Y
  Y
}

## Positive controls and nulls use different cell masks, so the cache key
## includes both the gene and mask kind.
get_score_cache <- function(gene_id, mask_key, mask_idx) {
  key <- paste0(mask_key, "::", gene_id)
  if (!is.null(score_cache_env[[key]])) return(score_cache_env[[key]])
  Y <- get_Y(gene_id)[mask_idx]
  Z_design <- if (identical(mask_key, "null")) {
    Z_design_pure_nt
  } else {
    Z_design_full[mask_idx, , drop = FALSE]
  }
  cache <- precompute_score_cache(Y, Z_design)
  score_cache_env[[key]] <- cache
  cache
}

cat("Pre-warming score caches for positive-control genes...\n")
t0 <- Sys.time()
for (g in pos_targets) {
  mask_idx <- which(pos_cell_masks[[g]])
  invisible(get_score_cache(g, paste0("pos:", g), mask_idx))
}
cat(sprintf("  %d positive-control caches in %.1fs\n", n_pos,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## --- Deterministically paired randomization draws --------------------------
## The seed is a pure function of the analysis coordinates.  Consequently the
## two permutation methods see the same label permutation and the two flip
## methods see the same Rademacher vector, even after their active sets diverge.
draw_permutation <- function(n, q, iteration, pair_idx) {
  set.seed(perturb_resample_seed(
    rep_index = rep_idx,
    max_reps = N_REPS_MAX,
    q = q,
    q_levels = Q_LEVELS,
    iteration = iteration,
    pair_idx = pair_idx,
    n_pairs = n_total,
    B_max = B_MAX,
    stream = 1L
  ))
  sample.int(n)
}

draw_signs <- function(n, q, iteration, pair_idx) {
  set.seed(perturb_resample_seed(
    rep_index = rep_idx,
    max_reps = N_REPS_MAX,
    q = q,
    q_levels = Q_LEVELS,
    iteration = iteration,
    pair_idx = pair_idx,
    n_pairs = n_total,
    B_max = B_MAX,
    stream = 2L
  ))
  rademacher_signs(n)
}

compact_adaptive_result <- function(x) {
  list(
    p_value = x$p_value,
    rejected = x$rejected,
    futility = x$futility,
    active = x$active,
    K = x$K,
    stop_round = x$stop_round,
    t_total = x$t_total,
    R_threshold = x$R_threshold
  )
}

## Post-hoc asymptotic-score fallback for degenerate pairs.  A pair whose
## efficient-score cache is singular has no trustworthy resampling reference,
## so if it could not be resampled away we substitute the asymptotic score
## p-value p0 = 2 * pnorm(|information_stat|) -- the same quantity the
## parametric method reports -- and re-decide its rejection at the adaptive
## run's realized BH threshold.  When no pair is degenerate (the usual case
## after gene resampling) `asym_idx` is empty and the result is unchanged.
apply_asymptotic_fallback <- function(res, asym_idx, p_asymptotic) {
  if (length(asym_idx) == 0L) return(res)
  res$p_value[asym_idx] <- p_asymptotic[asym_idx]
  res$active <- setdiff(res$active, asym_idx)
  res$futility <- setdiff(res$futility, asym_idx)
  reject_asym <- asym_idx[p_asymptotic[asym_idx] <= res$R_threshold]
  res$rejected <- sort(unique(c(setdiff(res$rejected, asym_idx), reject_asym)))
  res
}

## ---------------------------------------------------------------------------
## Per-replicate workhorse
## ---------------------------------------------------------------------------
run_rep <- function() {
  cat(sprintf("\n=== Rep %d (seed = %d) ===\n", rep_idx, rep_seed))
  set.seed(rep_seed)

  ## 1. Construct 50 positive + 950 null pairs.  Null response genes are
  ## sampled only from the prespecified non-positive pool.
  pairs <- vector("list", n_total)
  for (k in seq_len(n_pos)) {
    g <- pos_targets[[k]]
    mask_idx <- which(pos_cell_masks[[g]])
    pairs[[k]] <- list(
      pair_idx = k,
      type = "alt",
      gene = g,
      mask_key = paste0("pos:", g),
      mask_idx = mask_idx,
      X = pos_X_full[[g]][mask_idx],
      nt_subset = NA_character_
    )
  }
  ## Draw each null gene, resampling if its efficient-score cache is degenerate
  ## (singular weighted information on the pure-NT cells).  Such genes are ~0.1%
  ## of the pool, so a non-degenerate draw is found almost immediately; the cap
  ## bounds the rare pathological case, after which the pair falls back to the
  ## asymptotic score p-value (see use_asymptotic below).
  MAX_GENE_RESAMPLE <- 200L
  n_resampled <- 0L
  for (i in seq_len(n_null)) {
    g <- sample(null_gene_pool, 1L)
    tries <- 0L
    while (isTRUE(get_score_cache(g, "null", pure_nt_idx)$degenerate) &&
           tries < MAX_GENE_RESAMPLE) {
      g <- sample(null_gene_pool, 1L)
      tries <- tries + 1L
    }
    if (tries > 0L) n_resampled <- n_resampled + 1L
    nt_subset <- sample(nt_grna_ids, median_size, replace = FALSE)
    X_full <- as.integer(
      rowSums(nt_mat_full[, nt_subset, drop = FALSE]) >= 1L
    )
    k <- n_pos + i
    pairs[[k]] <- list(
      pair_idx = k,
      type = "null",
      gene = g,
      mask_key = "null",
      mask_idx = pure_nt_idx,
      X = X_full[pure_nt_idx],
      nt_subset = paste(nt_subset, collapse = ",")
    )
  }
  if (n_resampled > 0L) {
    cat(sprintf("Resampled %d null gene(s) with degenerate score caches.\n",
                n_resampled))
  }
  is_alt <- vapply(pairs, function(x) x$type == "alt", logical(1))
  null_genes <- vapply(pairs[!is_alt], `[[`, character(1), "gene")
  if (any(null_genes %in% pos_targets) ||
      !all(null_genes %in% null_gene_pool)) {
    stop("constructed null pairs violate the non-positive gene policy")
  }

  ## 2. One shared efficient-score cache per pair.
  cat("Building shared efficient-score caches...\n")
  t0 <- Sys.time()
  pair_caches <- vector("list", n_total)
  for (k in seq_len(n_total)) {
    pair <- pairs[[k]]
    gene_cache <- get_score_cache(
      pair$gene, pair$mask_key, pair$mask_idx
    )
    pair_score <- precompute_score_pair(pair$X, gene_cache)
    pair_caches[[k]] <- list(
      X = pair$X,
      n_cells = length(pair$X),
      gene_cache = gene_cache,
      pair_score = pair_score
    )
  }
  cache_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %d pair caches built in %.1fs\n", n_total, cache_seconds))

  ## Pairs whose score cache is still degenerate (resampling exhausted, or a
  ## degenerate positive control) use the asymptotic score p-value instead of a
  ## resampling reference.  Expected to be empty after null-gene resampling.
  use_asymptotic <- vapply(
    pair_caches, function(x) isTRUE(x$gene_cache$degenerate), logical(1)
  )
  asym_idx <- which(use_asymptotic)
  if (length(asym_idx) > 0L) {
    cat(sprintf(
      "  %d pair(s) degenerate after resampling; using asymptotic score p-value.\n",
      length(asym_idx)))
  }

  information_stat <- vapply(
    pair_caches,
    function(x) x$pair_score$information_stat,
    numeric(1)
  )
  empirical_stat <- vapply(
    pair_caches,
    function(x) x$pair_score$empirical_stat,
    numeric(1)
  )
  raw_score <- vapply(
    pair_caches,
    function(x) x$pair_score$score_sum,
    numeric(1)
  )

  observed <- list(
    parametric_score = information_stat,
    information_studentized_permutation = information_stat,
    studentized_permutation = empirical_stat,
    hgf_signflip_mc = raw_score,
    standardized_signflip = raw_score
  )
  stopifnot(identical(names(observed), PERTURB_METHOD_LEVELS))

  ## 3. Parametric normal calibration followed by ordinary BH.
  parametric_p <- 2 * pnorm(abs(information_stat), lower.tail = FALSE)
  parametric_rejected <- setNames(vector("list", length(Q_LEVELS)), Q_LEVELS)
  for (q in Q_LEVELS) {
    key <- as.character(q)
    parametric_rejected[[key]] <- which(
      p.adjust(parametric_p, method = "BH") <= q
    )
  }

  ## 4. Callback factories for the four adaptive randomization calibrations.
  make_permutation_callback <- function(q, empirical) {
    force(q)
    force(empirical)
    function(active_idx, t) {
      vapply(active_idx, function(j) {
        cache <- pair_caches[[j]]
        permutation <- draw_permutation(cache$n_cells, q, t, j)
        X_perm <- cache$X[permutation]
        if (empirical) {
          empirical_score_stat(X_perm, cache$gene_cache)
        } else {
          information_score_stat(X_perm, cache$gene_cache)
        }
      }, numeric(1))
    }
  }

  make_flip_callback <- function(q, standardized) {
    force(q)
    force(standardized)
    function(active_idx, t) {
      vapply(active_idx, function(j) {
        cache <- pair_caches[[j]]
        signs <- draw_signs(cache$n_cells, q, t, j)
        if (standardized) {
          standardized_signflip_stat(signs, cache$pair_score)
        } else {
          hgf_signflip_stat(signs, cache$pair_score)
        }
      }, numeric(1))
    }
  }

  methods <- lapply(PERTURB_METHOD_LEVELS, function(method) {
    list(observed = observed[[method]])
  })
  names(methods) <- PERTURB_METHOD_LEVELS
  methods$parametric_score$p_value <- parametric_p
  methods$parametric_score$rejected <- parametric_rejected

  timings <- list(cache = cache_seconds)
  for (method in PERTURB_RESAMPLING_METHODS) {
    methods[[method]]$by_q <- setNames(
      vector("list", length(Q_LEVELS)), Q_LEVELS
    )
    timings[[method]] <- setNames(numeric(length(Q_LEVELS)), Q_LEVELS)
  }

  for (q in Q_LEVELS) {
    key <- as.character(q)
    callbacks <- list(
      information_studentized_permutation =
        make_permutation_callback(q, empirical = FALSE),
      studentized_permutation =
        make_permutation_callback(q, empirical = TRUE),
      hgf_signflip_mc = make_flip_callback(q, standardized = FALSE),
      standardized_signflip = make_flip_callback(q, standardized = TRUE)
    )

    for (method in PERTURB_RESAMPLING_METHODS) {
      cat(sprintf("%s (q=%.2f, h=%d) ...\n",
                  PERTURB_METHOD_DISPLAY[[method]], q, H_STOP))
      t0 <- Sys.time()
      adaptive <- pcmt::adaptive_fdr_bh(
        T_orig = observed[[method]],
        permute_stat = callbacks[[method]],
        alpha = q,
        h = H_STOP,
        B_max = B_MAX,
        alternative = "two.sided"
      )
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      adaptive <- apply_asymptotic_fallback(adaptive, asym_idx, parametric_p)
      timings[[method]][[key]] <- elapsed
      methods[[method]]$by_q[[key]] <- compact_adaptive_result(adaptive)
      cat(sprintf(
        "  %.1fs; adaptive rounds=%d; rejected=%d (alt=%d, null=%d)\n",
        elapsed,
        adaptive$t_total,
        length(adaptive$rejected),
        sum(is_alt[adaptive$rejected]),
        sum(!is_alt[adaptive$rejected])
      ))
    }
  }

  pair_table <- data.frame(
    pair_idx = seq_len(n_total),
    gene = vapply(pairs, `[[`, character(1), "gene"),
    pair_type = ifelse(is_alt, "alt", "null"),
    is_alt = is_alt,
    nt_subset = vapply(pairs, `[[`, character(1), "nt_subset"),
    n_cells = vapply(pair_caches, `[[`, integer(1), "n_cells"),
    n_exposed = vapply(pair_caches, function(x) sum(x$X), integer(1)),
    theta_hat = vapply(
      pair_caches, function(x) x$gene_cache$theta_hat, numeric(1)
    ),
    score_sum = raw_score,
    information_denom_sq = vapply(
      pair_caches,
      function(x) x$pair_score$information_denom_sq,
      numeric(1)
    ),
    empirical_denom_sq = vapply(
      pair_caches,
      function(x) x$pair_score$empirical_denom_sq,
      numeric(1)
    ),
    information_stat = information_stat,
    empirical_stat = empirical_stat,
    score_eff_n = vapply(
      pair_caches, function(x) x$pair_score$score_eff_n, numeric(1)
    ),
    max_score_share = vapply(
      pair_caches, function(x) x$pair_score$max_score_share, numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    rep = rep_idx,
    seed = rep_seed,
    pairs = pair_table,
    methods = methods,
    timings = timings
  )
}

## --- Run and save ----------------------------------------------------------
t_global <- Sys.time()
result <- run_rep()
total_minutes <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
cat(sprintf("\n=== Rep %d completed in %.1f min ===\n", rep_idx, total_minutes))

output <- list(
  config = list(
    schema_version = PERTURB_SCHEMA_VERSION,
    method_levels = PERTURB_METHOD_LEVELS,
    n_reps_max = N_REPS_MAX,
    rep_idx = rep_idx,
    seed = rep_seed,
    n_pos = n_pos,
    n_null = n_null,
    n_total = n_total,
    median_size = median_size,
    n_cells_full = n_cells_full,
    n_pure_nt = n_pure_nt,
    null_gene_policy = PERTURB_NULL_GENE_POLICY,
    null_gene_pool_size = length(null_gene_pool),
    null_gene_pool = null_gene_pool,
    q_levels = Q_LEVELS,
    h = H_STOP,
    B_max = B_MAX,
    alternative = "two.sided"
  ),
  result = result
)

out_file <- file.path(
  results_dir,
  sprintf("score_fdr_rep%03d.rds", rep_idx)
)
saveRDS(output, out_file)
cat("Result saved to:", out_file, "\n")
