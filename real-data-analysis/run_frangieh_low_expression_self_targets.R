#!/usr/bin/env Rscript
## run_frangieh_low_expression_self_targets.R
##
## Exploratory low-expression self-target sensitivity analysis for the
## Frangieh et al. (2021) control condition. This is deliberately separate
## from the production 50-highest-expression FDR workflow: it derives all
## eligible self-targets, ranks them by raw expression in pure non-targeting
## cells, and selects the lowest N without using an X--Y association. The
## original file-backed ODM is supported, as is the scPerturb H5AD release.
## The latter contains log1p(CPM), which is inverted to integer UMI counts
## after an explicitly validated alignment to the saved Frangieh metadata.
## Barry et al.'s information-studentized permutation score and the
## empirically studentized score share the NB null fit and every label
## permutation. The primary test is two-sided, based on |T|.

suppressPackageStartupMessages({
  library(ondisc)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(script_dir, "score_helpers.R"))

usage <- function() {
  paste(
    "Usage: Rscript --vanilla run_frangieh_low_expression_self_targets.R",
    "  [--h5ad FILE] [--intermediate-dir DIR] [--output-dir DIR]",
    "  [--n-targets 50 | --all-targets]",
    "  [--genes GENE1,GENE2,...]  # optional subset of selected panel",
    "  [--min-nt-nonzero 50] [--min-treated 5] [--max-treated N]",
    "  [--permutations 499] [--seed 20260822] [--cores 1]",
    sep = "\n"
  )
}

parse_cli <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    if (token %in% c("--help", "--all-targets")) {
      out[[gsub("-", "_", sub("^--", "", token), fixed = TRUE)]] <- TRUE
      i <- i + 1L
      next
    }
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
      key <- gsub("-", "_", bits[[1L]], fixed = TRUE)
      out[[key]] <- paste(bits[-1L], collapse = "=")
      i <- i + 1L
      next
    }
    if (i == length(x)) stop("Missing value after ", token)
    key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

as_count <- function(x, label, minimum = 1L) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(label, " must be an integer >= ", minimum)
  }
  value
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(cli$help)) {
  cat(usage(), "\n")
  quit(save = "no", status = 0L)
}

h5ad_mode <- !is.null(cli$h5ad)
if (h5ad_mode) {
  if (!file.exists(cli$h5ad)) stop("H5AD file not found: ", cli$h5ad)
  if (!requireNamespace("HDF5Array", quietly = TRUE) ||
      !requireNamespace("rhdf5", quietly = TRUE)) {
    stop("--h5ad requires the HDF5Array and rhdf5 packages")
  }
  h5ad_file <- normalizePath(cli$h5ad)
} else {
  h5ad_file <- NA_character_
}

results_root <- Sys.getenv(
  "LOCAL_PDRT_DATA_DIR", unset = file.path(project_root, "results")
)
intermediate_dir <- if (!is.null(cli$intermediate_dir)) {
  cli$intermediate_dir
} else {
  file.path(results_root, "frangieh-2021", "intermediate-data")
}
required_inputs <- c(
  "subsetted_ODM.rds", "covariate_matrix.rds", "pure_nt.rds"
)
missing_inputs <- required_inputs[
  !file.exists(file.path(intermediate_dir, required_inputs))
]
if (length(missing_inputs)) {
  stop(
    "Missing Frangieh intermediate file(s) in ", intermediate_dir, ": ",
    paste(missing_inputs, collapse = ", ")
  )
}
intermediate_dir <- normalizePath(intermediate_dir)

output_dir <- if (!is.null(cli$output_dir)) {
  cli$output_dir
} else {
  file.path(results_root, "frangieh-2021",
            "low-expression-self-target-score-pilot")
}
if (dir.exists(output_dir) && length(list.files(output_dir))) {
  warning("Output directory is nonempty; named outputs will be replaced: ",
          output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)

n_targets_requested <- as_count(
  if (is.null(cli$n_targets)) 50 else cli$n_targets, "--n-targets"
)
select_all_targets <- isTRUE(cli$all_targets)
if (select_all_targets && !is.null(cli$n_targets)) {
  stop("Use only one of --n-targets and --all-targets")
}
min_nt_nonzero <- as_count(
  if (is.null(cli$min_nt_nonzero)) 50 else cli$min_nt_nonzero,
  "--min-nt-nonzero"
)
min_treated <- as_count(
  if (is.null(cli$min_treated)) 5 else cli$min_treated, "--min-treated"
)
max_treated <- if (is.null(cli$max_treated)) {
  Inf
} else {
  as_count(cli$max_treated, "--max-treated")
}
if (max_treated < min_treated) {
  stop("--max-treated must be at least --min-treated")
}
B <- as_count(
  if (is.null(cli$permutations)) 499 else cli$permutations,
  "--permutations", minimum = 19L
)
seed <- as_count(if (is.null(cli$seed)) 20260822 else cli$seed, "--seed")
cores <- as_count(if (is.null(cli$cores)) 1 else cli$cores, "--cores")

PERMUTATION_BATCH_SIZE <- 128L
ALPHA_LEVELS <- c(0.01, 0.05, 0.10)
MIN_TREATED_BATCHES <- 2L
MAX_TREATED_BATCH_SHARE <- 0.80

cat("Frangieh low-expression self-target score analysis\n")
cat("  backend:      ", if (h5ad_mode) "scPerturb H5AD" else "ondisc ODM",
    "\n", sep = "")
if (h5ad_mode) cat("  H5AD:        ", h5ad_file, "\n", sep = "")
cat("  intermediate: ", intermediate_dir, "\n", sep = "")
cat("  output:       ", output_dir, "\n", sep = "")
cat("  selection:    ", if (select_all_targets) {
  "all eligible targets"
} else {
  paste(n_targets_requested, "lowest pure-NT-expression eligible targets")
}, "\n", sep = "")
cat("  treated range:", min_treated, " to ",
    if (is.finite(max_treated)) max_treated else "Inf",
    if (h5ad_mode) " exact-one-guide cells\n" else " targeting cells\n",
    sep = " ")
cat("  permutations: ", B, " paired global draws\n", sep = "")

## ---- Load and validate the production intermediate schema ----------------
modalities <- readRDS(file.path(intermediate_dir, "subsetted_ODM.rds"))
covariate_matrix <- as.matrix(readRDS(
  file.path(intermediate_dir, "covariate_matrix.rds")
))
pure_nt <- as.logical(readRDS(file.path(intermediate_dir, "pure_nt.rds")))

if (!all(c("gene", "grna") %in% names(modalities))) {
  stop("subsetted_ODM.rds must contain gene and grna modalities")
}
n_cells <- ncol(modalities$gene)
if (ncol(modalities$grna) != n_cells ||
    nrow(covariate_matrix) != n_cells || length(pure_nt) != n_cells) {
  stop("Frangieh modality, covariate, and pure-NT cell dimensions disagree")
}
if (sum(pure_nt) < 100L) stop("Too few pure non-targeting cells")
if (any(!is.finite(covariate_matrix))) {
  stop("The saved covariate matrix contains non-finite values")
}

gene_features <- ondisc::get_feature_covariates(modalities$gene)
grna_features <- ondisc::get_feature_covariates(modalities$grna)
gene_cell_features <- ondisc::get_cell_covariates(modalities$gene)
grna_cell_features <- ondisc::get_cell_covariates(modalities$grna)
if (!all(c("target", "target_type") %in% names(grna_features))) {
  stop("The gRNA feature metadata lacks target/target_type")
}
if (!all(c("n_umis", "batch") %in% names(gene_cell_features)) ||
    !"n_nonzero" %in% names(grna_cell_features)) {
  stop("The saved Frangieh cell metadata lacks UMI, batch, or MOI fields")
}
cell_umi <- as.numeric(gene_cell_features$n_umis)
cell_batch <- as.character(gene_cell_features$batch)
cell_moi <- as.integer(grna_cell_features$n_nonzero)
if (length(cell_umi) != n_cells || length(cell_batch) != n_cells ||
    length(cell_moi) != n_cells || any(!is.finite(cell_umi)) ||
    any(cell_umi <= 0) || anyNA(cell_batch) || any(!nzchar(cell_batch)) ||
    anyNA(cell_moi) || any(cell_moi < 0L)) {
  stop("Invalid saved Frangieh cell metadata")
}
if (!"log_n_umis" %in% colnames(covariate_matrix) ||
    max(abs(covariate_matrix[, "log_n_umis"] - log1p(cell_umi))) > 1e-10) {
  stop("Saved covariates do not contain the expected log1p(UMI_count)")
}
gene_ids <- rownames(gene_features)
grna_features$grna_id <- rownames(grna_features)
targeting_grnas <- grna_features[
  grna_features$target_type == "gene" &
    grna_features$target %in% gene_ids,
  c("grna_id", "target"), drop = FALSE
]
targeting_grnas <- targeting_grnas[
  !is.na(targeting_grnas$target) & nzchar(targeting_grnas$target),
  , drop = FALSE
]
target_genes <- sort(unique(targeting_grnas$target))
if (!length(target_genes)) stop("No self-target response genes are available")
targeting_ids <- targeting_grnas$grna_id
grna_rows_by_target <- split(seq_along(targeting_ids),
                             targeting_grnas$target)

recover_h5ad_counts <- function(log_cpm, depths, context,
                                tolerance = 5e-3) {
  log_cpm <- as.matrix(log_cpm)
  depths <- as.numeric(depths)
  if (ncol(log_cpm) != length(depths) || any(!is.finite(log_cpm)) ||
      any(log_cpm < 0) || any(!is.finite(depths)) || any(depths <= 0)) {
    stop("Invalid log1p(CPM) block or UMI depth for ", context)
  }
  unrounded <- sweep(expm1(log_cpm), 2L, depths / 1e6, `*`)
  counts <- round(unrounded)
  max_error <- if (length(unrounded)) {
    max(abs(unrounded - counts))
  } else {
    0
  }
  if (!is.finite(max_error) || max_error > tolerance || any(counts < 0)) {
    stop(
      "log1p(CPM) inversion failed for ", context,
      sprintf(" (maximum rounding residual %.6g; tolerance %.6g)",
              max_error, tolerance)
    )
  }
  storage.mode(counts) <- "double"
  attr(counts, "max_rounding_error") <- max_error
  counts
}

if (h5ad_mode) {
  ## This scPerturb file uses the legacy AnnData categorical layout: integer
  ## codes in /obs/<field>, with strings in /obs/__categories/<field>.
  h5_listing <- rhdf5::h5ls(h5ad_file, recursive = TRUE)
  dataset_paths <- paste0(
    ifelse(h5_listing$group == "/", "", h5_listing$group),
    "/", h5_listing$name
  )[h5_listing$otype == "H5I_DATASET"]
  read_obs <- function(key, categorical = FALSE) {
    path <- paste0("/obs/", key)
    if (!path %in% dataset_paths) stop("Missing H5AD field: ", path)
    value <- rhdf5::h5read(h5ad_file, path)
    if (!categorical) return(value)
    categories_path <- paste0("/obs/__categories/", key)
    if (!categories_path %in% dataset_paths) {
      stop("Missing H5AD categorical levels: ", categories_path)
    }
    codes <- as.integer(value)
    categories <- as.character(rhdf5::h5read(h5ad_file, categories_path))
    ans <- rep(NA_character_, length(codes))
    ok <- !is.na(codes) & codes >= 0L
    if (any(codes[ok] >= length(categories))) {
      stop("Out-of-range H5AD categorical code for ", key)
    }
    ans[ok] <- categories[codes[ok] + 1L]
    ans
  }

  h5_condition <- as.character(read_obs("condition", categorical = TRUE))
  h5_umi <- suppressWarnings(as.numeric(read_obs("UMI_count")))
  h5_moi <- suppressWarnings(as.integer(read_obs("MOI")))
  h5_guide <- as.character(read_obs("sgRNA", categorical = TRUE))
  var_gene <- as.character(rhdf5::h5read(h5ad_file, "/var/gene"))
  X_delayed <- suppressWarnings(HDF5Array::H5ADMatrix(h5ad_file))
  if (!identical(as.integer(dim(X_delayed)),
                 c(length(var_gene), length(h5_condition)))) {
    stop("Unexpected H5AD matrix orientation or dimensions")
  }
  if (anyDuplicated(var_gene)) stop("H5AD var/gene must be unique")

  control_global_idx <- which(h5_condition == "Control")
  alignment_drop_positions <- c(27337L, 34511L, 56031L)
  if (length(control_global_idx) != n_cells + length(alignment_drop_positions) ||
      any(alignment_drop_positions > length(control_global_idx))) {
    stop("Unexpected Control-cell count; cannot apply audited alignment")
  }
  aligned_control_global_idx <- control_global_idx[-alignment_drop_positions]
  aligned_h5_umi <- h5_umi[aligned_control_global_idx]
  aligned_h5_moi <- h5_moi[aligned_control_global_idx]
  aligned_h5_guide <- h5_guide[aligned_control_global_idx]
  if (!identical(aligned_h5_umi, cell_umi)) {
    stop("H5AD/local UMI vectors do not align after the three audited drops")
  }
  if (!identical(aligned_h5_moi, cell_moi)) {
    stop("H5AD/local MOI vectors do not align after the three audited drops")
  }

  h5_guide_categories <- as.character(rhdf5::h5read(
    h5ad_file, "/obs/__categories/sgRNA"
  ))
  if (length(h5_guide_categories) != nrow(grna_features) ||
      !setequal(h5_guide_categories, rownames(grna_features))) {
    stop("The 818 H5AD sgRNA categories do not match local gRNA metadata")
  }
  exact_one_guide <- aligned_h5_moi == 1L
  if (any(is.na(aligned_h5_guide[exact_one_guide])) ||
      any(!aligned_h5_guide[exact_one_guide] %in% rownames(grna_features))) {
    stop("An MOI==1 cell lacks a recognized H5AD sgRNA assignment")
  }
  guide_type <- rep(NA_character_, n_cells)
  guide_target <- rep(NA_character_, n_cells)
  assigned <- exact_one_guide & !is.na(aligned_h5_guide)
  guide_type[assigned] <- grna_features[
    aligned_h5_guide[assigned], "target_type"
  ]
  guide_target[assigned] <- grna_features[
    aligned_h5_guide[assigned], "target"
  ]
  pure_nt <- exact_one_guide & guide_type == "non-targeting"
  is_targeting <- exact_one_guide & guide_type == "gene"
  pure_nt[is.na(pure_nt)] <- FALSE
  is_targeting[is.na(is_targeting)] <- FALSE
  if (sum(pure_nt) != 3770L || sum(is_targeting) != 26716L) {
    stop(
      "Unexpected exact-one-guide classification: pure NT=", sum(pure_nt),
      ", targeting=", sum(is_targeting), " (expected 3770 and 26716)"
    )
  }
  X_by_target <- lapply(target_genes, function(g) {
    as.integer(is_targeting & !is.na(guide_target) & guide_target == g)
  })
  names(X_by_target) <- target_genes

  response_index <- match(target_genes, var_gene)
  if (anyNA(response_index)) {
    stop("Some target response genes are absent from H5AD var/gene")
  }
  pure_nt_idx <- which(pure_nt)
  Y_nt <- recover_h5ad_counts(
    X_delayed[response_index, aligned_control_global_idx[pure_nt_idx],
              drop = FALSE],
    cell_umi[pure_nt_idx], "pure-NT target responses"
  )

  ## A whole-transcriptome audit validates both integer reconstruction and
  ## column totals, rather than checking only the response genes used below.
  audit_local_idx <- unique(as.integer(round(seq(
    1, n_cells, length.out = 10L
  ))))
  audit_counts <- recover_h5ad_counts(
    X_delayed[, aligned_control_global_idx[audit_local_idx], drop = FALSE],
    cell_umi[audit_local_idx], "10-cell whole-transcriptome audit"
  )
  count_inversion_max_error <- attr(audit_counts, "max_rounding_error")
  if (!identical(as.numeric(colSums(audit_counts)),
                 cell_umi[audit_local_idx])) {
    stop("Reconstructed whole-transcriptome UMI sums do not match metadata")
  }
  rm(audit_counts)
} else {
  pure_nt_idx <- which(pure_nt)
  Y_nt <- as.matrix(modalities$gene[[target_genes, pure_nt_idx]])
  if (!identical(dim(Y_nt), c(length(target_genes),
                              length(pure_nt_idx)))) {
    stop("Unexpected pure-NT self-target matrix dimensions")
  }
  if (any(!is.finite(Y_nt)) || any(Y_nt < 0) ||
      any(abs(Y_nt - round(Y_nt)) > 1e-8)) {
    stop("The self-target response block is not raw UMI counts")
  }
  storage.mode(Y_nt) <- "double"

  G_target <- modalities$grna[[targeting_ids, ]]
  if (nrow(G_target) != length(targeting_ids) || ncol(G_target) != n_cells) {
    stop("Unexpected targeting-gRNA matrix dimensions")
  }
  if (length(G_target@x) &&
      (any(!is.finite(G_target@x)) || any(G_target@x < 0))) {
    stop("The targeting-gRNA matrix contains invalid assignments")
  }
  X_by_target <- lapply(target_genes, function(g) {
    rows <- grna_rows_by_target[[g]]
    as.integer(Matrix::colSums(G_target[rows, , drop = FALSE]) >= 1L)
  })
  names(X_by_target) <- target_genes
  aligned_control_global_idx <- NULL
  response_index <- match(target_genes, gene_ids)
  count_inversion_max_error <- NA_real_
}

target_table <- data.frame(
  gene = target_genes,
  n_grnas = vapply(target_genes, function(g) {
    length(grna_rows_by_target[[g]])
  }, integer(1)),
  n_treated = vapply(X_by_target, sum, integer(1)),
  n_treated_batches = vapply(X_by_target, function(x) {
    length(unique(cell_batch[x == 1L]))
  }, integer(1)),
  max_treated_batch_share = vapply(X_by_target, function(x) {
    idx <- which(x == 1L)
    if (!length(idx)) return(NA_real_)
    max(table(cell_batch[idx])) / length(idx)
  }, numeric(1)),
  mean_count_pure_nt = rowMeans(Y_nt),
  nonzero_pure_nt = rowSums(Y_nt != 0),
  stringsAsFactors = FALSE
)
target_table$eligible <- with(
  target_table,
  n_treated >= min_treated & n_treated <= max_treated &
    n_treated_batches >= MIN_TREATED_BATCHES &
    max_treated_batch_share <= MAX_TREATED_BATCH_SHARE &
    nonzero_pure_nt >= min_nt_nonzero
)
eligible_idx <- which(target_table$eligible)
if (!length(eligible_idx)) stop("No eligible low-expression self-targets")
selection_order <- eligible_idx[order(
  target_table$mean_count_pure_nt[eligible_idx],
  target_table$gene[eligible_idx]
)]
selected_idx <- if (select_all_targets) {
  selection_order
} else {
  head(selection_order, min(n_targets_requested, length(selection_order)))
}
panel_selected_idx <- selected_idx
requested_genes <- if (is.null(cli$genes)) {
  character()
} else {
  unique(trimws(strsplit(cli$genes, ",", fixed = TRUE)[[1L]]))
}
requested_genes <- requested_genes[nzchar(requested_genes)]
if (length(requested_genes)) {
  missing_from_panel <- setdiff(
    requested_genes, target_table$gene[panel_selected_idx]
  )
  if (length(missing_from_panel)) {
    stop(
      "--genes must be a subset of the prespecified selected panel; missing: ",
      paste(missing_from_panel, collapse = ", ")
    )
  }
  selected_idx <- panel_selected_idx[
    target_table$gene[panel_selected_idx] %in% requested_genes
  ]
}
selected_seed_index <- match(selected_idx, selection_order)
target_table$selected <- seq_len(nrow(target_table)) %in% selected_idx
selected_table <- target_table[selected_idx, , drop = FALSE]
selected_genes <- selected_table$gene

cat("  cells:        ", n_cells, " total; ", sum(pure_nt),
    " pure NTC\n", sep = "")
if (h5ad_mode) {
  cat("  H5AD mode:    26716 exact-one-guide targeting cells; ",
      sprintf("count audit max error %.3g", count_inversion_max_error),
      "\n", sep = "")
}
cat("  targets:      ", nrow(target_table), " present; ",
    length(eligible_idx), " eligible; ", nrow(selected_table),
    " selected\n", sep = "")
cat("  selected NTC mean range: ",
    paste(format(range(selected_table$mean_count_pure_nt), digits = 5),
          collapse = " to "), "\n", sep = "")

## Materialize every selected response only once. In H5AD mode, columns map
## through the audited Control-cell alignment; in ODM mode they are local
## Frangieh cell indices. Either route produces the same integer-count schema.
selected_treated_idx <- sort(unique(unlist(lapply(
  X_by_target[selected_genes], function(x) which(x == 1L)
), use.names = FALSE)))
union_idx <- sort(unique(c(pure_nt_idx, selected_treated_idx)))
union_position <- integer(n_cells)
union_position[union_idx] <- seq_along(union_idx)
if (h5ad_mode) {
  Y_selected_union <- recover_h5ad_counts(
    X_delayed[response_index[selected_idx],
              aligned_control_global_idx[union_idx], drop = FALSE],
    cell_umi[union_idx], "selected response/cell union"
  )
} else {
  Y_selected_union <- as.matrix(
    modalities$gene[[selected_genes, union_idx]]
  )
  if (any(!is.finite(Y_selected_union)) || any(Y_selected_union < 0) ||
      any(abs(Y_selected_union - round(Y_selected_union)) > 1e-8)) {
    stop("Selected ODM response/cell union is not raw UMI counts")
  }
  storage.mode(Y_selected_union) <- "double"
}
if (!identical(dim(Y_selected_union),
               c(length(selected_genes), length(union_idx)))) {
  stop("Unexpected selected response/cell union dimensions")
}

seed_for_pair <- function(base_seed, pair_index) {
  modulus <- 2147483646
  as.integer((as.double(base_seed) + 10007 * pair_index) %% modulus + 1)
}

draw_global_label_matrix <- function(n, n_exposed, n_draws) {
  ans <- matrix(0, nrow = n, ncol = n_draws)
  for (b in seq_len(n_draws)) {
    ans[sample.int(n, n_exposed, replace = FALSE), b] <- 1
  }
  ans
}

efficient_score_batch <- function(X_matrix, cache) {
  weighted_X <- cache$w_hat * X_matrix
  projection_coef <- cache$ZWZ_inv %*%
    crossprod(cache$Z_design, weighted_X)
  X_tilde <- X_matrix - cache$Z_design %*% projection_coef
  score_contrib <- cache$r_hat * X_tilde
  score_sum <- colSums(score_contrib)
  information_denom_sq <- colSums(cache$w_hat * X_tilde^2)
  centered <- score_contrib - rep(colMeans(score_contrib),
                                  each = nrow(score_contrib))
  empirical_denom_sq <- colSums(centered^2)
  if (any(!is.finite(information_denom_sq)) ||
      any(!is.finite(empirical_denom_sq)) ||
      any(information_denom_sq <= 0) || any(empirical_denom_sq <= 0)) {
    stop("Invalid batched permutation-score variance")
  }
  list(
    information_stat = score_sum / sqrt(information_denom_sq),
    empirical_stat = score_sum / sqrt(empirical_denom_sq)
  )
}

run_one_target <- function(j) {
  gene <- selected_genes[[j]]
  X_full <- X_by_target[[gene]]
  cell_keep <- pure_nt | X_full == 1L
  cell_idx <- which(cell_keep)
  X <- X_full[cell_idx]
  Y <- as.numeric(Y_selected_union[j, union_position[cell_idx]])
  if (any(!is.finite(Y)) || any(Y < 0) ||
      any(abs(Y - round(Y)) > 1e-8)) {
    stop("Response is not raw UMI counts for target ", gene)
  }
  Z_design <- cbind(
    intercept = 1,
    covariate_matrix[cell_idx, , drop = FALSE]
  )
  if (qr(Z_design)$rank != ncol(Z_design)) {
    stop("Rank-deficient covariate design for target ", gene)
  }

  cache <- precompute_score_cache(Y, Z_design)
  observed <- precompute_score_pair(X, cache)
  if (isTRUE(cache$degenerate) || isTRUE(observed$degenerate)) {
    stop("Degenerate score fit for target ", gene)
  }

  hits_information_abs <- hits_empirical_abs <- 0L
  hits_information_negative <- hits_empirical_negative <- 0L
  hits_information_positive <- hits_empirical_positive <- 0L
  set.seed(seed_for_pair(seed, selected_seed_index[[j]]))
  for (batch_start in seq.int(1L, B, by = PERMUTATION_BATCH_SIZE)) {
    n_draws <- min(PERMUTATION_BATCH_SIZE, B - batch_start + 1L)
    X_perm <- draw_global_label_matrix(length(X), sum(X), n_draws)
    permuted <- efficient_score_batch(X_perm, cache)
    hits_information_abs <- hits_information_abs + sum(
      abs(permuted$information_stat) >= abs(observed$information_stat)
    )
    hits_empirical_abs <- hits_empirical_abs + sum(
      abs(permuted$empirical_stat) >= abs(observed$empirical_stat)
    )
    hits_information_negative <- hits_information_negative + sum(
      permuted$information_stat <= -abs(observed$information_stat)
    )
    hits_empirical_negative <- hits_empirical_negative + sum(
      permuted$empirical_stat <= -abs(observed$empirical_stat)
    )
    hits_information_positive <- hits_information_positive + sum(
      permuted$information_stat >= abs(observed$information_stat)
    )
    hits_empirical_positive <- hits_empirical_positive + sum(
      permuted$empirical_stat >= abs(observed$empirical_stat)
    )
  }

  nt_mean <- mean(Y[X == 0L])
  treated_mean <- mean(Y[X == 1L])
  row <- selected_table[j, , drop = FALSE]
  data.frame(
    pair_index = selected_seed_index[[j]],
    gene = gene,
    n_grnas = row$n_grnas,
    n_treated_batches = row$n_treated_batches,
    max_treated_batch_share = row$max_treated_batch_share,
    n_cells = length(X),
    n_pure_nt = sum(X == 0L),
    n_treated = sum(X),
    prevalence = mean(X),
    mean_count_pure_nt = nt_mean,
    mean_count_treated = treated_mean,
    log2_marginal_fold_change = log2(
      (treated_mean + 1e-8) / (nt_mean + 1e-8)
    ),
    fraction_nonzero_pure_nt = mean(Y[X == 0L] > 0),
    fraction_nonzero_treated = mean(Y[X == 1L] > 0),
    theta_hat = cache$theta_hat,
    dispersion_method = cache$dispersion_method,
    information_stat = observed$information_stat,
    empirical_stat = observed$empirical_stat,
    information_denom_sq = observed$information_denom_sq,
    empirical_denom_sq = observed$empirical_denom_sq,
    empirical_to_information_denom_ratio =
      observed$empirical_denom_sq / observed$information_denom_sq,
    score_eff_n = observed$score_eff_n,
    max_score_share = observed$max_score_share,
    hits_abs_barry = hits_information_abs,
    hits_abs_empirical = hits_empirical_abs,
    hits_negative_extreme_barry = hits_information_negative,
    hits_negative_extreme_empirical = hits_empirical_negative,
    hits_positive_extreme_barry = hits_information_positive,
    hits_positive_extreme_empirical = hits_empirical_positive,
    p_two_sided_barry = (1 + hits_information_abs) / (B + 1),
    p_two_sided_empirical = (1 + hits_empirical_abs) / (B + 1),
    stringsAsFactors = FALSE
  )
}

cat("Running ", length(selected_genes), " paired score tests on ", cores,
    " core(s)...\n", sep = "")
test_start <- Sys.time()
if (.Platform$OS.type != "windows" && cores > 1L) {
  result_list <- parallel::mclapply(
    seq_along(selected_genes), run_one_target,
    mc.cores = cores, mc.preschedule = TRUE
  )
} else {
  result_list <- lapply(seq_along(selected_genes), run_one_target)
}
if (any(vapply(result_list, inherits, logical(1), what = "try-error"))) {
  failed <- which(vapply(result_list, inherits, logical(1),
                         what = "try-error"))
  stop("Target(s) failed: ", paste(selected_genes[failed], collapse = ", "))
}
pair_results <- do.call(rbind, result_list)
test_seconds <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
cat(sprintf("  score tests complete in %.1fs\n", test_seconds))

summary_rows <- list()
paired_rows <- list()
for (a in seq_along(ALPHA_LEVELS)) {
  alpha <- ALPHA_LEVELS[[a]]
  reject_barry <- pair_results$p_two_sided_barry <= alpha
  reject_empirical <- pair_results$p_two_sided_empirical <= alpha
  summary_rows[[2L * a - 1L]] <- data.frame(
    method = "barry", method_display = "Barry et al.", alpha = alpha,
    n_tests = nrow(pair_results), n_detected = sum(reject_barry),
    detection_rate = mean(reject_barry)
  )
  summary_rows[[2L * a]] <- data.frame(
    method = "empirical",
    method_display = "Empirical-studentized permutation", alpha = alpha,
    n_tests = nrow(pair_results), n_detected = sum(reject_empirical),
    detection_rate = mean(reject_empirical)
  )
  empirical_only <- sum(reject_empirical & !reject_barry)
  barry_only <- sum(reject_barry & !reject_empirical)
  discordant <- empirical_only + barry_only
  paired_rows[[a]] <- data.frame(
    alpha = alpha,
    n_tests = nrow(pair_results),
    both = sum(reject_barry & reject_empirical),
    barry_only = barry_only,
    empirical_only = empirical_only,
    neither = sum(!reject_barry & !reject_empirical),
    empirical_minus_barry = mean(reject_empirical) - mean(reject_barry),
    exact_discordance_p = if (discordant > 0L) {
      stats::binom.test(empirical_only, discordant, p = 0.5)$p.value
    } else 1
  )
}
detection_summary <- do.call(rbind, summary_rows)
paired_summary <- do.call(rbind, paired_rows)

utils::write.csv(
  target_table,
  file.path(output_dir, "frangieh-low-expression-self-target-targets.csv"),
  row.names = FALSE
)
utils::write.csv(
  pair_results,
  file.path(output_dir, "frangieh-low-expression-self-target-pairs.csv"),
  row.names = FALSE
)
utils::write.csv(
  detection_summary,
  file.path(output_dir, "frangieh-low-expression-self-target-summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_summary,
  file.path(output_dir,
            "frangieh-low-expression-self-target-paired-summary.csv"),
  row.names = FALSE
)

result_object <- list(
  config = list(
    backend = if (h5ad_mode) "scperturb_h5ad" else "ondisc_odm",
    analysis_population = if (h5ad_mode) {
      "exactly_one_assigned_guide"
    } else {
      "production_odm_assignments"
    },
    seed = seed,
    n_targets_requested = n_targets_requested,
    n_targets_selected = nrow(selected_table),
    target_gene_filter = requested_genes,
    all_targets = select_all_targets,
    min_nt_nonzero = min_nt_nonzero,
    min_treated = min_treated,
    max_treated = max_treated,
    min_treated_batches = MIN_TREATED_BATCHES,
    max_treated_batch_share = MAX_TREATED_BATCH_SHARE,
    permutations = B,
    permutation_batch_size = PERMUTATION_BATCH_SIZE,
    cores = cores,
    alpha_levels = ALPHA_LEVELS,
    selection_mode = if (select_all_targets) {
      "all_eligible_ordered_by_pure_nt_expression"
    } else {
      "lowest_pure_nt_expression"
    },
    selection_uses_x_y_association = FALSE,
    primary_test = "two-sided absolute-statistic global label permutation",
    timings_seconds = list(score_tests = test_seconds)
  ),
  data_checks = list(
    n_cells = n_cells,
    n_pure_nt = sum(pure_nt),
    n_gene_features = nrow(gene_features),
    n_grna_features = nrow(grna_features),
    n_self_targets_present = nrow(target_table),
    n_eligible_self_targets = length(eligible_idx),
    raw_integer_nonnegative = TRUE,
    h5ad_control_rows_before_drops = if (h5ad_mode) {
      length(control_global_idx)
    } else {
      NA_integer_
    },
    h5ad_alignment_drop_positions = if (h5ad_mode) {
      alignment_drop_positions
    } else {
      integer()
    },
    h5ad_exact_one_guide_targeting = if (h5ad_mode) {
      sum(is_targeting)
    } else {
      NA_integer_
    },
    h5ad_count_inversion_max_error_10_cell_audit =
      count_inversion_max_error
  ),
  provenance = list(
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    intermediate_dir = intermediate_dir,
    h5ad_file = if (h5ad_mode) h5ad_file else NA_character_,
    h5ad_size_bytes = if (h5ad_mode) file.info(h5ad_file)$size else NA_real_,
    runner_md5 = unname(tools::md5sum(file.path(
      script_dir, "run_frangieh_low_expression_self_targets.R"
    ))),
    score_helpers_md5 = unname(tools::md5sum(file.path(
      script_dir, "score_helpers.R"
    ))),
    package_versions = vapply(
      c("ondisc", "Matrix", "MASS",
        if (h5ad_mode) c("HDF5Array", "rhdf5") else character()),
      function(pkg) as.character(utils::packageVersion(pkg)), character(1)
    ),
    session_info = sessionInfo()
  ),
  target_table = target_table,
  pair_results = pair_results,
  detection_summary = detection_summary,
  paired_summary = paired_summary
)
saveRDS(
  result_object,
  file.path(output_dir, "frangieh-low-expression-self-target-results.rds"),
  compress = "xz"
)

cat("\nTwo-sided detection summary:\n")
print(detection_summary, row.names = FALSE)
cat("\nPaired summary:\n")
print(paired_summary, row.names = FALSE)
cat("\nSaved results to: ", output_dir, "\n", sep = "")
cat(if (h5ad_mode) {
  paste0(
    "This is exact-one-guide empirical self-target sensitivity, ",
    "not formal power.\n"
  )
} else {
  "This is empirical self-target sensitivity, not formal power.\n"
})
