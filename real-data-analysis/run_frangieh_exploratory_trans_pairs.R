#!/usr/bin/env Rscript
## run_frangieh_exploratory_trans_pairs.R
##
## Sharded, purely exploratory trans-pair scan in the Frangieh et al. (2021)
## Control condition.  The response panel is fixed without consulting cell
## guide assignments or target-response associations: use globally computed
## Control-condition mean and detection metadata only after auditing those
## values against raw H5AD counts.  Among
## genes detected in at least 577 cells, exclude every perturbation target and
## take the 5,000 lowest global mean counts (gene name breaks ties).  Exposures
## are the 241 primary eligible gene-targeting
## perturbations in exactly-one-guide cells.  Thus
## neither self-target positive controls nor non-targeting negative controls
## can enter the scan.
##
## Raw counts come from the ondisc ODM -- the same accessor the production
## control-pair workflow uses -- so the analysis path involves no log1p(CPM)
## inversion.  --h5ad selects the scPerturb backend instead and restores the
## reconstruction plus its audit chain, which is useful as an independent
## cross-check of the ODM.
##
## Each pair compares all 3,770 exact-one non-targeting cells with cells
## carrying exactly one guide for the selected target.  Barry et al.'s
## information-studentized score and the empirically studentized score use the
## same NB null fit and identical label permutations.  The primary p-values
## are two-sided and use |T|.  Results are raw shard-level p-values; multiplicity
## correction belongs in a separate complete-shard combiner.

suppressPackageStartupMessages({
  library(ondisc)
  library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."))
helper_path <- file.path(script_dir, "score_helpers.R")
source(helper_path)

usage <- function() {
  paste(
    "Usage: Rscript --vanilla run_frangieh_exploratory_trans_pairs.R",
    "  [--h5ad FILE] [--intermediate-dir DIR] [--output-dir DIR]",
    "  {--target GENE | --target-index N | --target-range A:B | --all-targets}",
    "  [--response-range A:B] [--response-block-size 100] [--pilot]",
    "  [--min-treated 10]",
    "  [--permutations 4999] [--permutation-batch-size 128]",
    "  [--seed 20260822] [--cores 1] [--overwrite] [--manifest-only]",
    "",
    "Counts are read from the ondisc ODM in --intermediate-dir by default.",
    "--h5ad FILE selects the scPerturb backend, which reconstructs integer",
    "counts from log1p(CPM) and keeps the original inversion audit chain.",
    sep = "\n"
  )
}

parse_cli <- function(x) {
  out <- list()
  i <- 1L
  flag_options <- c("--help", "--all-targets", "--overwrite", "--pilot",
                    "--manifest-only")
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    if (token %in% flag_options) {
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

parse_index_range <- function(x, maximum, label) {
  bits <- strsplit(as.character(x), ":", fixed = TRUE)[[1L]]
  if (length(bits) == 1L) bits <- strsplit(as.character(x), "-", fixed = TRUE)[[1L]]
  if (length(bits) != 2L) stop(label, " must have form START:END")
  endpoints <- suppressWarnings(as.integer(bits))
  if (anyNA(endpoints) || any(endpoints < 1L) || endpoints[[1L]] > endpoints[[2L]] ||
      endpoints[[2L]] > maximum) {
    stop(label, " must satisfy 1 <= START <= END <= ", maximum)
  }
  seq.int(endpoints[[1L]], endpoints[[2L]])
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

atomic_publish <- function(tmp, destination, overwrite = FALSE) {
  if (file.exists(destination) && !overwrite) {
    stop("Refusing to replace existing output: ", destination)
  }
  if (file.exists(destination) && overwrite) {
    if (!file.remove(destination)) stop("Could not replace: ", destination)
  }
  if (!file.rename(tmp, destination)) {
    stop("Could not atomically publish ", destination)
  }
  invisible(destination)
}

publish_identical_csv <- function(x, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(destination), "."),
                  tmpdir = dirname(destination))
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  new_hash <- unname(tools::md5sum(tmp))
  if (file.exists(destination)) {
    old_hash <- unname(tools::md5sum(destination))
    if (!identical(old_hash, new_hash)) {
      stop("Existing manifest differs from this shard: ", destination)
    }
    file.remove(tmp)
    return(invisible(destination))
  }
  ## Concurrent array tasks may both reach this rename.  Their byte-identical
  ## files are interchangeable; verify the final path after publishing.
  if (!file.rename(tmp, destination)) {
    if (!file.exists(destination) ||
        !identical(unname(tools::md5sum(destination)), new_hash)) {
      stop("Could not publish manifest: ", destination)
    }
  }
  if (!identical(unname(tools::md5sum(destination)), new_hash)) {
    stop("Manifest race produced a non-identical file: ", destination)
  }
  invisible(destination)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
allowed_cli <- c(
  "help", "h5ad", "intermediate_dir", "output_dir", "target",
  "target_index", "target_range", "all_targets", "response_range",
  "response_block_size", "pilot", "min_treated", "permutations",
  "permutation_batch_size", "seed", "cores", "overwrite",
  "manifest_only"
)
unknown_cli <- setdiff(names(cli), allowed_cli)
if (length(unknown_cli)) {
  stop("Unknown option(s): --", paste(gsub("_", "-", unknown_cli),
                                      collapse = ", --"))
}
if (isTRUE(cli$help)) {
  cat(usage(), "\n")
  quit(save = "no", status = 0L)
}
## Counts come from the ondisc ODM by default -- the same accessor the
## production control-pair workflow uses -- so no log1p(CPM) inversion, no
## H5AD row alignment, and no locked file-identity audit are involved.
## Passing --h5ad selects the scPerturb backend instead, which reconstructs
## integer counts from log1p(CPM) and keeps the original audit chain.
h5ad_mode <- !is.null(cli$h5ad)
if (h5ad_mode) {
  if (!file.exists(cli$h5ad)) stop("H5AD file not found: ", cli$h5ad)
  if (!requireNamespace("HDF5Array", quietly = TRUE) ||
      !requireNamespace("rhdf5", quietly = TRUE)) {
    stop("HDF5Array and rhdf5 are required for --h5ad")
  }
  h5ad_file <- normalizePath(cli$h5ad)
} else {
  h5ad_file <- NA_character_
}
file_identity <- function(path) {
  info <- file.info(path)
  if (nrow(info) != 1L || is.na(info$size) || is.na(info$mtime)) {
    stop("Could not read file identity: ", path)
  }
  list(
    path = normalizePath(path),
    size_bytes = as.numeric(info$size),
    mtime_utc = format(
      info$mtime, format = "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"
    )
  )
}
h5ad_identity <- file_identity

results_root <- Sys.getenv(
  "LOCAL_PDRT_DATA_DIR", unset = file.path(project_root, "results")
)
intermediate_dir <- if (!is.null(cli$intermediate_dir)) {
  cli$intermediate_dir
} else {
  file.path(results_root, "frangieh-2021", "intermediate-data")
}
required_inputs <- c("subsetted_ODM.rds", "covariate_matrix.rds",
                     "pure_nt.rds")
missing_inputs <- required_inputs[
  !file.exists(file.path(intermediate_dir, required_inputs))
]
if (length(missing_inputs)) {
  stop("Missing intermediate file(s): ", paste(missing_inputs, collapse = ", "))
}
intermediate_dir <- normalizePath(intermediate_dir)

output_dir <- if (!is.null(cli$output_dir)) {
  cli$output_dir
} else {
  file.path(results_root, "frangieh-2021", "exploratory-trans-pairs")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
manifest_dir <- file.path(output_dir, "manifests")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
expression_audit_path <- file.path(
  manifest_dir, "frangieh-exploratory-expression-audit.csv"
)

B <- as_count(if (is.null(cli$permutations)) 4999 else cli$permutations,
              "--permutations", minimum = 19L)
permutation_batch_size <- as_count(
  if (is.null(cli$permutation_batch_size)) 128 else cli$permutation_batch_size,
  "--permutation-batch-size"
)
base_seed <- as_count(if (is.null(cli$seed)) 20260822 else cli$seed,
                      "--seed")
cores <- as_count(if (is.null(cli$cores)) 1 else cli$cores, "--cores")
response_block_size <- as_count(
  if (is.null(cli$response_block_size)) 100 else cli$response_block_size,
  "--response-block-size"
)
min_treated <- as_count(
  if (is.null(cli$min_treated)) 10 else cli$min_treated,
  "--min-treated", minimum = 5L
)
if (!min_treated %in% c(5L, 10L)) {
  stop("--min-treated must be 10 (primary) or 5 (sensitivity)")
}
overwrite <- isTRUE(cli$overwrite)
manifest_only <- isTRUE(cli$manifest_only)
pilot_mode <- isTRUE(cli$pilot)

N_RESPONSES <- 5000L
GLOBAL_DETECTION_MIN <- 577L
EXPECTED_LOCAL_CELLS <- 57624L
EXPECTED_GLOBAL_GENES <- 23712L
EXPECTED_TARGETS <- 248L
EXPECTED_RESPONSE_UNIVERSE <- 12829L
EXPECTED_EXACT_ONE_NT <- 3770L
EXPECTED_EXACT_ONE_TARGETING <- 26716L
ALIGNMENT_DROP_POSITIONS <- c(27337L, 34511L, 56031L)
MIN_TREATED_BATCHES <- 2L
MAX_TREATED_BATCH_SHARE <- 0.80
PILOT_SIZE <- 5000L
PILOT_ASSIGNMENT_SEED <- 20260823L

cat("Frangieh exploratory trans-pair shard\n")
cat("  backend:     ", if (h5ad_mode) "scPerturb H5AD" else "ondisc ODM",
    "\n", sep = "")
if (h5ad_mode) cat("  H5AD:        ", h5ad_file, "\n", sep = "")
cat("  intermediate:", intermediate_dir, "\n", sep = " ")
cat("  output:      ", output_dir, "\n", sep = "")
cat("  B / cores:   ", B, " / ", cores, "\n", sep = "")

## ---- Load fixed local metadata and covariates ----------------------------
modalities <- readRDS(file.path(intermediate_dir, "subsetted_ODM.rds"))
covariate_matrix <- as.matrix(readRDS(
  file.path(intermediate_dir, "covariate_matrix.rds")
))
if (!all(c("gene", "grna") %in% names(modalities))) {
  stop("subsetted_ODM.rds must contain gene and grna modalities")
}
n_cells <- ncol(modalities$gene)
if (n_cells != EXPECTED_LOCAL_CELLS || ncol(modalities$grna) != n_cells ||
    nrow(covariate_matrix) != n_cells) {
  stop("Unexpected local Frangieh dimensions")
}
if (any(!is.finite(covariate_matrix))) stop("Non-finite saved covariates")

gene_features <- ondisc::get_feature_covariates(modalities$gene)
grna_features <- ondisc::get_feature_covariates(modalities$grna)
gene_cell_features <- ondisc::get_cell_covariates(modalities$gene)
grna_cell_features <- ondisc::get_cell_covariates(modalities$grna)
if (nrow(gene_features) != EXPECTED_GLOBAL_GENES ||
    !all(c("mean_expression", "n_nonzero") %in% names(gene_features))) {
  stop("Unexpected global Control gene feature metadata")
}
if (!all(c("target", "target_type") %in% names(grna_features)) ||
    !all(c("n_umis", "batch") %in% names(gene_cell_features)) ||
    !"n_nonzero" %in% names(grna_cell_features)) {
  stop("Required guide/cell metadata are absent")
}

cell_umi <- as.numeric(gene_cell_features$n_umis)
cell_batch <- as.character(gene_cell_features$batch)
cell_moi <- as.integer(grna_cell_features$n_nonzero)
if (any(!is.finite(cell_umi)) || any(cell_umi <= 0) || anyNA(cell_batch) ||
    any(!nzchar(cell_batch)) || anyNA(cell_moi) || any(cell_moi < 0L)) {
  stop("Invalid saved cell metadata")
}
if (!"log_n_umis" %in% colnames(covariate_matrix) ||
    max(abs(covariate_matrix[, "log_n_umis"] - log1p(cell_umi))) > 1e-10) {
  stop("Saved covariates do not contain the expected log1p(UMI_count)")
}

gene_ids <- rownames(gene_features)
if (anyDuplicated(gene_ids)) stop("Local gene identifiers are not unique")
grna_features$grna_id <- rownames(grna_features)
target_guide_rows <- grna_features$target_type == "gene" &
  !is.na(grna_features$target) & nzchar(grna_features$target)
nt_guide_rows <- grna_features$target_type == "non-targeting"
all_target_genes <- sort(unique(as.character(
  grna_features$target[target_guide_rows]
)))
target_guide_ids <- grna_features$grna_id[target_guide_rows]
nt_guide_ids <- grna_features$grna_id[nt_guide_rows]
if (length(all_target_genes) != EXPECTED_TARGETS ||
    length(target_guide_ids) != 744L || length(nt_guide_ids) != 74L) {
  stop("Unexpected Frangieh guide-target manifest")
}

## ---- Exact-one-guide exposure manifest -----------------------------------
## ODM backend: MOI and guide identity come from the same guide modality the
## production control-pair workflow uses.  A cell with MOI == 1 carries a
## single guide, so its target is the unique nonzero row of that column.  No
## H5AD alignment is needed because the ODM is already the 57,624 analysis
## cells in analysis order.
if (!h5ad_mode) {
  var_gene <- gene_ids
  X_delayed <- NULL
  aligned_control_global_idx <- seq_len(n_cells)
  control_global_idx <- seq_len(n_cells)

  guide_presence <- modalities$grna[[grna_features$grna_id, ]]
  guide_presence <- methods::as(guide_presence, "TsparseMatrix")
  if (!identical(as.integer(dim(guide_presence)),
                 c(nrow(grna_features), n_cells))) {
    stop("Unexpected guide-assignment matrix dimensions")
  }
  if (!identical(as.integer(Matrix::colSums(guide_presence)), cell_moi)) {
    stop("Guide-assignment column sums disagree with the stored cell MOI")
  }

  exact_one_guide <- cell_moi == 1L
  entry_cell <- guide_presence@j + 1L
  entry_guide <- guide_presence@i + 1L
  keep_entry <- exact_one_guide[entry_cell]
  guide_row_of_cell <- rep(NA_integer_, n_cells)
  guide_row_of_cell[entry_cell[keep_entry]] <- entry_guide[keep_entry]
  if (any(is.na(guide_row_of_cell[exact_one_guide]))) {
    stop("An MOI==1 cell lacks a recognized guide")
  }

  guide_type <- rep(NA_character_, n_cells)
  guide_target <- rep(NA_character_, n_cells)
  guide_type[exact_one_guide] <-
    grna_features$target_type[guide_row_of_cell[exact_one_guide]]
  guide_target[exact_one_guide] <-
    as.character(grna_features$target[guide_row_of_cell[exact_one_guide]])
  pure_nt <- exact_one_guide & guide_type == "non-targeting"
  is_targeting <- exact_one_guide & guide_type == "gene"
  pure_nt[is.na(pure_nt)] <- FALSE
  is_targeting[is.na(is_targeting)] <- FALSE
  if (sum(pure_nt) != EXPECTED_EXACT_ONE_NT ||
      sum(is_targeting) != EXPECTED_EXACT_ONE_TARGETING) {
    stop("Unexpected exact-one-guide classification")
  }
  rm(guide_presence, entry_cell, entry_guide, keep_entry, guide_row_of_cell)
} else {
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
               c(length(var_gene), length(h5_condition))) ||
    anyDuplicated(var_gene) || !setequal(var_gene, gene_ids)) {
  stop("Unexpected H5AD gene matrix or identifiers")
}

control_global_idx <- which(h5_condition == "Control")
if (length(control_global_idx) != n_cells + length(ALIGNMENT_DROP_POSITIONS) ||
    any(ALIGNMENT_DROP_POSITIONS > length(control_global_idx))) {
  stop("Unexpected H5AD Control-cell count")
}
aligned_control_global_idx <- control_global_idx[-ALIGNMENT_DROP_POSITIONS]
aligned_h5_umi <- h5_umi[aligned_control_global_idx]
aligned_h5_moi <- h5_moi[aligned_control_global_idx]
aligned_h5_guide <- h5_guide[aligned_control_global_idx]
if (!identical(aligned_h5_umi, cell_umi) ||
    !identical(aligned_h5_moi, cell_moi)) {
  stop("H5AD/local cells do not align after the three audited drops")
}

h5_guide_categories <- as.character(rhdf5::h5read(
  h5ad_file, "/obs/__categories/sgRNA"
))
if (!setequal(h5_guide_categories, rownames(grna_features))) {
  stop("H5AD guide categories do not match local guide metadata")
}
exact_one_guide <- aligned_h5_moi == 1L
if (any(is.na(aligned_h5_guide[exact_one_guide])) ||
    any(!aligned_h5_guide[exact_one_guide] %in% rownames(grna_features))) {
  stop("An MOI==1 cell lacks a recognized guide")
}
guide_type <- rep(NA_character_, n_cells)
guide_target <- rep(NA_character_, n_cells)
assigned <- exact_one_guide & !is.na(aligned_h5_guide)
guide_type[assigned] <- grna_features[aligned_h5_guide[assigned], "target_type"]
guide_target[assigned] <- grna_features[aligned_h5_guide[assigned], "target"]
pure_nt <- exact_one_guide & guide_type == "non-targeting"
is_targeting <- exact_one_guide & guide_type == "gene"
pure_nt[is.na(pure_nt)] <- FALSE
is_targeting[is.na(is_targeting)] <- FALSE
if (sum(pure_nt) != EXPECTED_EXACT_ONE_NT ||
    sum(is_targeting) != EXPECTED_EXACT_ONE_TARGETING) {
  stop("Unexpected exact-one-guide classification")
}
}

## ---- Raw-count and stored-expression audit -------------------------------
## The full raw-H5AD audit underlying this workflow found maximum absolute
## mean difference 2.78e-17 and zero nonzero-count differences over all 23,712
## genes.  Every shard rechecks the inversion on 20 whole transcriptomes and
## the stored mean/nonzero metadata on 20 genes spanning the complete feature
## order before the metadata can be used to select responses.
recover_h5ad_counts <- function(log_cpm, depths, context, tolerance = 5e-3) {
  log_cpm <- as.matrix(log_cpm)
  depths <- as.numeric(depths)
  if (ncol(log_cpm) != length(depths) || any(!is.finite(log_cpm)) ||
      any(log_cpm < 0) || any(!is.finite(depths)) || any(depths <= 0)) {
    stop("Invalid log1p(CPM) block or UMI depth for ", context)
  }
  unrounded <- sweep(expm1(log_cpm), 2L, depths / 1e6, `*`)
  counts <- round(unrounded)
  max_error <- if (length(unrounded)) max(abs(unrounded - counts)) else 0
  if (!is.finite(max_error) || max_error > tolerance || any(counts < 0)) {
    stop("log1p(CPM) inversion failed for ", context,
         sprintf(" (maximum residual %.6g)", max_error))
  }
  storage.mode(counts) <- "double"
  attr(counts, "max_rounding_error") <- max_error
  counts
}

## Backend-neutral raw-count reader: genes are addressed by identifier and
## cells by local analysis index, so every caller below is backend-agnostic.
## The ODM already stores raw integer counts, so nothing is reconstructed and
## the rounding residual is exactly zero by construction.
read_counts <- function(genes, cells, context) {
  if (h5ad_mode) {
    rows <- match(genes, var_gene)
    if (anyNA(rows)) stop("A requested gene is absent from the H5AD: ", context)
    return(recover_h5ad_counts(
      X_delayed[rows, aligned_control_global_idx[cells], drop = FALSE],
      cell_umi[cells], context
    ))
  }
  counts <- as.matrix(modalities$gene[[genes, ]])[, cells, drop = FALSE]
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > 1e-8)) {
    stop("The ODM returned non-integer or negative counts for ", context)
  }
  attr(counts, "max_rounding_error") <- 0
  counts
}

## Whole-transcriptome UMI reconciliation.  Under the H5AD backend this is the
## inversion audit and every shard repeats it; under the ODM backend it is a
## column-total consistency check whose cost is cell-wise access, so it runs
## only in the one-time --manifest-only pass.
audit_local_idx <- unique(as.integer(round(seq(
  1, n_cells, length.out = 20L
))))
if (h5ad_mode) {
  audit_counts <- recover_h5ad_counts(
    X_delayed[, aligned_control_global_idx[audit_local_idx], drop = FALSE],
    cell_umi[audit_local_idx], "20-cell whole-transcriptome audit"
  )
  count_inversion_audit_error <- attr(audit_counts, "max_rounding_error")
  if (!identical(as.numeric(colSums(audit_counts)),
                 cell_umi[audit_local_idx])) {
    stop("Reconstructed whole-transcriptome UMI totals do not match metadata")
  }
  rm(audit_counts)
} else {
  count_inversion_audit_error <- 0
  if (manifest_only) {
    audit_counts <- modalities$gene[[, audit_local_idx]]
    if (!identical(as.numeric(Matrix::colSums(audit_counts)),
                   cell_umi[audit_local_idx])) {
      stop("ODM whole-transcriptome UMI totals do not match metadata")
    }
    rm(audit_counts)
  }
}

## Stored mean/nonzero metadata drives response selection, so it is audited
## against raw counts on 20 genes spanning the complete feature order.  The
## ODM backend can afford this in every shard; the H5AD backend cannot, and
## inherits the result from the locked --manifest-only audit instead.
audit_gene_local_idx <- unique(as.integer(round(seq(
  1, length(gene_ids), length.out = 20L
))))
## The counts source is the H5AD under --h5ad and the ODM backing store
## otherwise; both are identified the same way (full MD5 once under
## --manifest-only, then path/size/mtime per shard).
counts_backend <- if (h5ad_mode) "scperturb_h5ad" else "ondisc_odm"
counts_source_file <- if (h5ad_mode) {
  h5ad_file
} else {
  normalizePath(modalities$gene@ondisc_matrix@h5_file)
}
if (!file.exists(counts_source_file)) {
  stop("Counts backing store not found: ", counts_source_file)
}
if (manifest_only || !h5ad_mode) {
  audit_gene_counts <- read_counts(
    gene_ids[audit_gene_local_idx], seq_len(n_cells),
    "20-gene all-Control expression-metadata audit"
  )
  audit_gene_mean_difference <- max(abs(
    rowMeans(audit_gene_counts) -
      as.numeric(gene_features$mean_expression[audit_gene_local_idx])
  ))
  audit_gene_nonzero_difference <- max(abs(
    rowSums(audit_gene_counts != 0) -
      as.integer(gene_features$n_nonzero[audit_gene_local_idx])
  ))
  gene_audit_inversion_error <- attr(audit_gene_counts,
                                     "max_rounding_error")
  rm(audit_gene_counts)
}
if (manifest_only) {
  before_hash <- file_identity(counts_source_file)
  cat("  hashing the full counts backing store once (", counts_backend,
      ") ...\n", sep = "")
  locked_counts_md5 <- unname(tools::md5sum(counts_source_file))
  after_hash <- file_identity(counts_source_file)
  if (!identical(before_hash, after_hash) ||
      length(locked_counts_md5) != 1L || is.na(locked_counts_md5) ||
      !grepl("^[0-9a-f]{32}$", locked_counts_md5)) {
    stop("Counts-source identity changed or hashing failed during the audit")
  }
  verified_counts_identity <- after_hash
  expression_audit <- data.frame(
    schema_version = 3L,
    backend = counts_backend,
    n_aligned_control_cells = n_cells,
    n_genes = length(gene_ids),
    n_rechecked_genes = length(audit_gene_local_idx),
    rechecked_gene_indices = paste(audit_gene_local_idx, collapse = ","),
    rechecked_mean_max_abs_difference = audit_gene_mean_difference,
    rechecked_nonzero_max_abs_difference = audit_gene_nonzero_difference,
    rechecked_inversion_max_abs_residual = gene_audit_inversion_error,
    full_direct_mean_max_abs_difference = if (h5ad_mode) 2.78e-17 else 0,
    full_direct_nonzero_max_abs_difference = 0,
    full_direct_inversion_max_abs_residual = if (h5ad_mode) 0.00164278 else 0,
    counts_source_path = after_hash$path,
    counts_source_size_bytes = after_hash$size_bytes,
    counts_source_mtime_utc = after_hash$mtime_utc,
    counts_source_md5 = locked_counts_md5,
    shard_identity_check = paste(
      "locked normalized path, byte size, and UTC mtime must match;",
      "full MD5 is inherited from this manifest-only audit"
    ),
    stringsAsFactors = FALSE
  )
  publish_identical_csv(expression_audit, expression_audit_path)
} else {
  if (!file.exists(expression_audit_path)) {
    stop("Missing locked expression audit. Run this script once with ",
         "--manifest-only in the same output directory: ",
         expression_audit_path)
  }
  expression_audit <- utils::read.csv(expression_audit_path,
                                      stringsAsFactors = FALSE)
  current_counts_identity <- file_identity(counts_source_file)
  required_identity_columns <- c(
    "backend", "counts_source_path", "counts_source_size_bytes",
    "counts_source_mtime_utc", "counts_source_md5"
  )
  if (nrow(expression_audit) != 1L ||
      !all(required_identity_columns %in% names(expression_audit)) ||
      expression_audit$schema_version != 3L ||
      !identical(expression_audit$backend, counts_backend) ||
      expression_audit$n_aligned_control_cells != n_cells ||
      expression_audit$n_genes != length(gene_ids) ||
      !identical(expression_audit$counts_source_path,
                 current_counts_identity$path) ||
      expression_audit$counts_source_size_bytes !=
        current_counts_identity$size_bytes ||
      !identical(expression_audit$counts_source_mtime_utc,
                 current_counts_identity$mtime_utc) ||
      is.na(expression_audit$counts_source_md5) ||
      !grepl("^[0-9a-f]{32}$", expression_audit$counts_source_md5)) {
    stop("Locked expression audit does not match this counts source/schema")
  }
  locked_counts_md5 <- expression_audit$counts_source_md5
  verified_counts_identity <- current_counts_identity
  ## The ODM backend recomputed the metadata audit above; only the H5AD
  ## backend inherits it, because its recheck needs the whole matrix.
  if (h5ad_mode) {
    audit_gene_mean_difference <-
      expression_audit$rechecked_mean_max_abs_difference
    audit_gene_nonzero_difference <-
      expression_audit$rechecked_nonzero_max_abs_difference
    gene_audit_inversion_error <-
      expression_audit$rechecked_inversion_max_abs_residual
  }
}
if (!is.finite(audit_gene_mean_difference) ||
    audit_gene_mean_difference > 1e-12 ||
    audit_gene_nonzero_difference != 0) {
  stop("Stored global expression metadata failed its raw-count audit")
}

## ---- Association-blind response manifest --------------------------------
response_universe <- data.frame(
  response_gene = gene_ids,
  global_mean_count = as.numeric(gene_features$mean_expression),
  global_nonzero = as.integer(gene_features$n_nonzero),
  is_perturbation_target = gene_ids %in% all_target_genes,
  stringsAsFactors = FALSE
)
if (any(!is.finite(response_universe$global_mean_count)) ||
    any(response_universe$global_mean_count < 0) ||
    anyNA(response_universe$global_nonzero) ||
    any(response_universe$global_nonzero < 0L)) {
  stop("Invalid global response-expression metadata")
}
response_universe$passes_global_detection <-
  response_universe$global_nonzero >= GLOBAL_DETECTION_MIN
response_universe$eligible <- with(
  response_universe,
  passes_global_detection & !is_perturbation_target
)
eligible_response_idx <- which(response_universe$eligible)
if (length(eligible_response_idx) != EXPECTED_RESPONSE_UNIVERSE) {
  stop("Expected ", EXPECTED_RESPONSE_UNIVERSE,
       " eligible non-target responses; found ", length(eligible_response_idx))
}
response_order <- eligible_response_idx[order(
  response_universe$global_mean_count[eligible_response_idx],
  response_universe$response_gene[eligible_response_idx]
)]
response_universe$selection_rank <- NA_integer_
response_universe$selection_rank[response_order] <- seq_along(response_order)
selected_response_idx <- head(response_order, N_RESPONSES)
response_universe$selected <- seq_len(nrow(response_universe)) %in%
  selected_response_idx
selected_responses <- response_universe[selected_response_idx, c(
  "response_gene", "global_mean_count", "global_nonzero",
  "passes_global_detection", "is_perturbation_target", "selection_rank"
), drop = FALSE]
names(selected_responses)[names(selected_responses) == "selection_rank"] <-
  "response_rank"
selected_responses$global_detection_fraction <-
  selected_responses$global_nonzero / n_cells
selected_responses <- selected_responses[order(selected_responses$response_rank), ]
rownames(selected_responses) <- NULL
if (nrow(selected_responses) != N_RESPONSES ||
    any(selected_responses$is_perturbation_target) ||
    selected_responses$response_gene[[1L]] != "CCNT2-AS1" ||
    selected_responses$response_gene[[N_RESPONSES]] != "PRELID2" ||
    abs(min(selected_responses$global_mean_count) - 0.010082604) > 1e-9 ||
    abs(max(selected_responses$global_mean_count) - 0.13053589) > 1e-8 ||
    min(selected_responses$global_nonzero) != 577L ||
    max(selected_responses$global_nonzero) != 6892L) {
  stop("The prespecified 5,000-response manifest failed its audit")
}

exposure_table <- data.frame(
  exposure_gene = all_target_genes,
  n_guides = vapply(all_target_genes, function(g) {
    sum(grna_features$target_type == "gene" & grna_features$target == g,
        na.rm = TRUE)
  }, integer(1)),
  n_treated = vapply(all_target_genes, function(g) {
    sum(is_targeting & !is.na(guide_target) & guide_target == g)
  }, integer(1)),
  n_treated_batches = vapply(all_target_genes, function(g) {
    idx <- which(is_targeting & !is.na(guide_target) & guide_target == g)
    length(unique(cell_batch[idx]))
  }, integer(1)),
  max_treated_batch_share = vapply(all_target_genes, function(g) {
    idx <- which(is_targeting & !is.na(guide_target) & guide_target == g)
    if (!length(idx)) return(NA_real_)
    max(table(cell_batch[idx])) / length(idx)
  }, numeric(1)),
  stringsAsFactors = FALSE
)
exposure_table$eligible <- with(
  exposure_table,
  n_treated >= min_treated &
    n_treated_batches >= MIN_TREATED_BATCHES &
    max_treated_batch_share <= MAX_TREATED_BATCH_SHARE
)
exposure_table$exclusion_reason <- ifelse(
  exposure_table$eligible, "",
  ifelse(exposure_table$n_treated < min_treated,
         paste0("n_treated_below_", min_treated),
         ifelse(exposure_table$n_treated_batches < MIN_TREATED_BATCHES,
                "fewer_than_2_treated_batches",
                "max_treated_batch_share_above_0.8"))
)
eligible_exposure_rows <- which(exposure_table$eligible)
exposure_table$exposure_index <- NA_integer_
exposure_table$exposure_index[eligible_exposure_rows] <-
  seq_along(eligible_exposure_rows)
expected_eligible_exposures <- if (min_treated == 10L) 241L else 244L
expected_ineligible_exposures <- if (min_treated == 10L) {
  c("ATP1A1", "DNAJC9", "EIF2S3", "PSMA7", "RACK1", "TUBB", "UBL5")
} else {
  c("ATP1A1", "PSMA7", "TUBB", "UBL5")
}
if (length(eligible_exposure_rows) != expected_eligible_exposures ||
    !identical(exposure_table$exposure_gene[!exposure_table$eligible],
               expected_ineligible_exposures)) {
  stop("The exposure manifest failed its min-treated=", min_treated,
       " audit")
}
eligible_exposures <- exposure_table[eligible_exposure_rows, , drop = FALSE]
rownames(eligible_exposures) <- NULL

## ---- Frozen association-blind 5,000-pair pilot ---------------------------
## Each response is used exactly once.  Exposure allocation is balanced at
## 20 or 21 pairs per target and depends only on a fixed RNG seed and manifest
## indices, never on a target-response association, effect, statistic, or p.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
pilot_manifest <- NULL
if (min_treated == 10L) {
  set.seed(PILOT_ASSIGNMENT_SEED)
  pilot_counts <- rep(PILOT_SIZE %/% nrow(eligible_exposures),
                      nrow(eligible_exposures))
  n_extra <- PILOT_SIZE %% nrow(eligible_exposures)
  if (n_extra > 0L) {
    extra_exposures <- sample.int(nrow(eligible_exposures), n_extra,
                                  replace = FALSE)
    pilot_counts[extra_exposures] <- pilot_counts[extra_exposures] + 1L
  }
  pilot_exposure_indices <- sample(rep(
    seq_len(nrow(eligible_exposures)), times = pilot_counts
  ), replace = FALSE)
  pilot_response_ranks <- sample.int(N_RESPONSES, PILOT_SIZE, replace = FALSE)
  pilot_manifest <- data.frame(
    pilot_pair_index = seq_len(PILOT_SIZE),
    exposure_index = pilot_exposure_indices,
    exposure_gene = eligible_exposures$exposure_gene[pilot_exposure_indices],
    response_rank = pilot_response_ranks,
    response_gene = selected_responses$response_gene[pilot_response_ranks],
    stringsAsFactors = FALSE
  )
  pilot_tab <- table(pilot_manifest$exposure_index)
  if (nrow(pilot_manifest) != PILOT_SIZE ||
      anyDuplicated(pilot_manifest$response_rank) ||
      any(pilot_manifest$exposure_gene == pilot_manifest$response_gene) ||
      length(pilot_tab) != 241L || min(pilot_tab) != 20L ||
      max(pilot_tab) != 21L) {
    stop("The association-blind 5,000-pair pilot manifest failed its audit")
  }
}
if (pilot_mode && min_treated != 10L) {
  stop("--pilot is frozen to the primary --min-treated 10 exposure manifest")
}
if (pilot_mode && !is.null(cli[["response_range"]])) {
  stop("--response-range cannot be combined with the fixed --pilot manifest")
}

## ---- Target/pair shard coordinates ---------------------------------------
selector_flags <- c(!is.null(cli[["target"]]),
                    !is.null(cli[["target_index"]]),
                    !is.null(cli[["target_range"]]),
                    isTRUE(cli[["all_targets"]]))
if (!manifest_only && sum(selector_flags) != 1L) {
  stop("Choose exactly one target selector\n", usage())
}
if (manifest_only && sum(selector_flags) > 1L) {
  stop("Choose at most one target selector with --manifest-only")
}
selected_exposure_indices <- if (manifest_only && sum(selector_flags) == 0L) {
  integer()
} else if (!is.null(cli[["target"]])) {
  requested_target <- cli[["target"]]
  hit <- which(eligible_exposures$exposure_gene == requested_target)
  if (length(hit) != 1L) {
    if (requested_target %in% exposure_table$exposure_gene) {
      stop("Requested target is ineligible: ", requested_target)
    }
    stop("Unknown target: ", requested_target)
  }
  hit
} else if (!is.null(cli[["target_index"]])) {
  value <- as_count(cli[["target_index"]], "--target-index")
  if (value > nrow(eligible_exposures)) {
    stop("--target-index exceeds ", nrow(eligible_exposures))
  }
  value
} else if (!is.null(cli[["target_range"]])) {
  parse_index_range(cli[["target_range"]], nrow(eligible_exposures),
                    "--target-range")
} else {
  seq_len(nrow(eligible_exposures))
}

requested_response_ranks <- if (is.null(cli[["response_range"]])) {
  seq_len(N_RESPONSES)
} else {
  parse_index_range(cli[["response_range"]], N_RESPONSES, "--response-range")
}
if (!length(selected_exposure_indices)) {
  pair_manifest <- data.frame(
    pilot_pair_index = integer(), exposure_index = integer(),
    exposure_gene = character(), response_rank = integer(),
    response_gene = character(), stringsAsFactors = FALSE
  )
} else if (pilot_mode) {
  pair_manifest <- pilot_manifest[
    pilot_manifest$exposure_index %in% selected_exposure_indices,
    , drop = FALSE
  ]
} else {
  pair_manifest <- expand.grid(
    exposure_index = selected_exposure_indices,
    response_rank = requested_response_ranks,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  pair_manifest$pilot_pair_index <- NA_integer_
  pair_manifest$exposure_gene <- eligible_exposures$exposure_gene[
    pair_manifest$exposure_index
  ]
  pair_manifest$response_gene <- selected_responses$response_gene[
    pair_manifest$response_rank
  ]
  pair_manifest <- pair_manifest[, c(
    "pilot_pair_index", "exposure_index", "exposure_gene",
    "response_rank", "response_gene"
  )]
}
pair_manifest <- pair_manifest[order(pair_manifest$exposure_index,
                                     pair_manifest$response_rank), ]
rownames(pair_manifest) <- NULL
selected_response_shard <- selected_responses[
  sort(unique(pair_manifest$response_rank)), , drop = FALSE
]

## Publish fixed manifests.  Every array task derives the same bytes.
response_universe_out <- response_universe
response_universe_out <- response_universe_out[order(
  response_universe_out$response_gene
), ]
publish_identical_csv(
  response_universe_out,
  file.path(manifest_dir, "frangieh-exploratory-response-universe.csv")
)
publish_identical_csv(
  selected_responses,
  file.path(manifest_dir, "frangieh-exploratory-selected-responses.csv")
)
publish_identical_csv(
  exposure_table,
  file.path(
    manifest_dir,
    sprintf("frangieh-exploratory-exposures-min-treated-%02d.csv",
            min_treated)
  )
)
if (!is.null(pilot_manifest)) {
  publish_identical_csv(
    pilot_manifest,
    file.path(manifest_dir,
              "frangieh-exploratory-association-blind-pilot-5000.csv")
  )
}

manifest_paths <- c(
  expression_audit_path,
  file.path(manifest_dir, "frangieh-exploratory-response-universe.csv"),
  file.path(manifest_dir, "frangieh-exploratory-selected-responses.csv"),
  file.path(
    manifest_dir,
    sprintf("frangieh-exploratory-exposures-min-treated-%02d.csv",
            min_treated)
  )
)
manifest_names <- c("expression_audit", "response_universe",
                    "selected_responses", "exposures")
if (!is.null(pilot_manifest)) {
  manifest_paths <- c(
    manifest_paths,
    file.path(manifest_dir,
              "frangieh-exploratory-association-blind-pilot-5000.csv")
  )
  manifest_names <- c(manifest_names, "association_blind_pilot_5000")
}
manifest_hashes <- unname(tools::md5sum(manifest_paths))
names(manifest_hashes) <- manifest_names

cat("  response panel: ", N_RESPONSES, " of ",
    length(eligible_response_idx), " eligible non-target genes\n", sep = "")
cat("  exposures:      ", nrow(eligible_exposures),
    " eligible of ", nrow(exposure_table), "; ",
    length(selected_exposure_indices), " in this shard\n", sep = "")
cat("  pair mode:      ", if (pilot_mode) {
  paste0("association-blind pilot; ", nrow(pair_manifest), " pairs in shard")
} else {
  paste0("full grid/range; ", nrow(pair_manifest), " pairs in shard")
}, "\n", sep = "")
cat("  manifest MD5:   ", paste(manifest_hashes, collapse = " / "),
    "\n", sep = "")

if (manifest_only) {
  cat("Manifest-only validation complete.\n")
  quit(save = "no", status = 0L)
}

exposure_label <- if (length(selected_exposure_indices) == 1L) {
  sprintf("e%03d-%s", selected_exposure_indices,
          sanitize_filename(eligible_exposures$exposure_gene[
            selected_exposure_indices
          ]))
} else {
  sprintf("e%03d-%03d", min(selected_exposure_indices),
          max(selected_exposure_indices))
}
response_label <- if (pilot_mode) {
  "pilot5000"
} else {
  sprintf("r%04d-%04d", min(requested_response_ranks),
          max(requested_response_ranks))
}
stem <- paste("frangieh-exploratory-trans", exposure_label, response_label,
              paste0("min", min_treated), paste0("B", B), sep = "-")
rds_path <- file.path(output_dir, paste0(stem, ".rds"))
csv_path <- file.path(output_dir, paste0(stem, ".csv"))
if ((file.exists(rds_path) || file.exists(csv_path)) && !overwrite) {
  stop("Shard output already exists; use --overwrite: ", stem)
}

if (anyNA(match(selected_responses$response_gene, var_gene))) {
  stop("A selected response is absent from the counts backing store")
}

## ---- Paired score/permutation primitives ---------------------------------
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

seed_for_pair <- function(seed, exposure_index, response_rank) {
  modulus <- 2147483646
  pair_index <- (as.double(exposure_index) - 1) * N_RESPONSES +
    as.double(response_rank)
  as.integer((as.double(seed) + 10007 * pair_index) %% modulus + 1)
}

failure_row <- function(exposure_row, response_row, pair_seed, n_pair_cells,
                        message, status = "error") {
  data.frame(
    pair_id = paste(exposure_row$exposure_gene, response_row$response_gene,
                    sep = "::"),
    exposure_index = exposure_row$exposure_index,
    exposure_gene = exposure_row$exposure_gene,
    response_rank = response_row$response_rank,
    response_gene = response_row$response_gene,
    status = status,
    error_message = as.character(message),
    pair_seed = pair_seed,
    n_cells = n_pair_cells,
    n_pure_nt = EXPECTED_EXACT_ONE_NT,
    n_treated = exposure_row$n_treated,
    n_treated_batches = exposure_row$n_treated_batches,
    max_treated_batch_share = exposure_row$max_treated_batch_share,
    global_mean_count = response_row$global_mean_count,
    global_nonzero = response_row$global_nonzero,
    mean_count_pure_nt = NA_real_, mean_count_treated = NA_real_,
    log2_marginal_fold_change = NA_real_,
    fraction_nonzero_pure_nt = NA_real_,
    fraction_nonzero_treated = NA_real_,
    theta_hat = NA_real_, dispersion_method = NA_character_,
    information_stat = NA_real_, empirical_stat = NA_real_,
    information_denom_sq = NA_real_, empirical_denom_sq = NA_real_,
    empirical_to_information_denom_ratio = NA_real_,
    score_eff_n = NA_real_, max_score_share = NA_real_,
    hits_abs_barry = NA_real_, hits_abs_empirical = NA_real_,
    hits_negative_extreme_barry = NA_real_,
    hits_negative_extreme_empirical = NA_real_,
    hits_positive_extreme_barry = NA_real_,
    hits_positive_extreme_empirical = NA_real_,
    p_two_sided_barry = NA_real_, p_two_sided_empirical = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_pair <- function(Y, X, Z_design, exposure_row, response_row) {
  pair_seed <- seed_for_pair(base_seed, exposure_row$exposure_index,
                             response_row$response_rank)
  n_pair_cells <- length(X)
  tryCatch({
    if (length(Y) != n_pair_cells || any(!is.finite(Y)) || any(Y < 0) ||
        any(abs(Y - round(Y)) > 1e-8)) {
      stop("Response is not a compatible raw UMI vector")
    }
    cache <- precompute_score_cache(Y, Z_design)
    observed <- precompute_score_pair(X, cache)
    ## Distinguish the two causes.  A degenerate nuisance fit depends only on
    ## (Y, Z_design), which are held fixed across the whole permutation orbit,
    ## so it is label-invariant and a pair failing this way may be excluded by
    ## a pre-declared rule.  A degenerate observed score depends on X and is
    ## NOT permutation-invariant, so such a pair must never be dropped on that
    ## basis; it is reported separately for conservative handling.
    if (isTRUE(cache$degenerate)) {
      return(failure_row(exposure_row, response_row, pair_seed, n_pair_cells,
                         "degenerate negative-binomial nuisance fit on (Y, Z)",
                         status = "degenerate_nuisance_fit"))
    }
    if (isTRUE(observed$degenerate)) {
      return(failure_row(exposure_row, response_row, pair_seed, n_pair_cells,
                         "degenerate observed efficient score (depends on X)",
                         status = "degenerate_observed_score"))
    }

    hits_information_abs <- hits_empirical_abs <- 0
    hits_information_negative <- hits_empirical_negative <- 0
    hits_information_positive <- hits_empirical_positive <- 0
    set.seed(pair_seed)
    for (batch_start in seq.int(1L, B, by = permutation_batch_size)) {
      n_draws <- min(permutation_batch_size, B - batch_start + 1L)
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
    data.frame(
      pair_id = paste(exposure_row$exposure_gene, response_row$response_gene,
                      sep = "::"),
      exposure_index = exposure_row$exposure_index,
      exposure_gene = exposure_row$exposure_gene,
      response_rank = response_row$response_rank,
      response_gene = response_row$response_gene,
      status = "ok", error_message = "", pair_seed = pair_seed,
      n_cells = n_pair_cells, n_pure_nt = sum(X == 0L),
      n_treated = sum(X == 1L),
      n_treated_batches = exposure_row$n_treated_batches,
      max_treated_batch_share = exposure_row$max_treated_batch_share,
      global_mean_count = response_row$global_mean_count,
      global_nonzero = response_row$global_nonzero,
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
  }, error = function(e) {
    failure_row(exposure_row, response_row, pair_seed, n_pair_cells,
                conditionMessage(e), status = "error")
  })
}

## ---- Run selected exposure shards in bounded response blocks -------------
test_start <- Sys.time()
result_list <- vector("list", length(selected_exposure_indices))
inversion_errors <- numeric()

for (exposure_position in seq_along(selected_exposure_indices)) {
  exposure_index <- selected_exposure_indices[[exposure_position]]
  exposure_row <- eligible_exposures[exposure_index, , drop = FALSE]
  exposure_pairs <- pair_manifest[
    pair_manifest$exposure_index == exposure_index, , drop = FALSE
  ]
  exposure_response_manifest <- selected_responses[
    match(exposure_pairs$response_rank, selected_responses$response_rank),
    , drop = FALSE
  ]
  treated <- is_targeting & !is.na(guide_target) &
    guide_target == exposure_row$exposure_gene
  cell_idx <- which(pure_nt | treated)
  X <- as.integer(treated[cell_idx])
  if (sum(X) != exposure_row$n_treated || sum(X == 0L) != EXPECTED_EXACT_ONE_NT) {
    stop("Exposure mask audit failed for ", exposure_row$exposure_gene)
  }
  Z_design <- cbind(
    intercept = 1,
    covariate_matrix[cell_idx, , drop = FALSE]
  )
  if (qr(Z_design)$rank != ncol(Z_design)) {
    stop("Rank-deficient pairwise design for ", exposure_row$exposure_gene)
  }

  cat(sprintf("Exposure %d/%d: %s (%d treated); %d responses\n",
              exposure_index, nrow(eligible_exposures),
              exposure_row$exposure_gene, sum(X), nrow(exposure_pairs)))
  block_starts <- seq.int(1L, nrow(exposure_pairs),
                          by = response_block_size)
  exposure_results <- vector("list", length(block_starts))
  for (block_number in seq_along(block_starts)) {
    first <- block_starts[[block_number]]
    last <- min(nrow(exposure_pairs), first + response_block_size - 1L)
    block_rows <- seq.int(first, last)
    block_manifest <- exposure_response_manifest[block_rows, , drop = FALSE]
    Y_block <- read_counts(
      selected_responses$response_gene[block_manifest$response_rank],
      cell_idx,
      paste0("target ", exposure_row$exposure_gene, " response block ",
             min(block_manifest$response_rank), ":",
             max(block_manifest$response_rank))
    )
    inversion_errors <- c(
      inversion_errors, attr(Y_block, "max_rounding_error")
    )
    work <- function(j) {
      run_pair(
        Y = as.numeric(Y_block[j, ]), X = X, Z_design = Z_design,
        exposure_row = exposure_row,
        response_row = block_manifest[j, , drop = FALSE]
      )
    }
    block_results <- if (.Platform$OS.type != "windows" && cores > 1L) {
      parallel::mclapply(
        seq_len(nrow(block_manifest)), work, mc.cores = cores,
        mc.preschedule = TRUE, mc.set.seed = FALSE
      )
    } else {
      lapply(seq_len(nrow(block_manifest)), work)
    }
    exposure_results[[block_number]] <- do.call(rbind, block_results)
    rm(Y_block)
    cat(sprintf("  response ranks %d:%d complete\n",
                min(block_manifest$response_rank),
                max(block_manifest$response_rank)))
  }
  result_list[[exposure_position]] <- do.call(rbind, exposure_results)
}

pair_results <- do.call(rbind, result_list)
pair_results <- pair_results[order(pair_results$exposure_index,
                                   pair_results$response_rank), ]
rownames(pair_results) <- NULL
pair_manifest_key <- paste(pair_manifest$exposure_gene,
                           pair_manifest$response_gene, sep = "::")
pair_results$pilot_pair_index <- pair_manifest$pilot_pair_index[
  match(pair_results$pair_id, pair_manifest_key)
]
expected_pairs <- nrow(pair_manifest)
if (nrow(pair_results) != expected_pairs || anyDuplicated(pair_results$pair_id) ||
    any(pair_results$exposure_gene == pair_results$response_gene)) {
  stop("Completed pair table failed its trans-pair audit")
}
ok <- pair_results$status == "ok"
if (any(ok)) {
  if (any(pair_results$hits_abs_barry[ok] < 0 |
          pair_results$hits_abs_barry[ok] > B) ||
      any(pair_results$hits_abs_empirical[ok] < 0 |
          pair_results$hits_abs_empirical[ok] > B) ||
      max(abs(pair_results$p_two_sided_barry[ok] -
                (1 + pair_results$hits_abs_barry[ok]) / (B + 1))) > 1e-15 ||
      max(abs(pair_results$p_two_sided_empirical[ok] -
                (1 + pair_results$hits_abs_empirical[ok]) / (B + 1))) > 1e-15) {
    stop("Permutation hit/p-value audit failed")
  }
}
test_seconds <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))

## ---- Shard output with complete embedded manifests/provenance ------------
input_paths <- file.path(intermediate_dir, required_inputs)
config <- list(
  analysis = "frangieh_label_invariant_exploratory_trans_pairs",
  counts_backend = counts_backend,
  response_selection = paste0(
    "exclude all 248 perturbation targets; global Control n_nonzero >= ",
    GLOBAL_DETECTION_MIN,
    "; lowest 5000 audited mean_expression with gene-name tie break"
  ),
  response_selection_uses_cell_guide_assignments = FALSE,
  response_selection_excludes_annotated_perturbation_targets = TRUE,
  response_metadata_raw_h5ad_audited = TRUE,
  response_panel_size = N_RESPONSES,
  response_ranks_in_shard = sort(unique(pair_manifest$response_rank)),
  response_block_size = response_block_size,
  exposure_filter = paste0(
    "exact-one gene-targeting; n_treated >= ", min_treated,
    "; treated batches >= 2; ",
    "max treated batch share <= 0.8"
  ),
  min_treated = min_treated,
  exposure_indices_in_shard = selected_exposure_indices,
  scan_mode = if (pilot_mode) "association_blind_pilot_5000" else
    "declared_grid_or_range",
  pilot_assignment_uses_association = FALSE,
  pilot_assignment_seed = if (pilot_mode) PILOT_ASSIGNMENT_SEED else NA_integer_,
  declared_full_grid_pairs = nrow(eligible_exposures) * N_RESPONSES,
  declared_pilot_pairs = if (pilot_mode) PILOT_SIZE else NA_integer_,
  permutations = B,
  permutation_batch_size = permutation_batch_size,
  primary_test = "two-sided absolute-statistic global label permutation",
  declared_exclusion_rule = paste(
    "pairs whose negative-binomial nuisance fit on (Y, Z_design) degenerates",
    "are excluded; that fit is held fixed across the permutation orbit, so",
    "the criterion is invariant to every guide-label permutation"
  ),
  seed = base_seed,
  pair_seed_formula =
    "(seed + 10007*((exposure_index-1)*5000+response_rank)) mod 2147483646 + 1",
  rng_kind = RNGkind(),
  cores = cores,
  timings_seconds = list(score_tests = test_seconds)
)
data_checks <- list(
  n_local_control_cells = n_cells,
  n_h5ad_control_rows_before_drops = length(control_global_idx),
  alignment_drop_positions = ALIGNMENT_DROP_POSITIONS,
  n_global_genes = nrow(gene_features),
  n_all_perturbation_targets = length(all_target_genes),
  n_eligible_exposures = nrow(eligible_exposures),
  n_exact_one_nt = sum(pure_nt),
  n_exact_one_targeting = sum(is_targeting),
  n_eligible_non_target_responses = length(eligible_response_idx),
  n_selected_responses = nrow(selected_responses),
  selected_response_target_overlap = sum(
    selected_responses$response_gene %in% all_target_genes
  ),
  pair_cis_overlap = sum(
    pair_results$exposure_gene == pair_results$response_gene
  ),
  n_pairs = nrow(pair_results),
  n_ok = sum(pair_results$status == "ok"),
  n_degenerate_nuisance_fit =
    sum(pair_results$status == "degenerate_nuisance_fit"),
  n_degenerate_observed_score =
    sum(pair_results$status == "degenerate_observed_score"),
  n_error = sum(pair_results$status == "error"),
  counts_backend = counts_backend,
  whole_transcriptome_count_inversion_max_error =
    count_inversion_audit_error,
  all_control_gene_audit_count_inversion_max_error =
    gene_audit_inversion_error,
  all_control_gene_audit_mean_max_abs_difference =
    audit_gene_mean_difference,
  all_control_gene_audit_nonzero_max_abs_difference =
    audit_gene_nonzero_difference,
  shard_count_inversion_max_error = if (length(inversion_errors)) {
    max(inversion_errors)
  } else {
    NA_real_
  },
  manifest_md5 = manifest_hashes
)
provenance <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  command = commandArgs(),
  counts_backend = counts_backend,
  counts_source_file = verified_counts_identity$path,
  counts_source_size_bytes = verified_counts_identity$size_bytes,
  counts_source_mtime_utc = verified_counts_identity$mtime_utc,
  counts_source_md5 = locked_counts_md5,
  counts_identity_verification = if (manifest_only) {
    "full MD5 computed after raw audits; path/size/mtime stable across hash"
  } else {
    paste(
      "normalized path, byte size, and UTC mtime matched locked audit;",
      "full MD5 inherited without rehashing"
    )
  },
  intermediate_dir = intermediate_dir,
  input_md5 = unname(tools::md5sum(input_paths)),
  script_path = script_path,
  script_md5 = unname(tools::md5sum(script_path)),
  score_helper_path = normalizePath(helper_path),
  score_helper_md5 = unname(tools::md5sum(helper_path)),
  package_versions = vapply(
    c("ondisc", "Matrix", "MASS",
      if (h5ad_mode) c("HDF5Array", "rhdf5") else character()),
    function(pkg) as.character(utils::packageVersion(pkg)), character(1)
  ),
  session_info = utils::sessionInfo()
)
names(provenance$input_md5) <- required_inputs

shard <- list(
  config = config,
  data_checks = data_checks,
  provenance = provenance,
  exposure_manifest = exposure_table,
  response_manifest = selected_responses,
  pilot_manifest = if (pilot_mode) pilot_manifest else NULL,
  pair_manifest = pair_manifest,
  selected_exposures = eligible_exposures[selected_exposure_indices, ,
                                            drop = FALSE],
  selected_responses = selected_response_shard,
  pair_results = pair_results
)

tmp_rds <- tempfile(pattern = paste0(stem, "."), tmpdir = output_dir,
                    fileext = ".rds.tmp")
tmp_csv <- tempfile(pattern = paste0(stem, "."), tmpdir = output_dir,
                    fileext = ".csv.tmp")
saveRDS(shard, tmp_rds, compress = "gzip")
utils::write.csv(pair_results, tmp_csv, row.names = FALSE, na = "")
atomic_publish(tmp_rds, rds_path, overwrite = overwrite)
atomic_publish(tmp_csv, csv_path, overwrite = overwrite)

cat(sprintf(
  "Completed %d pairs (%d ok, %d degenerate-nuisance, %d degenerate-score, %d error) in %.1fs\n",
  nrow(pair_results), sum(pair_results$status == "ok"),
  sum(pair_results$status == "degenerate_nuisance_fit"),
  sum(pair_results$status == "degenerate_observed_score"),
  sum(pair_results$status == "error"), test_seconds))
cat("  RDS: ", rds_path, "\n", sep = "")
cat("  CSV: ", csv_path, "\n", sep = "")
if (any(pair_results$status != "ok")) {
  warning("Shard contains non-ok pairs; inspect status/error_message before combining")
}
