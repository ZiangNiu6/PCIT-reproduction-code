#!/usr/bin/env Rscript
# data_preprocessing.R
#
# Preprocess Frangieh 2021 (control condition) for the PCIT reproduction.
# Loads gene expression and gRNA assignment from ondisc format,
# selects 1000 lowly-expressed genes, builds covariate matrix,
# and saves intermediate data for downstream analysis.
#
# Usage: Rscript data_preprocessing.R <num_genes>
#   e.g., Rscript data_preprocessing.R 1000

library(ondisc)
library(dplyr)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
force_rerun <- "--force" %in% args
args <- args[args != "--force"]
if (length(args) != 1L) {
  stop("Usage: Rscript data_preprocessing.R <num_genes> [--force]")
}
num_genes <- as.numeric(args[[1L]])
if (length(num_genes) != 1L || !is.finite(num_genes) ||
    num_genes != as.integer(num_genes) || num_genes < 51L) {
  stop("num_genes must be an integer of at least 51")
}
num_genes <- as.integer(num_genes)

# --- Get script directory ---
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_all[grep("--file=", args_all)]))
script_dir <- dirname(script_path)
source(file.path(script_dir, "..", "reproduction-helpers.R"))
project_root <- reproduction_root(script_dir)
results_root <- reproduction_results_root(project_root)
intermediate_dir <- Sys.getenv(
  "FRANGIEH_INTERMEDIATE_DIR",
  unset = file.path(results_root, "frangieh-2021", "intermediate-data")
)

source(file.path(script_dir, "workflow-io.R"))
if (num_genes != 1000L) stop("The canonical reproduction requires num_genes=1000")
existing <- file.exists(file.path(intermediate_dir, "current-generation.rds")) ||
  any(file.exists(file.path(intermediate_dir, frangieh_intermediate_files)))
if (existing && !force_rerun) {
  frangieh_validate_intermediates(intermediate_dir, script_dir = script_dir)
  cat("Reusing the complete validated intermediate bundle.\n")
  quit(save = "no", status = 0L)
}

# --- Paths ---
frangieh_dir <- Sys.getenv("LOCAL_FRANGIEH_2021_DATA_DIR")
if (frangieh_dir == "") {
  stop(
    "Set LOCAL_FRANGIEH_2021_DATA_DIR to the Frangieh 2021 data root."
  )
}
data_dir <- file.path(frangieh_dir, "processed", "control")

gene_odm_fp <- file.path(data_dir, "gene", "gene_expression_matrix.odm")
gene_metadata_fp <- file.path(data_dir, "gene", "gene_expression_metadata.rds")
grna_odm_fp <- file.path(data_dir, "grna_assignment", "grna_assignments_ungrouped.odm")
grna_metadata_fp <- file.path(data_dir, "grna_assignment", "grna_assignments_ungrouped_metadata.rds")

cat("Loading ondisc matrices...\n")
gene_odm <- read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)
grna_odm <- read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)

cat(sprintf("  Gene: %d genes x %d cells\n", nrow(gene_odm), ncol(gene_odm)))
cat(sprintf("  gRNA: %d gRNAs x %d cells\n", nrow(grna_odm), ncol(grna_odm)))

# --- Filter cells with zero UMIs ---
multimodal_odm <- multimodal_ondisc_matrix(
  covariate_ondisc_matrix_list = list(grna = grna_odm, gene = gene_odm)
)
# gRNA assignment matrix has no n_umis (binary); filter on gene expression only
ok_cells_gene <- (gene_odm |> ondisc::get_cell_covariates() |> pull(n_umis)) != 0L
multimodal_odm_sub <- multimodal_odm[, ok_cells_gene]

n_cells <- ncol(multimodal_odm_sub@modalities$gene)
cat(sprintf("After filtering: %d cells\n", n_cells))

# --- Build gRNA target data frame ---
grna_feature_df <- grna_odm |>
  ondisc::get_feature_covariates() |>
  tibble::rownames_to_column(var = "grna_id") |>
  select(grna_id, target, target_type) |>
  tibble::as_tibble()

# Create grna_group: for gene-targeting, group by target; for non-targeting, use grna_id
grna_feature_df <- grna_feature_df |>
  mutate(grna_group = if_else(target_type == "non-targeting", grna_id, target))

# --- Identify non-targeting gRNA groups ---
nt_grna_groups <- grna_feature_df |>
  filter(target_type == "non-targeting")
cat(sprintf("Non-targeting gRNAs: %d\n", nrow(nt_grna_groups)))

# --- Identify positive control pairs ---
gene_full_list <- multimodal_odm_sub@modalities$gene |> ondisc::get_feature_ids()
targeting_grnas <- grna_feature_df |>
  filter(target_type == "gene", target %in% gene_full_list)
positive_control_targets <- unique(targeting_grnas$target)
cat(sprintf("Positive control gene targets: %d\n", length(positive_control_targets)))

# --- Select top 50 highly-expressed positive controls + non-positive null pool ---
# Design: alts = top N_POS_CTRL_TOP positive controls by mean expression.
#         nulls (random gene x random NT-gRNA in run_combined_fdr.R) drawn
#         from a sampled pool of non-positive-control nonzero genes.
# Total gene panel = top positives + sampled non-positive genes = num_genes.
N_POS_CTRL_TOP <- 50L

gene_feat <- gene_odm |> ondisc::get_feature_covariates()
gene_feat$gene_id <- rownames(gene_feat)

# Filter to genes with non-zero expression
nonzero_genes <- gene_feat |> filter(mean_expression > 0)
cat(sprintf("Genes with non-zero expression: %d\n", nrow(nonzero_genes)))

# Top-N_POS_CTRL_TOP positive controls by mean expression
pos_ctrl_present <- positive_control_targets[positive_control_targets %in% nonzero_genes$gene_id]
cat(sprintf("Positive control targets present in expression matrix: %d\n",
            length(pos_ctrl_present)))
if (length(pos_ctrl_present) < N_POS_CTRL_TOP) {
  stop(sprintf("Need at least %d positive controls in data; have %d.",
               N_POS_CTRL_TOP, length(pos_ctrl_present)))
}

pos_ctrl_with_expr <- nonzero_genes |>
  filter(gene_id %in% pos_ctrl_present) |>
  arrange(desc(mean_expression))
top_pos_ctrl <- pos_ctrl_with_expr$gene_id[seq_len(N_POS_CTRL_TOP)]
cat(sprintf("Selected top %d positive controls by mean expression\n", N_POS_CTRL_TOP))
cat(sprintf("  Top-1 mean expression:    %.3f\n",
            pos_ctrl_with_expr$mean_expression[1]))
cat(sprintf("  Top-%d mean expression:   %.3f\n", N_POS_CTRL_TOP,
            pos_ctrl_with_expr$mean_expression[N_POS_CTRL_TOP]))

# Sample the remaining (num_genes - top positives) slots from non-positive-control
# genes whose mean_expression >= NULL_MIN_MEAN_EXPR.  These form the universe
# of null gene candidates in the combined experiment (filtered to avoid
# ultra-sparse genes where score-test calibration breaks down).
NULL_MIN_MEAN_EXPR <- 0.005
set.seed(2024)
non_pos_pool <- nonzero_genes |>
  filter(!gene_id %in% positive_control_targets,
         mean_expression >= NULL_MIN_MEAN_EXPR) |>
  pull(gene_id)
cat(sprintf("Non-positive null pool with mean_expression >= %.3f: %d genes\n",
            NULL_MIN_MEAN_EXPR, length(non_pos_pool)))
n_to_sample  <- num_genes - N_POS_CTRL_TOP
if (n_to_sample > length(non_pos_pool)) {
  stop(sprintf("Non-positive null pool has %d genes; need %d.  Lower NULL_MIN_MEAN_EXPR.",
               length(non_pos_pool), n_to_sample))
}
sampled_genes <- sample(non_pos_pool, n_to_sample, replace = FALSE)
gene_list <- c(top_pos_ctrl, sampled_genes)

cat(sprintf("Selected %d genes total (%d top-expressed positives + %d sampled non-positives)\n",
            length(gene_list), N_POS_CTRL_TOP, length(sampled_genes)))

# --- (Gene and gRNA matrices accessed via subsetted_ODM in downstream scripts) ---

# --- Build covariate matrix directly ---
cell_cov_df <- multimodal_odm_sub |>
  ondisc::get_cell_covariates()

# Build design matrix: batch (factor -> dummies), log(n_umis)
batch_vec <- cell_cov_df$gene_batch
n_umis <- cell_cov_df$gene_n_umis

# Construct covariate data frame
cov_df <- data.frame(
  log_n_umis = log(n_umis + 1),
  batch = factor(batch_vec)
)

# Create design matrix (drop intercept)
covariate_matrix <- model.matrix(~ log_n_umis + batch, data = cov_df)[, -1, drop = FALSE]
cat(sprintf("Covariate matrix: %d cells x %d covariates\n",
            nrow(covariate_matrix), ncol(covariate_matrix)))
cat(sprintf("  Columns: %s\n", paste(colnames(covariate_matrix), collapse = ", ")))

# Positive control pairs (only for genes in our selected set)
positive_control_grna_df <- targeting_grnas |>
  filter(target %in% gene_list)

# --- Build negative control list ---
negative_control_list <- list(
  grna = nt_grna_groups |> select(grna_id, grna_group),
  ## Only non-positive genes are eligible as synthetic-null responses.
  gene_id = sampled_genes
)

# --- Build positive control list ---
positive_control_list <- positive_control_grna_df |>
  select(grna_id, grna_target = target, grna_group)

# --- Generate 100 random null splits of NT gRNAs ---
# Each split: 5 "treated" NT gRNAs vs 69 "control" NT gRNAs
N_SPLITS <- 100L
N_TREATED <- 5L
nt_grna_ids <- nt_grna_groups$grna_id

set.seed(2025)
null_splits <- lapply(seq_len(N_SPLITS), function(s) {
  treated <- sample(nt_grna_ids, N_TREATED, replace = FALSE)
  control <- setdiff(nt_grna_ids, treated)
  list(treated = treated, control = control)
})
cat(sprintf("Generated %d null splits (%d treated vs %d control NT gRNAs each)\n",
            N_SPLITS, N_TREATED, length(nt_grna_ids) - N_TREATED))

# --- Pre-compute NT gRNA assignment matrix (binary, cells x 74 NT gRNAs) ---
nt_grna_mat <- multimodal_odm_sub@modalities$grna[[nt_grna_ids, ]] |>
  Matrix::t()
colnames(nt_grna_mat) <- nt_grna_ids
# gRNA assignment matrix is already binary (0/1); cell is assigned if value >= 1
nt_grna_assigned <- (nt_grna_mat >= 1)

# Identify cells assigned to at least one NT gRNA
has_any_nt <- Matrix::rowSums(nt_grna_assigned) > 0
cat(sprintf("Cells with at least one NT gRNA: %d / %d\n",
            sum(has_any_nt), length(has_any_nt)))

# --- Build per-positive-control cell mask + pure-NT control set --------------
# Original PCIT v2 cohort retained by explicit user instruction. This includes
# multi-guide cells and excludes only measured-target guides from pure_nt; it
# intentionally does not adopt the single-guide cohort used by Barry et al.
# A scientific correction requires a separate reviewed numerical rerun.
# Historical PCIT cell-set construction:
#   * Historical pure_nt = cells carrying NT gRNAs and NO measured-target
#                          gene-targeting gRNAs.
#                         These are the universal control population.
#   * Per-positive-g    = pure-NT cells (controls) UNION cells carrying any
#                         gRNA targeting g (treated).
# For null pairs (NT-subset X), we use the pure-NT cell set.
# For positive pair g, we use pos_cell_masks[[g]].
all_targeting_grna_ids <- targeting_grnas$grna_id
cat(sprintf("Loading targeting-gRNA assignment matrix (%d targeting gRNAs)...\n",
            length(all_targeting_grna_ids)))
targeting_grna_mat <- multimodal_odm_sub@modalities$grna[[all_targeting_grna_ids, ]] |>
  Matrix::t()
colnames(targeting_grna_mat) <- all_targeting_grna_ids
targeting_grna_assigned <- (targeting_grna_mat >= 1)

has_any_targeting <- Matrix::rowSums(targeting_grna_assigned) > 0
pure_nt           <- has_any_nt & !has_any_targeting
cat(sprintf("Pure-NT cells (NT only, no targeting): %d / %d\n",
            sum(pure_nt), length(pure_nt)))

# Per-positive-control cell mask: pure-NT controls + cells targeting g
pos_cell_masks <- vector("list", length(top_pos_ctrl))
names(pos_cell_masks) <- top_pos_ctrl
for (g in top_pos_ctrl) {
  g_grna_ids <- positive_control_grna_df |>
    filter(target == !!g) |> pull(grna_id)
  has_g_targeting <- Matrix::rowSums(targeting_grna_assigned[, g_grna_ids, drop = FALSE]) > 0
  pos_cell_masks[[g]] <- pure_nt | has_g_targeting
}
n_treated_per_pos <- vapply(top_pos_ctrl, function(g) {
  g_grna_ids <- positive_control_grna_df |>
    filter(target == !!g) |> pull(grna_id)
  as.integer(sum(Matrix::rowSums(targeting_grna_assigned[, g_grna_ids, drop = FALSE]) > 0))
}, integer(1))
cat(sprintf("Cells per positive-control pair (mean +/- SD): %.0f +/- %.0f (treated: %.0f +/- %.0f)\n",
            mean(vapply(pos_cell_masks, sum, integer(1))),
            sd  (vapply(pos_cell_masks, sum, integer(1))),
            mean(n_treated_per_pos), sd(n_treated_per_pos)))

# --- Publish an immutable intermediate bundle ------------------------------
# A same-filesystem pointer rename publishes all ten files together. A reader
# already using an older generation keeps that complete immutable generation.
generations <- file.path(intermediate_dir, ".generations")
dir.create(generations, showWarnings = FALSE, recursive = TRUE)
generation <- paste0("generation-", format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"), "-", Sys.getpid())
staging <- tempfile(pattern = ".staging-", tmpdir = generations)
dir.create(staging)
objects <- list(
  subsetted_ODM = multimodal_odm_sub@modalities,
  covariate_matrix = covariate_matrix,
  negative_control_pairs = negative_control_list,
  positive_control_pairs = positive_control_list,
  gene_list = gene_list, null_splits = null_splits,
  nt_grna_assigned = nt_grna_assigned, has_any_nt = has_any_nt,
  pure_nt = pure_nt, pos_cell_masks = pos_cell_masks
)
for (name in names(objects)) saveRDS(objects[[name]], file.path(staging, paste0(name, ".rds")))
provenance <- list(
  schema = "frangieh_intermediate_bundle_v1",
  files = frangieh_file_manifest(file.path(staging, frangieh_intermediate_files)),
  raw_metadata = frangieh_file_manifest(c(gene_metadata_fp, grna_metadata_fp)),
  raw_backing_files = frangieh_file_manifest(c(gene_odm_fp, grna_odm_fp), hash = FALSE),
  code = frangieh_file_manifest(c(file.path(script_dir, "data_preprocessing.R"), file.path(script_dir, "workflow-io.R"))),
  packages = setNames(vapply(c("ondisc", "dplyr", "Matrix"), function(x) as.character(utils::packageVersion(x)), character(1)), c("ondisc", "dplyr", "Matrix")),
  session = utils::sessionInfo(),
  cohort_policy = "original_v2_cohort_retained_without_scientific_changes"
)
saveRDS(provenance, file.path(staging, "intermediate-provenance.rds"))
frangieh_validate_intermediates(staging, script_dir = script_dir)
published <- file.path(generations, generation)
if (!file.rename(staging, published)) stop("Could not publish intermediate generation")
frangieh_atomic_rds(generation, file.path(intermediate_dir, "current-generation.rds"), replace = force_rerun)
intermediate_dir <- published

cat(sprintf("\nIntermediate data saved to %s\n", intermediate_dir))
cat(sprintf("  subsetted_ODM.rds\n"))
cat(sprintf("  covariate_matrix.rds (%d x %d)\n", nrow(covariate_matrix), ncol(covariate_matrix)))
cat(sprintf("  negative_control_pairs.rds (%d NT gRNAs, %d genes)\n",
            nrow(nt_grna_groups), length(negative_control_list$gene_id)))
cat(sprintf("  positive_control_pairs.rds (%d pairs)\n", nrow(positive_control_list)))
cat(sprintf("  gene_list.rds (%d genes)\n", length(gene_list)))
cat(sprintf("  null_splits.rds (%d splits, %d treated per split)\n", N_SPLITS, N_TREATED))
cat(sprintf("  nt_grna_assigned.rds (%d cells x %d NT gRNAs)\n",
            nrow(nt_grna_assigned), ncol(nt_grna_assigned)))
cat("Done!\n")
