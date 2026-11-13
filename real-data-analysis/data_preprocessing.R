#!/usr/bin/env Rscript
# Frangieh control-condition preprocessing matching the pinned 2024 Genome
# Biology implementation. PCIT retains its own 50-positive/950-null benchmark;
# this is not a reproduction of the paper's full association-testing benchmark.
# Usage: Rscript data_preprocessing.R 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !identical(args[[1L]], "1000")) {
  stop("Usage: Rscript data_preprocessing.R 1000")
}
num_genes <- 1000L
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "..", "reproduction-helpers.R"))
project_root <- reproduction_root(script_dir)
source(file.path(script_dir, "method_config.R"))
source(file.path(script_dir, "gb2024-preprocessing.R"))
source(file.path(script_dir, "workflow-io.R"))
results_root <- reproduction_results_root(project_root)
intermediate_dir <- Sys.getenv("FRANGIEH_INTERMEDIATE_DIR", unset =
  file.path(results_root, "frangieh-2021", PERTURB_INTERMEDIATE_SUBDIR))
intermediate_dir <- reproduction_absolute_path(intermediate_dir, project_root)
if (file.exists(intermediate_dir) &&
    (!dir.exists(intermediate_dir) || length(list.files(intermediate_dir, all.files = TRUE, no.. = TRUE)))) {
  stop("Preprocessing requires a fresh empty intermediate directory: ", intermediate_dir)
}

suppressPackageStartupMessages({ library(ondisc); library(dplyr) })

raw_root <- Sys.getenv("LOCAL_FRANGIEH_2021_DATA_DIR")
if (!nzchar(raw_root)) stop("Set LOCAL_FRANGIEH_2021_DATA_DIR to the processed Frangieh input root.")
data_dir <- file.path(path.expand(raw_root), "processed", "control")
gene_odm_fp <- file.path(data_dir, "gene", "gene_expression_matrix.odm")
gene_metadata_fp <- file.path(data_dir, "gene", "gene_expression_metadata.rds")
grna_odm_fp <- file.path(data_dir, "grna_assignment", "grna_assignments_ungrouped.odm")
grna_metadata_fp <- file.path(data_dir, "grna_assignment", "grna_assignments_ungrouped_metadata.rds")
gene_odm <- read_odm(gene_odm_fp, gene_metadata_fp)
grna_odm <- read_odm(grna_odm_fp, grna_metadata_fp)
prepared <- gb2024_preprocess_modalities(gene_odm, grna_odm)
modalities <- prepared$modalities
gene_features <- prepared$gene_features
grna_features <- prepared$grna_features
covariate_matrix <- prepared$covariate_matrix
n_cells <- ncol(modalities$gene)

# The published preprocessing determines eligibility. Ranking 50 positives and
# drawing a fixed 950-gene null universe remain PCIT-specific benchmark choices.
# Preserve the original stored mean-expression ranking and null mean threshold.
N_POS_CTRL_TOP <- 50L
NULL_MIN_MEAN_EXPR <- 0.005
all_measured_targets <- intersect(
  grna_features$target[grna_features$target_type == "gene"], gene_features$gene_id)
eligible_positive_targets <- intersect(prepared$eligible_target_groups, gene_features$gene_id)
positive_features <- gene_features |>
  filter(gene_id %in% eligible_positive_targets, mean_expression > 0) |>
  arrange(desc(mean_expression))
if (nrow(positive_features) < N_POS_CTRL_TOP) stop("Fewer than 50 eligible positive-control targets after GB2024 preprocessing")
top_pos_ctrl <- positive_features$gene_id[seq_len(N_POS_CTRL_TOP)]
non_pos_pool <- gene_features |>
  filter(!gene_id %in% all_measured_targets, mean_expression >= NULL_MIN_MEAN_EXPR) |>
  pull(gene_id)
if (length(non_pos_pool) < 950L) stop("Fewer than 950 eligible non-targeted null genes")
set.seed(2024)
sampled_genes <- sample(non_pos_pool, 950L, replace = FALSE)
gene_list <- c(top_pos_ctrl, sampled_genes)

# Retain ALL guides for an eligible target group, as in the historical grouped
# analysis; the >=10-cell rule determines target eligibility, not row deletion.
positive_control_list <- grna_features |>
  filter(target_type == "gene", target %in% top_pos_ctrl) |>
  select(grna_id, grna_target = target, grna_group)
negative_control_list <- list(
  grna = grna_features |>
    filter(grna_id %in% prepared$eligible_nt_ids) |>
    select(grna_id, grna_group),
  gene_id = sampled_genes
)
nt_grna_ids <- prepared$all_nt_ids
nt_grna_mat <- Matrix::t(modalities$grna[[nt_grna_ids, ]])
colnames(nt_grna_mat) <- nt_grna_ids
nt_grna_assigned <- nt_grna_mat >= 1
has_any_nt <- Matrix::rowSums(nt_grna_assigned) > 0

# Cell QC was evaluated on all original guide assignments, before any feature
# eligibility/normalization. A retained NT cell therefore has exactly one NT
# guide and cannot contain an unmeasured, rare, or any other targeting guide.
all_target_ids <- grna_features$grna_id[grna_features$target_type == "gene"]
targeting_grna_assigned <- Matrix::t(modalities$grna[[all_target_ids, ]]) >= 1
colnames(targeting_grna_assigned) <- all_target_ids
has_any_targeting <- Matrix::rowSums(targeting_grna_assigned) > 0
pure_nt <- as.logical(has_any_nt & !has_any_targeting)
stopifnot(all(Matrix::rowSums(nt_grna_assigned)[pure_nt] == 1L))
pos_cell_masks <- setNames(lapply(top_pos_ctrl, function(g) {
  ids <- positive_control_list$grna_id[positive_control_list$grna_target == g]
  has_target <- Matrix::rowSums(targeting_grna_assigned[, ids, drop = FALSE]) > 0
  as.logical(pure_nt | has_target)
}), top_pos_ctrl)
n_treated_per_pos <- vapply(top_pos_ctrl, function(g) {
  ids <- positive_control_list$grna_id[positive_control_list$grna_target == g]
  as.integer(sum(Matrix::rowSums(targeting_grna_assigned[, ids, drop = FALSE]) > 0))
}, integer(1))

# Retained only for compatibility with the established intermediate format.
# Production null exposures use three eligible NT guides, not these splits.
set.seed(2025)
null_splits <- lapply(seq_len(100L), function(s) {
  treated <- sample(prepared$eligible_nt_ids, 5L, replace = FALSE)
  list(treated = treated, control = setdiff(prepared$eligible_nt_ids, treated))
})
contract <- prepared$contract
contract$counts$pure_nt_cells <- as.integer(sum(pure_nt))
contract$benchmark <- list(
  design = "pcit_50_positive_950_null", positive_targets = top_pos_ctrl,
  null_gene_pool = sampled_genes, eligible_null_pool_size = length(non_pos_pool),
  positive_ranking = "original_stored_mean_expression",
  null_min_mean_expression = NULL_MIN_MEAN_EXPR,
  null_response_exclusion = "all_measured_gene_target_annotations",
  pairwise_effective_sample_size_filter = FALSE,
  note = "Published ESS>=7 pairwise QC and the complete GB2024 benchmark are not adopted in this preprocessing-only alignment."
)
contract$positive_treated_cells <- n_treated_per_pos
contract$eligible_nt_ids <- prepared$eligible_nt_ids
contract$retained_cell_barcodes <- prepared$retained_cell_barcodes
contract$id_mapping <- prepared$id_mapping
frangieh_validate_preprocessing_contract(contract)

# Publish complete immutable intermediate generations.
generations <- file.path(intermediate_dir, ".generations")
dir.create(generations, recursive = TRUE, showWarnings = FALSE)
generation <- paste0("generation-", format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"), "-", Sys.getpid())
staging <- tempfile(pattern = ".staging-", tmpdir = generations)
dir.create(staging)
objects <- list(
  subsetted_ODM = modalities, covariate_matrix = covariate_matrix,
  negative_control_pairs = negative_control_list,
  positive_control_pairs = positive_control_list, gene_list = gene_list,
  null_splits = null_splits, nt_grna_assigned = nt_grna_assigned,
  has_any_nt = has_any_nt, pure_nt = pure_nt, pos_cell_masks = pos_cell_masks,
  `preprocessing-contract` = contract
)
for (name in names(objects)) saveRDS(objects[[name]], file.path(staging, paste0(name, ".rds")))
write.csv(prepared$cohort_audit, file.path(staging, "cohort-audit.csv"), row.names = FALSE)
write.csv(data.frame(target = names(n_treated_per_pos), treated_cells = unname(n_treated_per_pos)),
          file.path(staging, "positive-target-counts.csv"), row.names = FALSE)
provenance <- list(
  schema = "frangieh_intermediate_bundle_v2",
  files = frangieh_file_manifest(file.path(staging, frangieh_intermediate_files)),
  raw_metadata = frangieh_file_manifest(c(gene_metadata_fp, grna_metadata_fp)),
  raw_backing_files = frangieh_file_manifest(c(gene_odm_fp, grna_odm_fp), hash = FALSE),
  code = frangieh_file_manifest(frangieh_preprocessing_code_files(script_dir)),
  packages = setNames(vapply(c("ondisc", "dplyr", "Matrix"), function(x) as.character(utils::packageVersion(x)), character(1)), c("ondisc", "dplyr", "Matrix")),
  session = utils::sessionInfo(), cohort_policy = PERTURB_COHORT_POLICY
)
saveRDS(provenance, file.path(staging, "intermediate-provenance.rds"))
frangieh_validate_intermediates(staging, script_dir = script_dir)
published <- file.path(generations, generation)
if (!file.rename(staging, published)) stop("Could not publish intermediate generation")
frangieh_atomic_rds(generation, file.path(intermediate_dir, "current-generation.rds"))
cat(sprintf("GB2024 preprocessing complete: %d cells, %d QC genes, %d NT controls.\n",
            n_cells, nrow(gene_features), sum(pure_nt)))
cat(sprintf("PCIT benchmark: 50 positives, 950 null-gene candidates; %d eligible NT guides.\n",
            length(prepared$eligible_nt_ids)))
cat("Intermediates: ", published, "\n", sep = "")
print(prepared$cohort_audit, row.names = FALSE)
