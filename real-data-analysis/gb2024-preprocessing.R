## Literal preprocessing contract of Barry et al., Genome Biology (2024).
## Historical manuscript code (not current SCEPTRE defaults):
## https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/005c28b9ba30e7d39c3ee18d51111722f49d32ec/R_scripts/uniform_processing_lowmoi.R#L67-L107
## https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/005c28b9ba30e7d39c3ee18d51111722f49d32ec/R_scripts/uniform_processing_lowmoi.R#L147-L180
## Covariates are log(total RNA UMIs), log(number of detected RNA features):
## the same historical file, lines 36--50; neither +1 nor batch/phase is used.
##
## IMPORTANT historical implementation detail: ondisc commit 5adcc53, called
## "v1.1.0" by the manuscript installer, has DESCRIPTION Version: 0.0.5.
## Its cell-subset operation leaves feature covariates unchanged:
## https://github.com/timothy-barry/ondisc/blob/5adcc537664aada08c3e79f08f2a944a0d31598d/R/s4_subset_extract.R#L118-L128
## Its feature-QC helper divides the ORIGINAL stored n_nonzero by the CURRENT
## retained-cell count and applies >=, without recomputing the numerator:
## https://github.com/timothy-barry/ondisc/blob/5adcc537664aada08c3e79f08f2a944a0d31598d/R/qc_helpers.R#L49-L55
## Accordingly this helper deliberately implements historical_stored_counts.
## It must not be described as recomputed detection within retained cells.
## No function here draws random numbers or changes an RNG configuration.

gb2024_assert_ids <- function(ids, label) {
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " must be unique, non-missing character identifiers")
  }
  normalized <- gsub("_", "-", ids, fixed = TRUE)
  if (anyDuplicated(normalized)) stop(label, " collide after historical '_' to '-' normalization")
  normalized
}

gb2024_order_hash <- function(ids) {
  path <- tempfile("gb2024-order-")
  on.exit(unlink(path), add = TRUE)
  writeBin(serialize(enc2utf8(as.character(ids)), NULL, version = 2L), path)
  unname(tools::md5sum(path))
}

## guide_assignments must include ALL original guide rows, before any target,
## guide-frequency, or response-feature selection. Entries are original binary
## author assignments, not raw guide UMIs to be thresholded anew.
gb2024_cell_qc <- function(guide_assignments, cell_covariates = NULL,
                           stored_assignment_counts = NULL) {
  dimensions <- dim(guide_assignments)
  if (length(dimensions) != 2L || any(dimensions < 1L)) stop("Guide assignments must be a nonempty matrix")
  values <- if (inherits(guide_assignments, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(guide_assignments)) guide_assignments@x else 1L
  } else as.vector(guide_assignments)
  if (anyNA(values) || any(!is.finite(values)) || any(!values %in% c(0, 1))) {
    stop("Use original binary guide assignments; do not pass raw guide UMI counts")
  }
  presence <- guide_assignments != 0
  n_assigned <- as.integer(if (inherits(presence, "Matrix")) Matrix::colSums(presence) else colSums(presence))
  if (!is.null(stored_assignment_counts) &&
      (!is.numeric(stored_assignment_counts) || length(stored_assignment_counts) != dimensions[[2L]] ||
       anyNA(stored_assignment_counts) || any(stored_assignment_counts != n_assigned))) {
    stop("Original guide assignment counts disagree with the all-guide assignment matrix")
  }
  single <- n_assigned == 1L
  passed <- rep(TRUE, dimensions[[2L]])
  passed_qc_column <- "none"
  if (!is.null(cell_covariates)) {
    if (!is.data.frame(cell_covariates) || nrow(cell_covariates) != dimensions[[2L]]) stop("Cell covariates must align with guide columns")
    candidates <- which(grepl("passed_qc", names(cell_covariates)))
    if (length(candidates)) {
      passed_qc_column <- names(cell_covariates)[candidates[[1L]]]
      passed <- cell_covariates[[candidates[[1L]]]]
      if ((!is.logical(passed) && !is.numeric(passed)) || anyNA(passed) ||
          any(!passed %in% c(FALSE, TRUE))) stop("The original passed_qc flag must contain non-missing binary values")
      passed <- as.logical(passed)
    }
  }
  keep <- single & passed
  if (!any(keep)) stop("No cells pass historical single-guide/original-author QC")
  entries <- if (inherits(presence, "Matrix")) Matrix::which(presence, arr.ind = TRUE) else which(presence, arr.ind = TRUE)
  assigned_row <- integer(dimensions[[2L]])
  assigned_row[entries[, 2L]] <- entries[, 1L]
  list(keep = keep, n_assigned = n_assigned, assigned_row = assigned_row[keep],
       passed_qc_column = passed_qc_column,
       audit = data.frame(stage = c("original_cells", "exactly_one_original_guide",
                                    "original_author_qc", "retained_cells"),
                          n_cells = c(dimensions[[2L]], sum(single), sum(passed), sum(keep))))
}

gb2024_historical_feature_qc <- function(gene_features, grna_features, n_retained_cells) {
  if (length(n_retained_cells) != 1L || !is.finite(n_retained_cells) ||
      n_retained_cells < 1L || n_retained_cells != as.integer(n_retained_cells)) stop("Invalid retained-cell count")
  for (features in list(gene_features, grna_features)) {
    if (!is.data.frame(features) || !"n_nonzero" %in% names(features) ||
        !is.numeric(features$n_nonzero) || anyNA(features$n_nonzero) ||
        any(!is.finite(features$n_nonzero)) || any(features$n_nonzero < 0) ||
        any(features$n_nonzero != as.integer(features$n_nonzero))) stop("Invalid original stored feature detection counts")
  }
  if (!all(c("target", "target_type") %in% names(grna_features))) stop("Missing original guide annotations")
  ## Values above one are possible: the numerator is intentionally PRE-QC.
  gene_fraction <- gene_features$n_nonzero / n_retained_cells
  eligible_guides <- grna_features$n_nonzero >= 10L
  target <- as.character(grna_features$target)
  list(gene_keep = gene_fraction >= 0.005, gene_fraction = gene_fraction,
       eligible_guides = eligible_guides,
       eligible_target_groups = unique(target[eligible_guides & target != "non-targeting"]))
}

gb2024_covariates <- function(n_umis, n_nonzero, cell_barcodes = NULL) {
  if (!is.numeric(n_umis) || !is.numeric(n_nonzero) || !length(n_umis) ||
      length(n_umis) != length(n_nonzero) || anyNA(n_umis) || anyNA(n_nonzero) ||
      any(!is.finite(n_umis)) || any(!is.finite(n_nonzero)) ||
      any(n_umis <= 0) || any(n_nonzero <= 0)) {
    stop("Historical log covariates require positive full-transcriptome n_umis and n_nonzero; no extra cell filter or +1 is applied")
  }
  result <- cbind(log_n_umis = log(n_umis), log_n_nonzero = log(n_nonzero))
  if (!is.null(cell_barcodes)) {
    if (length(cell_barcodes) != nrow(result)) stop("Covariate cell identifiers do not align")
    rownames(result) <- cell_barcodes
  }
  result
}

gb2024_preprocess_modalities <- function(gene_odm, grna_odm) {
  if (!requireNamespace("ondisc", quietly = TRUE)) stop("ondisc is required")
  original_cells <- ondisc::get_cell_barcodes(gene_odm)
  if (!identical(original_cells, ondisc::get_cell_barcodes(grna_odm))) stop("Gene and guide cell orders must agree before historical QC")
  original_gene_ids <- ondisc::get_feature_ids(gene_odm)
  original_grna_ids <- ondisc::get_feature_ids(grna_odm)
  gene_ids <- gb2024_assert_ids(original_gene_ids, "Gene IDs")
  grna_ids <- gb2024_assert_ids(original_grna_ids, "Guide IDs")
  gene_features <- ondisc::get_feature_covariates(gene_odm)
  grna_features <- ondisc::get_feature_covariates(grna_odm)
  if (!identical(rownames(gene_features), original_gene_ids) ||
      !identical(rownames(grna_features), original_grna_ids)) stop("Feature metadata and matrix identifiers do not align")
  if (!"mean_expression" %in% names(gene_features) || anyNA(gene_features$mean_expression) ||
      any(!is.finite(gene_features$mean_expression)) || any(gene_features$mean_expression < 0)) stop("Missing original mean-expression metadata")
  if (!all(c("target", "target_type") %in% names(grna_features)) ||
      anyNA(grna_features$target) || anyNA(grna_features$target_type) ||
      any(!as.character(grna_features$target_type) %in% c("gene", "non-targeting"))) stop("Invalid original guide annotations")
  mm <- ondisc::multimodal_ondisc_matrix(list(gene = gene_odm, grna = grna_odm))
  global_cell_covariates <- ondisc::get_cell_covariates(mm)
  grna_cells <- ondisc::get_cell_covariates(grna_odm)
  if (!"n_nonzero" %in% names(grna_cells)) stop("Missing original assignment-count metadata")
  assignments <- grna_odm[[seq_along(original_grna_ids), ]]
  cells <- gb2024_cell_qc(assignments, global_cell_covariates, grna_cells$n_nonzero)
  retained_cells <- original_cells[cells$keep]
  features <- gb2024_historical_feature_qc(gene_features, grna_features, length(retained_cells))
  if (!any(features$gene_keep)) stop("No responses pass historical feature QC")

  ## Keep ALL guide rows even when a guide has fewer than 10 assignments. That
  ## threshold selects eligible target groups/pair lists, not guide-row deletion.
  gene <- gene_odm[features$gene_keep, cells$keep]
  grna <- grna_odm[, cells$keep]
  gene@ondisc_matrix@feature_ids <- gsub("_", "-", gene@ondisc_matrix@feature_ids, fixed = TRUE)
  grna@ondisc_matrix@feature_ids <- gsub("_", "-", grna@ondisc_matrix@feature_ids, fixed = TRUE)
  rownames(gene@feature_covariates) <- gene_ids[features$gene_keep]
  rownames(grna@feature_covariates) <- grna_ids
  assigned_grna_id <- grna_ids[cells$assigned_row]
  grna@cell_covariates$assigned_grna <- assigned_grna_id
  gene_features <- gene_features[features$gene_keep, , drop = FALSE]
  rownames(gene_features) <- gene_ids[features$gene_keep]
  gene_features$gene_id <- gene_ids[features$gene_keep]
  gene_features$detection_fraction_for_qc <- features$gene_fraction[features$gene_keep]
  rownames(grna_features) <- grna_ids
  grna_features$grna_id <- grna_ids
  grna_features$target <- as.character(grna_features$target)
  grna_features$target_type <- as.character(grna_features$target_type)
  ## NT guide IDs remain individually addressable for PCIT's null construction.
  grna_features$grna_group <- ifelse(grna_features$target_type == "non-targeting",
                                     grna_ids, grna_features$target)
  gene_cells <- ondisc::get_cell_covariates(gene_odm)
  if (!all(c("n_umis", "n_nonzero") %in% names(gene_cells))) stop("Missing original full-transcriptome cell covariates")
  covariates <- gb2024_covariates(gene_cells$n_umis[cells$keep],
                                 gene_cells$n_nonzero[cells$keep], retained_cells)
  is_nt <- grna_features$target_type == "non-targeting"
  assigned_target <- grna_features$target[cells$assigned_row]
  assigned_target_type <- grna_features$target_type[cells$assigned_row]
  contract <- list(
    schema = "frangieh_gb2024_preprocessing_v1",
    policy = "gb2024_single_guide_historical_stored_counts_v1",
    feature_stat_policy = "historical_stored_counts",
    cell_rule = "exactly_one_original_assigned_guide_before_feature_qc",
    reference_commit = "005c28b9ba30e7d39c3ee18d51111722f49d32ec",
    ondisc_commit = "5adcc537664aada08c3e79f08f2a944a0d31598d",
    gene_fraction_min = 0.005, guide_cells_min = 10L,
    gene_statistic = "original_metadata_n_nonzero_divided_by_retained_cell_count",
    guide_rule = "original_metadata_count_selects_target_groups_keep_all_guides_within_group",
    covariates = colnames(covariates),
    passed_qc_column = cells$passed_qc_column, id_cleanup = "underscore_to_dash",
    counts = list(input_cells = length(original_cells), retained_cells = length(retained_cells),
                  pure_nt_cells = sum(assigned_target_type == "non-targeting"),
                  input_genes = length(original_gene_ids), retained_genes = nrow(gene_features),
                  input_guides = length(original_grna_ids), eligible_guides = sum(features$eligible_guides),
                  eligible_target_groups = length(features$eligible_target_groups)),
    order_hashes = c(cells = gb2024_order_hash(retained_cells),
                     genes = gb2024_order_hash(gene_features$gene_id),
                     guides = gb2024_order_hash(grna_ids)))
  list(modalities = list(gene = gene, grna = grna), gene_features = gene_features,
       grna_features = grna_features, eligible_target_groups = features$eligible_target_groups,
       eligible_nt_ids = grna_ids[features$eligible_guides & is_nt], all_nt_ids = grna_ids[is_nt],
       covariate_matrix = covariates, contract = contract, cohort_audit = cells$audit,
       retained_cell_barcodes = retained_cells, assigned_grna_id = assigned_grna_id,
       assigned_target = assigned_target, assigned_target_type = assigned_target_type,
       id_mapping = list(gene = data.frame(original_id = original_gene_ids, feature_id = gene_ids,
                                           retained = features$gene_keep),
                         grna = data.frame(original_id = original_grna_ids, feature_id = grna_ids)))
}
