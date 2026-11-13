#!/usr/bin/env Rscript
## Synthetic contract checks only; no external data or numerical replicates.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "gb2024-preprocessing.R"))

gb_expect_error <- function(expression) {
  if (!inherits(tryCatch(force(expression), error = identity), "error")) stop("Expected a contract validation error")
}

run_gb2024_checks <- function() {
  # Columns include single NT, single targeting, multiple NT, mixed NT/target,
  # two guides for one target (including an unmeasured target), empty, and QC-fail.
  assignments <- matrix(0L, 4L, 8L)
  assignments[1L, 1L] <- assignments[3L, 2L] <- assignments[4L, 6L] <- assignments[2L, 8L] <- 1L
  assignments[c(1L, 2L), 3L] <- 1L
  assignments[c(1L, 3L), 4L] <- 1L
  assignments[c(3L, 4L), 5L] <- 1L
  flags <- data.frame(gene_passed_qc = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
                      grna_passed_qc = rep(FALSE, 8L))
  set.seed(2024)
  rng <- .Random.seed
  kind <- RNGkind()
  qc <- gb2024_cell_qc(assignments, flags, colSums(assignments))
  stopifnot(identical(which(qc$keep), c(1L, 2L, 8L)),
            identical(qc$passed_qc_column, "gene_passed_qc"))
  stopifnot(identical(which(gb2024_cell_qc(assignments)$keep), c(1L, 2L, 6L, 8L)))
  # A caller must not prefilter guide rows: dropping one guide would falsely
  # turn the two-guide, same-target column into a single-guide cell.
  stopifnot(!qc$keep[[5L]], gb2024_cell_qc(assignments[-4L, ])$keep[[5L]])
  wrong_counts <- colSums(assignments); wrong_counts[[1L]] <- 0L
  gb_expect_error(gb2024_cell_qc(assignments, stored_assignment_counts = wrong_counts))
  invalid <- assignments; invalid[1L, 1L] <- 2L
  gb_expect_error(gb2024_cell_qc(invalid))
  invalid_flags <- flags; invalid_flags$gene_passed_qc[[1L]] <- NA
  gb_expect_error(gb2024_cell_qc(assignments, invalid_flags))

  original_genes <- data.frame(n_nonzero = c(0L, 1L, 2L, 201L))
  original_guides <- data.frame(n_nonzero = c(9L, 10L, 9L, 10L),
                                target = c("A", "A", "B", "non-targeting"),
                                target_type = c("gene", "gene", "gene", "non-targeting"))
  feature_qc <- gb2024_historical_feature_qc(original_genes, original_guides, 200L)
  stopifnot(identical(feature_qc$gene_keep, c(FALSE, TRUE, TRUE, TRUE)),
            identical(feature_qc$gene_fraction[[2L]], 0.005),
            feature_qc$gene_fraction[[4L]] > 1,
            identical(feature_qc$eligible_guides, c(FALSE, TRUE, FALSE, TRUE)),
            identical(feature_qc$eligible_target_groups, "A"))
  covariates <- gb2024_covariates(c(1, 10), c(1, 2))
  stopifnot(identical(colnames(covariates), c("log_n_umis", "log_n_nonzero")),
            identical(unname(covariates), cbind(log(c(1, 10)), log(c(1, 2)))))
  gb_expect_error(gb2024_covariates(c(0, 1), c(1, 1)))
  gb_expect_error(gb2024_assert_ids(c("GENE_A", "GENE-A"), "Synthetic IDs"))
  stopifnot(identical(gb2024_assert_ids("GENE_A", "Synthetic IDs"), "GENE-A"),
            !identical(gb2024_order_hash(c("a", "b")), gb2024_order_hash(c("b", "a"))),
            identical(.Random.seed, rng), identical(RNGkind(), kind))

  # Real ondisc indexing fixture: two genes detected ONLY in excluded cells
  # must still survive the literal original-count / retained-N feature rule.
  if (!requireNamespace("ondisc", quietly = TRUE)) stop("ondisc is required for the synthetic integration fixture")
  directory <- tempfile("gb2024-preprocess-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  cells <- paste0("C", seq_len(206L))
  gene_ids <- c("G_MAIN", "G_EDGE_1", "G_EDGE_2", "G_ZERO")
  genes <- matrix(0L, 4L, 206L)
  genes[1L, ] <- 1L; genes[2L, 205L] <- 1L; genes[3L, 201L] <- 1L
  gene <- ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = genes, barcodes = cells, features_df = data.frame(gene_name = gene_ids),
    odm_fp = file.path(directory, "genes.odm"), metadata_fp = file.path(directory, "genes.rds"))
  guide_ids <- c("NT_ONE", "NT_RARE", "G_MAIN_1", "G_MAIN_2", "UNMEASURED_1", "UNMEASURED_2")
  guide_matrix <- matrix(0L, 6L, 206L)
  guide_matrix[3L, 1:150] <- 1L; guide_matrix[4L, 151:160] <- 1L
  guide_matrix[1L, 161:199] <- 1L; guide_matrix[2L, 200L] <- 1L
  guide_matrix[c(1L, 2L), 201L] <- 1L
  guide_matrix[c(1L, 3L), 202L] <- 1L
  guide_matrix[c(5L, 6L), 203L] <- 1L
  guide_matrix[5L, 204L] <- 1L; guide_matrix[3L, 205L] <- 1L
  grna <- ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = guide_matrix, barcodes = cells, features_df = data.frame(grna_name = guide_ids),
    odm_fp = file.path(directory, "guides.odm"), metadata_fp = file.path(directory, "guides.rds"))
  grna@feature_covariates$target <- c("non-targeting", "non-targeting", "G_MAIN", "G_MAIN", "UNMEASURED", "UNMEASURED")
  grna@feature_covariates$target_type <- c("non-targeting", "non-targeting", rep("gene", 4L))
  gene@cell_covariates$passed_qc <- c(rep(TRUE, 203L), FALSE, FALSE, FALSE)
  grna@cell_covariates$passed_qc <- rep(FALSE, 206L) # first matching gene flag wins
  set.seed(17)
  rng <- .Random.seed
  result <- gb2024_preprocess_modalities(gene, grna)
  stopifnot(identical(result$retained_cell_barcodes, cells[1:200]),
            identical(result$gene_features$gene_id, c("G-MAIN", "G-EDGE-1", "G-EDGE-2")),
            all(result$gene_features$detection_fraction_for_qc[2:3] == 0.005),
            identical(result$gene_features$mean_expression, gene@feature_covariates$mean_expression[1:3]),
            all(as.numeric(result$modalities$gene[["G-EDGE-1", ]]) == 0),
            nrow(result$modalities$grna) == 6L,
            identical(result$eligible_target_groups, "G_MAIN"),
            identical(result$eligible_nt_ids, "NT-ONE"),
            identical(result$all_nt_ids, c("NT-ONE", "NT-RARE")),
            result$contract$counts$retained_cells == 200L,
            result$contract$counts$pure_nt_cells == 40L,
            identical(result$contract$feature_stat_policy, "historical_stored_counts"),
            identical(result$id_mapping$grna$original_id, guide_ids),
            identical(result$id_mapping$grna$feature_id, gsub("_", "-", guide_ids, fixed = TRUE)),
            identical(result$grna_features$target, grna@feature_covariates$target),
            identical(.Random.seed, rng), identical(RNGkind(), kind))
  cat("GB2024 literal historical preprocessing tests passed; no external data or analyses used.\n")
}
run_gb2024_checks()
