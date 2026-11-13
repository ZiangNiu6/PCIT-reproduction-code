#!/usr/bin/env Rscript
# Build a plain-R cache for a separate q=.10 run of the exact canonical pairs.
# No canonical input/output is changed. This script does not draw resamples.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
suppressPackageStartupMessages({library(ondisc); library(Matrix); library(digest); library(jsonlite)})

main <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(script_arg) == 1L)
  script <- normalizePath(sub("^--file=", "", script_arg))
  repo <- normalizePath(file.path(dirname(script), "../../.."))
  options <- list(
    `results-dir` = file.path(repo, "results/frangieh-q10"),
    `canonical-dir` = file.path(repo, "reproduction-code/results/frangieh-2021/score-calibration-gb2024-single-guide-v3"),
    `intermediate-dir` = file.path(repo, "reproduction-code/results/frangieh-2021/intermediate-data-gb2024-single-guide-v3"),
    `effects-cache` = file.path(repo, "results/frangieh-main-effects/effect-inputs.rds"),
    `data-root` = ""
  )
  args <- commandArgs(TRUE)
  while (length(args)) {
    key <- sub("^--", "", args[[1L]])
    if (length(args) < 2L || !key %in% names(options)) stop("Expected one of --", paste(names(options), collapse = ", --"), " followed by a path")
    options[[key]] <- args[[2L]]
    args <- args[-c(1L, 2L)]
  }
  canonical <- normalizePath(options[["canonical-dir"]])
  intermediate <- normalizePath(options[["intermediate-dir"]])
  audit_path <- file.path(canonical, "input-audit.json")
  audit <- fromJSON(audit_path, simplifyVector = FALSE)
  raw_root <- normalizePath(if (nzchar(options[["data-root"]])) options[["data-root"]] else audit$root)
  out <- options[["results-dir"]]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out)
  stopifnot(!startsWith(out, paste0(canonical, "/")), out != canonical)
  lock <- file.path(out, ".cache-preparation.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Cache preparation already locked: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  cache_path <- file.path(out, "cache.rds")
  manifest_path <- file.path(out, "cache-manifest.json")
  sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  record <- function(path) {
    before <- file.info(path)
    if (!is.finite(before$size) || before$isdir) stop("Missing file: ", path)
    digest <- sha(path)
    after <- file.info(path)
    stopifnot(identical(before$size, after$size), identical(before$mtime, after$mtime))
    list(path = normalizePath(path), bytes = as.numeric(after$size), sha256 = digest,
         mtime_utc = format(after$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  }
  records <- function(paths) setNames(lapply(paths, record), names(paths))
  workflow <- file.path(repo, "reproduction-code/real-data-analysis")
  Sys.setenv(FRANGIEH_ANALYSIS_SCHEMA = "gb2024", LOCAL_FRANGIEH_2021_DATA_DIR = raw_root)
  source(file.path(workflow, "method_config.R"))
  source(file.path(workflow, "workflow-io.R"))
  source(file.path(workflow, "score_helpers.R"))
  generation <- frangieh_validate_intermediates(intermediate, script_dir = workflow)
  contract <- frangieh_read_preprocessing_contract(generation)
  pointer_path <- file.path(intermediate, "current-generation.rds")
  rep_paths <- setNames(file.path(canonical, sprintf("score_fdr_rep%03d.rds", 1:20)), sprintf("rep%03d", 1:20))
  intermediate_paths <- setNames(list.files(generation, full.names = TRUE), list.files(generation))
  code_paths <- c(preparer = script, setNames(file.path(workflow, c("score_helpers.R", "method_config.R", "workflow-io.R", "gb2024-preprocessing.R", "run_combined_fdr.R", "result-validation.R")),
                                               c("score_helpers.R", "method_config.R", "workflow-io.R", "gb2024-preprocessing.R", "run_combined_fdr.R", "result-validation.R")))
  raw_paths <- setNames(file.path(raw_root, vapply(audit$files, `[[`, character(1), "path")), vapply(audit$files, `[[`, character(1), "path"))
  cat("Checking original raw input hashes and canonical replicate receipts...\n")
  raw_records <- records(raw_paths)
  for (entry in audit$files) stopifnot(raw_records[[entry$path]]$bytes == entry$bytes,
                                      identical(raw_records[[entry$path]]$sha256, entry$sha256))
  effects_path <- normalizePath(options[["effects-cache"]])
  effects_manifest_path <- file.path(dirname(effects_path), "effect-inputs-manifest.json")
  effects_manifest <- fromJSON(effects_manifest_path, simplifyVector = FALSE)
  stopifnot(identical(effects_manifest$status, "complete"), identical(sha(effects_path), effects_manifest$output$sha256))
  config <- list(schema = "frangieh_q10_plain_score_cache_v1", source_schema = PERTURB_SCHEMA_VERSION,
                 pointer = record(pointer_path), intermediate_inputs = records(intermediate_paths),
                 canonical_replicates = records(rep_paths), input_audit = record(audit_path),
                 raw_inputs = raw_records, effects_cache = record(effects_path),
                 effects_manifest = record(effects_manifest_path), code = records(code_paths),
                 score_absolute_tolerance = 1e-7,
                 packages = as.list(setNames(vapply(c("MASS", "Matrix", "ondisc", "pcmt"), function(p) as.character(packageVersion(p)), character(1)), c("MASS", "Matrix", "ondisc", "pcmt"))),
                 R = as.character(getRversion()))
  config_hash <- digest(toJSON(config, auto_unbox = TRUE, digits = NA, null = "null"), algo = "sha256", serialize = FALSE)
  if (file.exists(cache_path) || file.exists(manifest_path)) {
    if (!file.exists(cache_path) || !file.exists(manifest_path)) stop("Incomplete existing cache; choose a new output directory")
    old <- fromJSON(manifest_path, simplifyVector = FALSE)
    stopifnot(identical(old$status, "complete"), identical(old$config_sha256, config_hash),
              old$output$bytes == file.info(cache_path)$size, identical(old$output$sha256, sha(cache_path)))
    cat("Verified existing immutable plain-R cache without rewriting it.\n")
    return(invisible(TRUE))
  }
  started <- Sys.time()
  objects <- setNames(lapply(file.path(generation, frangieh_intermediate_files), readRDS), frangieh_intermediate_files)
  modalities <- frangieh_rebind_modalities(objects[["subsetted_ODM.rds"]], raw_root)
  covariates <- objects[["covariate_matrix.rds"]]
  Z_design_full <- cbind(intercept = 1, as.matrix(covariates))
  pure_nt <- objects[["pure_nt.rds"]]
  pure_nt_idx <- which(pure_nt)
  Z_design_pure_nt <- Z_design_full[pure_nt_idx, , drop = FALSE]
  positive <- objects[["positive_control_pairs.rds"]]
  pos_targets <- unique(positive$grna_target)
  masks <- objects[["pos_cell_masks.rds"]]
  cell_ids <- get_cell_barcodes(modalities$gene)
  nt_grna_ids <- objects[["negative_control_pairs.rds"]]$grna$grna_id
  null_genes <- perturb_null_gene_pool(objects[["gene_list.rds"]], objects[["negative_control_pairs.rds"]]$gene_id, pos_targets, 950L)
  stopifnot(length(cell_ids) == 30486L, length(pure_nt_idx) == 3770L, length(pos_targets) == 50L,
            identical(cell_ids, rownames(covariates)), identical(cell_ids, get_cell_barcodes(modalities$grna)),
            identical(contract$retained_cell_barcodes, cell_ids), length(nt_grna_ids) == 73L,
            !anyDuplicated(cell_ids), !anyDuplicated(pos_targets))
  nt_mat_full <- as.matrix(Matrix::t(modalities$grna[[nt_grna_ids, ]]))
  nt_mat_full <- nt_mat_full >= 1L
  dimnames(nt_mat_full) <- list(cell_ids, nt_grna_ids)
  effects <- readRDS(effects_path)
  stopifnot(identical(effects$schema, "frangieh_mean_effect_inputs_v1"), identical(effects$targets, pos_targets),
            identical(effects$preprocessing_contract, contract), identical(effects$Z, covariates),
            identical(colnames(effects$Y), cell_ids), identical(effects$pure_nt, pure_nt))
  pos_X_full <- setNames(lapply(pos_targets, function(g) {
    guides <- positive$grna_id[positive$grna_target == g]
    x <- as.integer(Matrix::colSums(modalities$grna[[guides, ]]) >= 1L)
    stopifnot(identical(unname(effects$X[[g]]), x), identical(unname(masks[[g]]), unname(pure_nt | x == 1L)),
              sum(x) == contract$positive_treated_cells[[g]])
    x
  }), pos_targets)
  null_Y <- as.matrix(modalities$gene[[null_genes, pure_nt_idx]])
  stopifnot(identical(dim(null_Y), c(950L, 3770L)), all(is.finite(null_Y)), all(null_Y >= 0), all(null_Y == round(null_Y)))
  gene_caches <- list()
  for (i in seq_along(null_genes)) {
    g <- null_genes[[i]]
    gene_caches[[paste0("null::", g)]] <- precompute_score_cache(as.numeric(null_Y[i, ]), Z_design_pure_nt)
    if (i %% 100L == 0L) cat("Fitted null-gene caches: ", i, "/950\n", sep = "")
  }
  for (g in pos_targets) gene_caches[[paste0("pos:", g, "::", g)]] <- precompute_score_cache(as.numeric(effects$Y[g, masks[[g]]]), Z_design_full[masks[[g]], , drop = FALSE])
  stopifnot(length(gene_caches) == 1000L)
  canonical_pairs <- canonical_configs <- setNames(vector("list", 20L), names(rep_paths))
  fields <- c("theta_hat", "score_sum", "information_denom_sq", "empirical_denom_sq", "information_stat", "empirical_stat", "score_eff_n", "max_score_share")
  maximum_absolute_difference <- setNames(rep(0, length(fields)), fields)
  input_md5 <- frangieh_file_manifest(file.path(generation, frangieh_intermediate_files))
  used_keys <- character()
  for (r in seq_along(rep_paths)) {
    saved <- readRDS(rep_paths[[r]])
    frangieh_validate_provenance(saved$provenance, workflow, generation)
    stopifnot(identical(saved$config$schema_version, PERTURB_SCHEMA_VERSION), saved$config$rep_idx == r,
              identical(saved$config$preprocessing_contract, contract), identical(saved$provenance$inputs, input_md5),
              identical(saved$config$covariate_names, colnames(covariates)))
    pairs <- saved$result$pairs
    stopifnot(nrow(pairs) == 1000L, identical(pairs$pair_idx, 1:1000), sum(pairs$is_alt) == 50L,
              identical(pairs$gene[pairs$is_alt], pos_targets), all(pairs$gene[!pairs$is_alt] %in% null_genes))
    for (i in seq_len(nrow(pairs))) {
      g <- pairs$gene[[i]]
      if (pairs$is_alt[[i]]) {
        key <- paste0("pos:", g, "::", g)
        x <- pos_X_full[[g]][masks[[g]]]
      } else {
        key <- paste0("null::", g)
        selected <- strsplit(pairs$nt_subset[[i]], ",", fixed = TRUE)[[1L]]
        stopifnot(length(selected) == saved$config$median_size, !anyDuplicated(selected), all(selected %in% nt_grna_ids))
        x <- as.integer(rowSums(nt_mat_full[pure_nt_idx, selected, drop = FALSE]) >= 1L)
      }
      gc <- gene_caches[[key]]
      ps <- precompute_score_pair(x, gc)
      stopifnot(length(x) == pairs$n_cells[[i]], sum(x) == pairs$n_exposed[[i]])
      observed <- c(theta_hat = gc$theta_hat, vapply(ps[fields[-1L]], as.numeric, numeric(1)))
      expected <- unlist(pairs[i, fields], use.names = TRUE)
      delta <- abs(observed - expected)
      if (any(!is.finite(delta)) || any(delta > 1e-7)) stop("Observed score mismatch: replicate ", r, ", pair ", i, ", ", g, "; max absolute difference ", max(delta))
      maximum_absolute_difference <- pmax(maximum_absolute_difference, delta)
      used_keys <- c(used_keys, key)
    }
    canonical_pairs[[r]] <- pairs
    canonical_configs[[r]] <- saved$config
    cat("Validated all 1000 canonical pair scores in replicate ", r, "/20\n", sep = "")
  }
  pair_validation <- list(pairs_checked = 20000L, fields_checked_per_pair = fields,
                          absolute_tolerance = 1e-7, maximum_absolute_difference = as.list(maximum_absolute_difference),
                          cache_fit_methods = as.list(table(vapply(gene_caches, `[[`, character(1), "dispersion_method"))),
                          degenerate_gene_caches = names(gene_caches)[vapply(gene_caches, `[[`, logical(1), "degenerate")],
                          used_degenerate_gene_caches = intersect(unique(used_keys), names(gene_caches)[vapply(gene_caches, `[[`, logical(1), "degenerate")]))
  stopifnot(!length(pair_validation$used_degenerate_gene_caches))
  cache <- list(schema = "frangieh_q10_plain_score_cache_v1", targets = pos_targets,
                cell_ids = cell_ids, Z_design_full = Z_design_full, Z_design_pure_nt = Z_design_pure_nt,
                pure_nt = pure_nt, pure_nt_idx = pure_nt_idx, nt_grna_ids = nt_grna_ids, nt_mat_full = nt_mat_full,
                positive_control_pairs = positive, pos_targets = pos_targets,
                pos_X_full = pos_X_full, pos_cell_masks = masks, canonical_pairs = canonical_pairs,
                canonical_configs = canonical_configs, gene_caches = gene_caches, pair_validation = pair_validation,
                config = config, config_sha256 = config_hash, preprocessing_contract = contract)
  assert_plain <- function(object) {
    if (isS4(object) || is.environment(object) || is.function(object) || typeof(object) == "externalptr") stop("Cache contains an active/non-plain object")
    if (is.list(object)) for (element in object) assert_plain(element)
    invisible(TRUE)
  }
  assert_plain(cache)
  for (entries in list(config$canonical_replicates, config$intermediate_inputs, config$code)) {
    for (entry in entries) stopifnot(identical(sha(entry$path), entry$sha256))
  }
  stopifnot(identical(sha(pointer_path), config$pointer$sha256), identical(sha(effects_path), config$effects_cache$sha256))
  stopifnot(identical(sha(audit_path), config$input_audit$sha256))
  for (entry in raw_records) {
    info <- file.info(entry$path)
    stopifnot(as.numeric(info$size) == entry$bytes,
              identical(format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"), entry$mtime_utc))
  }
  tmp <- tempfile(".q10-cache-", tmpdir = out)
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(cache, tmp, compress = "gzip", version = 3L)
  if (!file.rename(tmp, cache_path)) stop("Could not atomically publish cache")
  manifest <- list(schema = "frangieh_q10_plain_score_cache_manifest_v1", status = "complete",
                   completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                   config = config, config_sha256 = config_hash,
                   config_hash_encoding = "SHA256 of jsonlite::toJSON(config,auto_unbox=TRUE,digits=NA,null='null')",
                   output = record(cache_path), fields = names(cache), pair_validation = pair_validation,
                   dimensions = list(cells = length(cell_ids), pure_nt = length(pure_nt_idx), positives = length(pos_targets),
                                     null_gene_caches = length(null_genes), gene_caches = length(gene_caches), nt_guides = length(nt_grna_ids), replicates = 20L),
                   scope = "Plain numeric/logical matrices and lists; no active ODM objects; exact original cohorts/pairs; no resamples generated; canonical files unmodified")
  temp_json <- tempfile(".q10-cache-manifest-", tmpdir = out)
  on.exit(unlink(temp_json), add = TRUE)
  writeLines(toJSON(manifest, auto_unbox = TRUE, digits = NA, null = "null", pretty = TRUE), temp_json, useBytes = TRUE)
  if (!file.rename(temp_json, manifest_path)) stop("Could not atomically publish manifest")
  cat(toJSON(list(status = "complete", cache = cache_path, bytes = manifest$output$bytes,
                 sha256 = manifest$output$sha256, pair_validation = pair_validation), auto_unbox = TRUE, pretty = TRUE), "\n")
}

main()
