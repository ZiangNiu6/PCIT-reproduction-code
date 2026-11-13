## Filesystem, preprocessing-contract and provenance checks. Scientific
## preprocessing is versioned separately from the five unchanged score methods.

frangieh_intermediate_files <- c(
  "subsetted_ODM.rds", "covariate_matrix.rds", "negative_control_pairs.rds",
  "positive_control_pairs.rds", "gene_list.rds", "null_splits.rds",
  "nt_grna_assigned.rds", "has_any_nt.rds", "pure_nt.rds", "pos_cell_masks.rds",
  "preprocessing-contract.rds"
)

frangieh_validate_preprocessing_contract <- function(contract) {
  expected <- list(
    schema = "frangieh_gb2024_preprocessing_v1",
    feature_stat_policy = "historical_stored_counts",
    cell_rule = "exactly_one_original_assigned_guide_before_feature_qc",
    reference_commit = "005c28b9ba30e7d39c3ee18d51111722f49d32ec",
    ondisc_commit = "5adcc537664aada08c3e79f08f2a944a0d31598d",
    covariates = c("log_n_umis", "log_n_nonzero"),
    id_cleanup = "underscore_to_dash"
  )
  if (!is.list(contract)) stop("Missing GB2024 preprocessing contract")
  for (field in names(expected)) {
    if (!identical(contract[[field]], expected[[field]])) {
      stop("Incompatible GB2024 preprocessing contract field: ", field)
    }
  }
  policy <- if (!is.null(contract[["policy"]])) contract[["policy"]] else contract[["id"]]
  if (!identical(policy, "gb2024_single_guide_historical_stored_counts_v1") ||
      (!is.null(contract[["id"]]) && !identical(contract[["id"]], policy)) ||
      !is.numeric(contract$gene_fraction_min) || length(contract$gene_fraction_min) != 1L ||
      !is.finite(contract$gene_fraction_min) || contract$gene_fraction_min != 0.005 ||
      !is.numeric(contract$guide_cells_min) || length(contract$guide_cells_min) != 1L ||
      !is.finite(contract$guide_cells_min) || contract$guide_cells_min != 10L) {
    stop("Invalid GB2024 preprocessing policy or feature-QC thresholds")
  }
  count_fields <- c("input_cells", "retained_cells", "pure_nt_cells", "input_genes",
                    "retained_genes", "input_guides", "eligible_guides", "eligible_target_groups")
  if (!is.list(contract$counts) || !all(count_fields %in% names(contract$counts)) ||
      !all(vapply(contract$counts[count_fields], function(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
          value >= 0 && value == as.integer(value)
      }, logical(1L))) ||
      contract$counts$retained_cells < 1L ||
      contract$counts$retained_cells > contract$counts$input_cells ||
      contract$counts$pure_nt_cells > contract$counts$retained_cells ||
      contract$counts$retained_genes > contract$counts$input_genes ||
      contract$counts$eligible_guides > contract$counts$input_guides) {
    stop("Invalid GB2024 preprocessing audit counts")
  }
  invisible(contract)
}

frangieh_read_preprocessing_contract <- function(directory) {
  path <- file.path(directory, "preprocessing-contract.rds")
  if (!file.exists(path)) stop("GB2024 preprocessing contract is missing: ", path)
  contract <- readRDS(path)
  frangieh_validate_preprocessing_contract(contract)
  contract
}

frangieh_preprocessing_code_files <- function(script_dir) {
  files <- c("data_preprocessing.R", "workflow-io.R", "method_config.R", "gb2024-preprocessing.R")
  file.path(script_dir, files)
}

frangieh_atomic_rds <- function(object, path, replace = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(path, ".publish-lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    stop("Publication is locked by another process: ", lock_dir,
         ". Retry after that writer finishes; inspect abandoned locks before removing them.")
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  if (file.exists(path) && !replace) stop("Refusing to replace existing file: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not atomically publish: ", path)
  invisible(path)
}

frangieh_intermediate_directory <- function(directory) {
  pointer <- file.path(directory, "current-generation.rds")
  if (!file.exists(pointer)) return(directory)
  generation <- readRDS(pointer)
  if (!is.character(generation) || length(generation) != 1L ||
      !grepl("^generation-[A-Za-z0-9-]+$", generation)) {
    stop("Invalid Frangieh intermediate generation pointer: ", pointer)
  }
  result <- file.path(directory, ".generations", generation)
  if (!dir.exists(result)) stop("Missing intermediate generation: ", result)
  result
}

frangieh_file_manifest <- function(paths, hash = TRUE) {
  if (is.null(names(paths))) names(paths) <- basename(paths)
  if (any(!file.exists(paths))) stop("Missing provenance input: ", paste(paths[!file.exists(paths)], collapse = ", "))
  info <- file.info(paths)
  data.frame(name = names(paths), size = as.numeric(info$size),
             md5 = if (hash) unname(tools::md5sum(paths)) else NA_character_,
             stringsAsFactors = FALSE, row.names = NULL)
}

frangieh_code_files <- function(script_dir) {
  files <- c("run_combined_fdr.R", "method_config.R", "score_helpers.R",
             "workflow-io.R", "result-validation.R", "gb2024-preprocessing.R")
  file.path(script_dir, files)
}

frangieh_provenance <- function(intermediate_dir, script_dir, kind = "native_run") {
  if (!kind %in% c("native_run", "imported_verified_q10")) stop("Invalid provenance kind")
  packages <- c("pcmt", "ondisc", "MASS", "dplyr", "Matrix")
  list(schema = PERTURB_PROVENANCE_SCHEMA, kind = kind,
       role = if (kind == "native_run") "native_execution" else "current_compatible_implementation_and_inputs",
       inputs = frangieh_file_manifest(file.path(intermediate_dir, frangieh_intermediate_files)),
       code = frangieh_file_manifest(frangieh_code_files(script_dir)),
       packages = setNames(vapply(packages, function(package) as.character(utils::packageVersion(package)), character(1)), packages),
       R = as.character(getRversion()), session = utils::sessionInfo(),
       cohort_policy = PERTURB_COHORT_POLICY,
       preprocessing_contract = frangieh_read_preprocessing_contract(intermediate_dir))
}

frangieh_archived_rds <- function(bytes) {
  if (!is.raw(bytes) || !length(bytes)) stop("Archived selected execution must contain raw RDS bytes")
  connection <- rawConnection(bytes, "rb")
  if (length(bytes) >= 2L && identical(bytes[1:2], as.raw(c(31, 139)))) connection <- gzcon(connection)
  on.exit(close(connection), add = TRUE)
  readRDS(connection)
}

frangieh_validate_adoption <- function(adoption, result = NULL) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to verify selected-execution provenance")
  if (!is.list(adoption) || !identical(adoption$schema, "frangieh_selected_q10_adoption_v1") ||
      !is.list(adoption$archive) || !is.raw(adoption$archive$original_rds)) stop("Invalid selected-execution adoption record")
  bytes <- adoption$archive$original_rds
  record <- adoption$source_rds
  receipt <- adoption$original_receipt
  sha <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
  if (!is.list(record) || !is.numeric(record$bytes) || record$bytes != length(bytes) ||
      !identical(record$sha256, sha) || !is.list(receipt) || !identical(receipt$status, "complete") ||
      length(receipt$outputs) != 1L || receipt$outputs[[1L]]$bytes != length(bytes) ||
      !identical(receipt$outputs[[1L]]$sha256, sha)) stop("Archived selected execution does not match its original receipt")
  original <- frangieh_archived_rds(bytes)
  cfg <- original$config
  if (!identical(cfg$schema_version, "frangieh_score5_q10_extension_v1") || isTRUE(cfg$pilot) ||
      !identical(cfg$q, PERTURB_HEADLINE_Q) || cfg$h != 15L || cfg$B_max != 25000L ||
      !identical(cfg$method_levels, PERTURB_METHOD_LEVELS) ||
      !identical(cfg$RNGkind, PERTURB_RNG_KIND) || cfg$seed_reference_B_max != 25000L ||
      length(cfg$seed_reference_q_levels) != 2L ||
      match(cfg$seed_reference_q, cfg$seed_reference_q_levels) != 2L ||
      !identical(receipt$fingerprint, original$provenance$fingerprint) ||
      !identical(receipt$fingerprint, digest::digest(list(config = cfg, inputs = original$provenance$inputs,
        sources = original$provenance$sources), algo = "sha256"))) stop("Invalid original selected q10 execution metadata")
  for (field in c("config", "inputs", "sources")) {
    expected <- if (field == "config") cfg else original$provenance[[field]]
    if (!identical(jsonlite::toJSON(receipt[[field]], auto_unbox = TRUE, digits = NA, null = "null"),
                   jsonlite::toJSON(expected, auto_unbox = TRUE, digits = NA, null = "null"))) {
      stop("Original receipt and selected execution metadata disagree: ", field)
    }
  }
  snapshots <- adoption$archive$source_snapshots
  if (!is.list(snapshots) || is.null(names(snapshots)) || anyDuplicated(names(snapshots))) stop("Missing immutable execution-source snapshots")
  for (source in original$provenance$sources) {
    bytes <- snapshots[[basename(source$path)]]
    if (!is.raw(bytes) || length(bytes) != source$bytes ||
        !identical(digest::digest(bytes, algo = "sha256", serialize = FALSE), source$sha256)) {
      stop("Archived execution source does not match its original receipt: ", basename(source$path))
    }
  }
  expected_result <- digest::digest(original$result, algo = "sha256")
  if (!identical(adoption$result_sha256, expected_result) ||
      (!is.null(result) && !identical(result, original$result))) stop("Adopted numerical result differs from the selected execution")
  invisible(original)
}

frangieh_adopt_verified_q10 <- function(original_object, source_rds_record, original_receipt,
                                       execution_archive_manifest, current_provenance, current_config) {
  adoption <- list(schema = "frangieh_selected_q10_adoption_v1", source_rds = source_rds_record,
                   original_receipt = original_receipt, archive = execution_archive_manifest,
                   result_sha256 = digest::digest(original_object$result, algo = "sha256"))
  archived <- frangieh_validate_adoption(adoption, original_object$result)
  if (!identical(archived, original_object)) stop("Original execution object and archived RDS differ")
  current_provenance$kind <- "imported_verified_q10"
  current_provenance$role <- "current_compatible_implementation_and_inputs"
  current_provenance$adoption <- adoption
  list(config = current_config, result = original_object$result, provenance = current_provenance)
}

frangieh_validate_provenance <- function(provenance, script_dir, intermediate_dir = NULL, result = NULL) {
  if (!is.list(provenance) || !identical(provenance$schema, PERTURB_PROVENANCE_SCHEMA)) {
    stop("Invalid Frangieh provenance schema")
  }
  if (!is.character(provenance$kind) || length(provenance$kind) != 1L ||
      !provenance$kind %in% c("native_run", "imported_verified_q10") ||
      !identical(provenance$role, if (provenance$kind == "native_run") "native_execution" else "current_compatible_implementation_and_inputs")) stop("Invalid provenance kind or role")
  for (field in c("inputs", "code")) {
    manifest <- provenance[[field]]
    if (!is.data.frame(manifest) || !all(c("name", "size", "md5") %in% names(manifest)) ||
        anyNA(manifest) || anyDuplicated(manifest$name) ||
        any(!grepl("^[a-f0-9]{32}$", manifest$md5))) stop("Invalid provenance ", field, " manifest")
  }
  if (!identical(provenance$inputs$name, frangieh_intermediate_files)) stop("Invalid provenance input coverage")
  # File/code fingerprints describe the execution. They are metadata, not a
  # requirement that today's source tree or raw inputs match a saved result.
  if (!is.character(provenance$packages) || anyNA(provenance$packages) ||
      !identical(names(provenance$packages), c("pcmt", "ondisc", "MASS", "dplyr", "Matrix")) ||
      any(!nzchar(provenance$packages)) || !is.character(provenance$R) ||
      length(provenance$R) != 1L || is.na(provenance$R) || !nzchar(provenance$R) ||
      !identical(provenance$cohort_policy, PERTURB_COHORT_POLICY)) {
    stop("Invalid package, R, or cohort provenance")
  }
  frangieh_validate_preprocessing_contract(provenance$preprocessing_contract)
  if (!is.null(intermediate_dir) && !identical(
    provenance$preprocessing_contract, frangieh_read_preprocessing_contract(intermediate_dir)
  )) stop("Recorded preprocessing contract differs from the intermediate bundle")
  if (provenance$kind == "imported_verified_q10") frangieh_validate_adoption(provenance$adoption, result)
  else if (!is.null(provenance$adoption)) stop("Native execution must not claim imported execution metadata")
  invisible(TRUE)
}

frangieh_validate_intermediates <- function(directory, script_dir = NULL) {
  directory <- frangieh_intermediate_directory(directory)
  missing <- frangieh_intermediate_files[!file.exists(file.path(directory, frangieh_intermediate_files))]
  if (length(missing)) stop("Incomplete Frangieh intermediate bundle: ", paste(missing, collapse = ", "))
  objects <- setNames(lapply(file.path(directory, frangieh_intermediate_files), readRDS), frangieh_intermediate_files)
  contract <- frangieh_read_preprocessing_contract(directory)
  z <- objects[["covariate_matrix.rds"]]
  genes <- objects[["gene_list.rds"]]
  positive <- objects[["positive_control_pairs.rds"]]
  negative <- objects[["negative_control_pairs.rds"]]
  pure_nt <- objects[["pure_nt.rds"]]
  masks <- objects[["pos_cell_masks.rds"]]
  n <- nrow(z)
  if (!is.matrix(z) || !is.numeric(z) || any(!is.finite(z)) ||
      !is.character(genes) || length(genes) != 1000L || anyDuplicated(genes) ||
      !is.data.frame(positive) || !all(c("grna_id", "grna_target", "grna_group") %in% names(positive)) ||
      length(unique(positive$grna_target)) != 50L ||
      !is.logical(pure_nt) || length(pure_nt) != n || anyNA(pure_nt) ||
      !is.list(masks) || !setequal(names(masks), unique(positive$grna_target)) ||
      any(!vapply(masks, function(x) is.logical(x) && length(x) == n && !anyNA(x) && all(x[pure_nt]), logical(1)))) {
    stop("Invalid Frangieh intermediate dimensions or cohort metadata")
  }
  if (!is.list(negative) || length(negative$gene_id) != 950L ||
      anyDuplicated(negative$gene_id) || !all(negative$gene_id %in% genes) ||
      any(negative$gene_id %in% positive$grna_target) ||
      !is.data.frame(negative$grna) || !all(c("grna_id", "grna_group") %in% names(negative$grna))) {
    stop("Invalid Frangieh negative-control intermediate metadata")
  }
  modalities <- objects[["subsetted_ODM.rds"]]
  if (!is.list(modalities) || !all(c("gene", "grna") %in% names(modalities)) ||
      ncol(modalities$gene) != n || ncol(modalities$grna) != n) stop("Invalid saved ondisc modality dimensions")
  if (nrow(objects[["nt_grna_assigned.rds"]]) != n ||
      length(objects[["has_any_nt.rds"]]) != n) stop("Invalid saved NT assignment dimensions")
  if (!identical(colnames(z), contract$covariates) ||
        !identical(unname(pure_nt), unname(objects[["has_any_nt.rds"]])) ||
        n != contract$counts$retained_cells || sum(pure_nt) != contract$counts$pure_nt_cells) {
    stop("GB2024 covariates or single-guide NT controls disagree with the contract")
  }
  manifest_path <- file.path(directory, "intermediate-provenance.rds")
  if (!file.exists(manifest_path)) stop("GB2024 intermediates require intermediate-provenance.rds")
  manifest <- readRDS(manifest_path)
  if (!is.list(manifest) || !identical(manifest$schema, PERTURB_INTERMEDIATE_SCHEMA) ||
      !identical(manifest$cohort_policy, PERTURB_COHORT_POLICY)) {
    stop("Invalid Frangieh intermediate provenance schema or cohort policy")
  }
  # These metadata identify the original QC/count dataset used by the saved
  # ODM subsets, independently of the code that produced the intermediates.
  raw_root <- Sys.getenv("LOCAL_FRANGIEH_2021_DATA_DIR")
  if (nzchar(raw_root)) {
    raw_metadata <- file.path(raw_root, "processed", "control",
      c("gene/gene_expression_metadata.rds", "grna_assignment/grna_assignments_ungrouped_metadata.rds"))
    if (!identical(manifest$raw_metadata, frangieh_file_manifest(raw_metadata))) {
      stop("Frangieh raw metadata do not match the dataset used for these intermediates")
    }
  }
  invisible(directory)
}

## Rebind only the backing-store path. Check the original full-file identity and
## preserve every saved subset, covariate, assignment, and cell/feature order.
frangieh_rebind_modalities <- function(modalities, raw_root = Sys.getenv("LOCAL_FRANGIEH_2021_DATA_DIR")) {
  for (modality in c("gene", "grna")) {
    saved <- modalities[[modality]]
    if (!nzchar(raw_root)) {
      if (!file.exists(saved@ondisc_matrix@h5_file)) stop("Saved ODM backing path is unavailable. Set LOCAL_FRANGIEH_2021_DATA_DIR to relocate it.")
      next
    }
    relative <- if (modality == "gene") c("gene", "gene_expression_matrix.odm", "gene_expression_metadata.rds") else
      c("grna_assignment", "grna_assignments_ungrouped.odm", "grna_assignments_ungrouped_metadata.rds")
    base <- file.path(raw_root, "processed", "control", relative[[1L]])
    fresh <- ondisc::read_odm(odm_fp = file.path(base, relative[[2L]]), metadata_fp = file.path(base, relative[[3L]]))
    for (field in c("underlying_dimension", "odm_id", "logical_mat")) {
      if (!identical(slot(saved@ondisc_matrix, field), slot(fresh@ondisc_matrix, field))) {
        stop("Relocated ", modality, " ODM does not match saved input identity (", field, ")")
      }
    }
    underlying_indices <- function(matrix, subset_field, dimension, id_field) {
      subset <- slot(matrix, subset_field)
      indices <- if (length(subset) == 1L && is.na(subset)) {
        seq_len(matrix@underlying_dimension[[dimension]])
      } else subset
      if (!is.numeric(indices) || anyNA(indices) || anyDuplicated(indices) ||
          any(indices != as.integer(indices)) || any(indices < 1L) ||
          any(indices > matrix@underlying_dimension[[dimension]]) ||
          length(indices) != length(slot(matrix, id_field))) {
        stop("Invalid saved/raw ODM subset mapping for ", id_field)
      }
      as.integer(indices)
    }
    for (axis in list(
      list(subset = "feature_subset", dimension = 1L, ids = "feature_ids"),
      list(subset = "cell_subset", dimension = 2L, ids = "cell_barcodes")
    )) {
      saved_indices <- underlying_indices(saved@ondisc_matrix, axis$subset, axis$dimension, axis$ids)
      fresh_indices <- underlying_indices(fresh@ondisc_matrix, axis$subset, axis$dimension, axis$ids)
      matched <- match(saved_indices, fresh_indices)
      if (anyNA(matched)) stop("Relocated ", modality, " ODM omits saved ", axis$ids)
      expected <- slot(fresh@ondisc_matrix, axis$ids)[matched]
      if (identical(axis$ids, "feature_ids")) {
        expected <- gsub("_", "-", expected, fixed = TRUE)
      }
      if (!identical(slot(saved@ondisc_matrix, axis$ids), expected)) {
        stop("Relocated ", modality, " ODM does not match saved input identity (", axis$ids, ")")
      }
    }
    modalities[[modality]]@ondisc_matrix@h5_file <- normalizePath(file.path(base, relative[[2L]]))
  }
  modalities
}
