## Reusable structural validation; no analysis is executed or output written.
schema_error <- function(path, ...) {
  stop(basename(path), ": incompatible ", PERTURB_SCHEMA_VERSION,
       " result: ", paste0(..., collapse = ""), call. = FALSE)
}

check_rejection_indices <- function(x, n_pairs, path, label) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != as.integer(x)) || any(x < 1L | x > n_pairs) ||
      anyDuplicated(x)) {
    schema_error(path, label, " must contain unique pair indices in 1:",
                 n_pairs)
  }
  invisible(TRUE)
}

check_p_values <- function(x, n_pairs, path, label) {
  if (!is.numeric(x) || length(x) != n_pairs || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0 | x > 1)) {
    schema_error(path, label, " must be a finite numeric vector of length ",
                 n_pairs, " with values in [0, 1]")
  }
  invisible(TRUE)
}

validate_rep <- function(obj, path, expected_rep = NULL, intermediate_dir = NULL) {
  if (!is.list(obj) || !all(c("config", "result") %in% names(obj))) {
    schema_error(path, "top level must contain config and result")
  }
  cfg <- obj$config
  res <- obj$result
  if (!identical(cfg$schema_version, PERTURB_SCHEMA_VERSION)) {
    schema_error(
      path, "schema_version is ", deparse(cfg$schema_version), "; expected ",
      deparse(PERTURB_SCHEMA_VERSION),
      ". Do not mix legacy and corrected result files"
    )
  }
  if (!identical(cfg$method_levels, PERTURB_METHOD_LEVELS)) {
    schema_error(path, "config$method_levels does not match method_config.R")
  }
  required_config <- c(
    "n_reps_max", "rep_idx", "seed", "n_pos", "n_null", "n_total",
    "median_size", "n_cells_full", "n_pure_nt", "q_levels", "h",
    "B_max", "alternative", "null_gene_policy", "null_gene_pool_size",
    "null_gene_pool"
  )
  if (!all(required_config %in% names(cfg))) {
    schema_error(path, "config is missing: ",
                 paste(setdiff(required_config, names(cfg)), collapse = ", "))
  }
  if (identical(PERTURB_ANALYSIS_SCHEMA, "gb2024")) {
    required_gb <- c("preprocessing_policy", "preprocessing_contract", "covariate_names", "n_nt_exposure_guides")
    if (!all(required_gb %in% names(cfg)) ||
        !identical(cfg$preprocessing_policy, PERTURB_COHORT_POLICY)) {
      schema_error(path, "missing or incompatible GB2024 preprocessing policy")
    }
    frangieh_validate_preprocessing_contract(cfg$preprocessing_contract)
    contract <- cfg$preprocessing_contract
    if (!identical(cfg$covariate_names, contract$covariates) ||
        !is.numeric(cfg$n_cells_full) || length(cfg$n_cells_full) != 1L ||
        !is.finite(cfg$n_cells_full) || cfg$n_cells_full != contract$counts$retained_cells ||
        !is.numeric(cfg$n_pure_nt) || length(cfg$n_pure_nt) != 1L ||
        !is.finite(cfg$n_pure_nt) || cfg$n_pure_nt != contract$counts$pure_nt_cells ||
        !is.numeric(cfg$n_nt_exposure_guides) || length(cfg$n_nt_exposure_guides) != 1L ||
        !is.finite(cfg$n_nt_exposure_guides) || cfg$n_nt_exposure_guides != as.integer(cfg$n_nt_exposure_guides) ||
        cfg$n_nt_exposure_guides < 1L || cfg$n_nt_exposure_guides > contract$counts$eligible_guides) {
      schema_error(path, "GB2024 cohort sizes or covariates disagree with the preprocessing contract")
    }
  }
  if (length(cfg$h) != 1L || !is.numeric(cfg$h) || cfg$h < 1L ||
      length(cfg$B_max) != 1L || !is.numeric(cfg$B_max) || cfg$B_max < 1L ||
      !identical(cfg$alternative, "two.sided")) {
    schema_error(path, "config has invalid h, B_max, or alternative")
  }
  if (!identical(cfg$null_gene_policy, PERTURB_NULL_GENE_POLICY) ||
      !is.character(cfg$null_gene_pool) || anyNA(cfg$null_gene_pool) ||
      any(cfg$null_gene_pool == "") || anyDuplicated(cfg$null_gene_pool) ||
      !identical(
        as.integer(cfg$null_gene_pool_size),
        as.integer(length(cfg$null_gene_pool))
      )) {
    schema_error(path, "config has an invalid non-positive null-gene pool")
  }

  q_levels <- as.numeric(cfg$q_levels)
  if (length(q_levels) == 0L || anyNA(q_levels) ||
      any(q_levels <= 0 | q_levels >= 1) || anyDuplicated(q_levels)) {
    schema_error(path, "config$q_levels must be unique probabilities in (0, 1)")
  }
  q_keys <- perturb_q_key(q_levels)
  canonical <- list(n_reps_max = 20L, n_pos = 50L, n_null = 950L,
                    n_total = 1000L, h = 15L, B_max = 25000L)
  if (identical(PERTURB_ANALYSIS_SCHEMA, "legacy-v2")) canonical$median_size <- 3L
  for (field in names(canonical)) {
    value <- cfg[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != canonical[[field]]) schema_error(path, "noncanonical config$", field)
  }
  if (identical(PERTURB_ANALYSIS_SCHEMA, "gb2024") &&
      (!is.numeric(cfg$median_size) || length(cfg$median_size) != 1L ||
       !is.finite(cfg$median_size) || cfg$median_size != as.integer(cfg$median_size) ||
       cfg$median_size < 1L || cfg$median_size > cfg$n_nt_exposure_guides)) {
    schema_error(path, "invalid target-matched NT subset size")
  }
  if (!identical(q_levels, c(0.01, 0.05))) schema_error(path, "noncanonical q grid")
  if (!is.numeric(cfg$rep_idx) || length(cfg$rep_idx) != 1L ||
      !is.finite(cfg$rep_idx) || cfg$rep_idx != as.integer(cfg$rep_idx) ||
      !cfg$rep_idx %in% 1:20 || !identical(as.integer(cfg$seed), 1000L + as.integer(cfg$rep_idx))) {
    schema_error(path, "invalid replicate index or deterministic seed")
  }
  filename_rep <- suppressWarnings(as.integer(sub("^score_fdr_rep([0-9]{3})[.]rds$", "\\1", basename(path))))
  if (is.na(filename_rep) || filename_rep != cfg$rep_idx ||
      (!is.null(expected_rep) && cfg$rep_idx != expected_rep)) schema_error(path, "replicate does not match its expected filename/index")
  frangieh_validate_provenance(obj$provenance, script_dir, intermediate_dir)
  if (identical(PERTURB_ANALYSIS_SCHEMA, "gb2024") &&
      !identical(obj$provenance$preprocessing_contract, cfg$preprocessing_contract)) {
    schema_error(path, "result configuration and provenance preprocessing contracts differ")
  }

  required_result <- c("rep", "seed", "pairs", "methods", "timings")
  if (!all(required_result %in% names(res))) {
    schema_error(path, "result is missing: ",
                 paste(setdiff(required_result, names(res)), collapse = ", "))
  }
  if (length(res$rep) != 1L || !is.numeric(res$rep) || is.na(res$rep)) {
    schema_error(path, "result$rep must be one non-missing number")
  }
  if (!identical(as.integer(res$rep), as.integer(cfg$rep_idx))) {
    schema_error(path, "result$rep and config$rep_idx disagree")
  }
  if (length(res$seed) != 1L || !is.numeric(res$seed) || is.na(res$seed) ||
      !identical(as.integer(res$seed), as.integer(cfg$seed))) {
    schema_error(path, "result$seed and config$seed must agree")
  }

  required_pair_columns <- c(
    "pair_idx", "gene", "pair_type", "is_alt", "nt_subset", "n_cells",
    "n_exposed", "theta_hat", "score_sum", "information_denom_sq",
    "empirical_denom_sq", "information_stat", "empirical_stat",
    "score_eff_n", "max_score_share"
  )
  pairs <- res$pairs
  if (!is.data.frame(pairs) ||
      !all(required_pair_columns %in% names(pairs))) {
    schema_error(path, "result$pairs is missing: ",
                 paste(setdiff(required_pair_columns, names(pairs)),
                       collapse = ", "))
  }
  n_pairs <- nrow(pairs)
  if (n_pairs == 0L || !identical(as.integer(pairs$pair_idx), seq_len(n_pairs))) {
    schema_error(path, "pairs$pair_idx must equal 1:nrow(pairs)")
  }
  if (!is.logical(pairs$is_alt) || anyNA(pairs$is_alt)) {
    schema_error(path, "pairs$is_alt must be a non-missing logical vector")
  }
  if (anyNA(pairs$pair_type) ||
      !all(as.character(pairs$pair_type) == ifelse(pairs$is_alt, "alt", "null"))) {
    schema_error(path, "pairs$pair_type and pairs$is_alt disagree")
  }
  if (!is.null(cfg$n_total) && as.integer(cfg$n_total) != n_pairs) {
    schema_error(path, "config$n_total and nrow(result$pairs) disagree")
  }
  if (sum(pairs$is_alt) != as.integer(cfg$n_pos) ||
      sum(!pairs$is_alt) != as.integer(cfg$n_null)) {
    schema_error(path, "configured alternative/null counts disagree with pairs")
  }
  if (identical(PERTURB_ANALYSIS_SCHEMA, "gb2024")) {
    if (anyNA(pairs$n_cells) || any(!is.finite(pairs$n_cells)) ||
        any(pairs$n_cells != as.integer(pairs$n_cells)) ||
        any(pairs$n_cells < 1L | pairs$n_cells > cfg$n_cells_full) ||
        anyNA(pairs$n_exposed) || any(!is.finite(pairs$n_exposed)) ||
        any(pairs$n_exposed != as.integer(pairs$n_exposed)) ||
        any(pairs$n_exposed < 1L | pairs$n_exposed > pairs$n_cells) ||
        any(pairs$n_cells[!pairs$is_alt] != cfg$n_pure_nt) ||
        any(pairs$n_cells[pairs$is_alt] - cfg$n_pure_nt != pairs$n_exposed[pairs$is_alt])) {
      schema_error(path, "pair populations disagree with the GB2024 single-guide control cohort")
    }
  }
  alt_genes <- unique(as.character(pairs$gene[pairs$is_alt]))
  null_genes <- as.character(pairs$gene[!pairs$is_alt])
  if (anyNA(null_genes) || any(null_genes == "") ||
      any(null_genes %in% alt_genes) ||
      !all(null_genes %in% cfg$null_gene_pool)) {
    schema_error(
      path,
      "null pairs must use only configured non-positive response genes"
    )
  }

  methods <- res$methods
  if (!is.list(methods) || !identical(names(methods), PERTURB_METHOD_LEVELS)) {
    schema_error(path, "result$methods must have the configured methods in order")
  }
  for (method in PERTURB_METHOD_LEVELS) {
    method_result <- methods[[method]]
    if (!is.list(method_result) || !is.numeric(method_result$observed) ||
        length(method_result$observed) != n_pairs) {
      schema_error(path, "methods$", method,
                   "$observed must be numeric with one value per pair")
    }
    if (method == "parametric_score") {
      check_p_values(method_result$p_value, n_pairs, path,
                     "methods$parametric_score$p_value")
      if (!is.list(method_result$rejected) ||
          !setequal(names(method_result$rejected), q_keys)) {
        schema_error(path, "methods$parametric_score$rejected must be keyed by q")
      }
      for (key in q_keys) {
        check_rejection_indices(
          method_result$rejected[[key]], n_pairs, path,
          paste0("methods$parametric_score$rejected[[", key, "]]")
        )
      }
    } else {
      if (!is.list(method_result$by_q) ||
          !setequal(names(method_result$by_q), q_keys)) {
        schema_error(path, "methods$", method, "$by_q must be keyed by q")
      }
      for (key in q_keys) {
        by_q <- method_result$by_q[[key]]
        required_by_q <- c(
          "p_value", "rejected", "futility", "active", "K",
          "stop_round", "t_total", "R_threshold"
        )
        if (!is.list(by_q) || !all(required_by_q %in% names(by_q))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]] is missing adaptive-result fields")
        }
        check_p_values(by_q$p_value, n_pairs, path,
                       paste0("methods$", method, "$by_q[[", key,
                              "]]$p_value"))
        check_rejection_indices(
          by_q$rejected, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$rejected")
        )
        check_rejection_indices(
          by_q$futility, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$futility")
        )
        check_rejection_indices(
          by_q$active, n_pairs, path,
          paste0("methods$", method, "$by_q[[", key, "]]$active")
        )
        partition <- c(by_q$rejected, by_q$futility, by_q$active)
        if (!identical(sort(as.integer(partition)), seq_len(n_pairs))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]] rejection/futility/active sets do not partition pairs")
        }
        if (!is.numeric(by_q$K) || length(by_q$K) != n_pairs ||
            anyNA(by_q$K) || any(by_q$K != as.integer(by_q$K)) ||
            any(by_q$K < 0L | by_q$K > as.integer(cfg$h))) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$K must contain one loss count in 0:h per pair")
        }
        if (!is.numeric(by_q$stop_round) || length(by_q$stop_round) != n_pairs) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$stop_round must have one value per pair")
        }
        if (!all(is.na(by_q$stop_round[by_q$active])) ||
            anyNA(by_q$stop_round[c(by_q$rejected, by_q$futility)])) {
          schema_error(path, "methods$", method, "$by_q[[", key,
                       "]]$stop_round must be NA exactly for active pairs")
        }
      }
    }
  }

  timings <- res$timings
  if (!is.list(timings) || !is.numeric(timings$cache) ||
      length(timings$cache) != 1L || !is.finite(timings$cache)) {
    schema_error(path, "result$timings$cache must be one finite number")
  }
  for (method in PERTURB_RESAMPLING_METHODS) {
    method_timings <- timings[[method]]
    if (!is.numeric(method_timings) ||
        !setequal(names(method_timings), q_keys) ||
        anyNA(method_timings) || any(!is.finite(method_timings))) {
      schema_error(path, "result$timings$", method,
                   " must be a finite numeric vector keyed by q")
    }
  }

  obj
}


frangieh_result_files <- function(results_dir) {
  files <- sort(list.files(results_dir, pattern = "^score_fdr_rep", full.names = TRUE))
  expected <- sprintf("score_fdr_rep%03d.rds", 1:20)
  unexpected <- files[!basename(files) %in% expected]
  if (length(unexpected)) stop("Unexpected replicate filenames: ", paste(basename(unexpected), collapse = ", "))
  files
}

frangieh_validate_collection <- function(files, require_complete = TRUE, intermediate_dir = NULL) {
  objects <- lapply(files, function(path) validate_rep(readRDS(path), path, intermediate_dir = intermediate_dir))
  ids <- unname(vapply(objects, function(x) as.integer(x$result$rep), integer(1)))
  if (anyDuplicated(ids)) stop("Duplicate replicate IDs in score result files")
  if (require_complete && !identical(sort(ids), 1:20)) stop("Canonical outputs require exactly replicate IDs 1:20; found: ", paste(sort(ids), collapse = ", "))
  if (length(objects) > 1L) {
    reference <- objects[[1L]]$config[frangieh_analysis_config_fields]
    for (i in seq_along(objects)[-1L]) {
      if (!identical(objects[[i]]$config[frangieh_analysis_config_fields], reference)) stop(basename(files[[i]]), ": analysis configuration differs between replicates")
    }
    recorded <- which(vapply(objects, function(x) !is.null(x$provenance), logical(1)))
    if (length(recorded) > 1L) {
      signatures <- lapply(objects[recorded], function(x) x$provenance[c("inputs", "code", "packages", "R", "cohort_policy", "preprocessing_contract")])
      if (any(!vapply(signatures, identical, logical(1), signatures[[1L]]))) stop("Inconsistent provenance across new replicate outputs")
    }
  }
  objects
}

frangieh_analysis_config_fields <- c(
  "schema_version", "method_levels", "n_reps_max", "n_pos", "n_null",
  "n_total", "median_size", "n_cells_full", "n_pure_nt", "q_levels",
  "h", "B_max", "alternative", "null_gene_policy", "null_gene_pool_size", "null_gene_pool"
)
if (identical(PERTURB_ANALYSIS_SCHEMA, "gb2024")) {
  frangieh_analysis_config_fields <- c(
    frangieh_analysis_config_fields, "preprocessing_policy", "preprocessing_contract",
    "covariate_names", "n_nt_exposure_guides"
  )
}
