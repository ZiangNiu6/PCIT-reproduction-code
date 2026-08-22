#!/usr/bin/env Rscript
## Strict combiner for the standalone Frangieh exploratory trans-pair shards.
## It verifies a complete declared pilot or full-grid universe before writing
## anything.  Outputs contain raw permutation p-values only: no BH/q-values
## are computed because a fixed, modest B has inadequate extreme-tail
## resolution for multiplicity claims over 5,000 or 1.2M hypotheses.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_path <- normalizePath(sub("^--file=", "", script_arg))

usage <- function() {
  paste(
    "Usage: Rscript --vanilla combine_frangieh_exploratory_trans_pairs.R",
    "  --input-dir DIR --output-dir DIR --mode pilot|full",
    "  --min-treated 10|5 --permutations B --seed SEED [--overwrite]",
    sep = "\n"
  )
}

parse_cli <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    if (token %in% c("--help", "--overwrite")) {
      out[[gsub("-", "_", sub("^--", "", token), fixed = TRUE)]] <- TRUE
      i <- i + 1L
      next
    }
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
      out[[gsub("-", "_", bits[[1L]], fixed = TRUE)]] <-
        paste(bits[-1L], collapse = "=")
      i <- i + 1L
      next
    }
    if (i == length(x)) stop("Missing value after ", token)
    out[[gsub("-", "_", sub("^--", "", token), fixed = TRUE)]] <-
      x[[i + 1L]]
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

same_named_atomic <- function(x, y) {
  identical(names(x), names(y)) && identical(unname(x), unname(y))
}

atomic_publish <- function(tmp, destination, overwrite = FALSE) {
  if (file.exists(destination) && !overwrite) {
    stop("Refusing to replace existing output: ", destination)
  }
  if (file.exists(destination) && overwrite && !file.remove(destination)) {
    stop("Could not replace: ", destination)
  }
  if (!file.rename(tmp, destination)) {
    stop("Could not atomically publish: ", destination)
  }
  invisible(destination)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
allowed_cli <- c("help", "input_dir", "output_dir", "mode", "min_treated",
                 "permutations", "seed", "overwrite")
unknown_cli <- setdiff(names(cli), allowed_cli)
if (length(unknown_cli)) {
  stop("Unknown option(s): --", paste(gsub("_", "-", unknown_cli),
                                      collapse = ", --"))
}
if (isTRUE(cli[["help"]])) {
  cat(usage(), "\n")
  quit(save = "no", status = 0L)
}
required_cli <- c("input_dir", "output_dir", "mode", "min_treated",
                  "permutations", "seed")
missing_cli <- required_cli[!vapply(
  required_cli, function(key) !is.null(cli[[key]]), logical(1)
)]
if (length(missing_cli)) {
  stop("Missing required option(s): --",
       paste(gsub("_", "-", missing_cli), collapse = ", --"), "\n", usage())
}

input_dir <- normalizePath(cli[["input_dir"]])
output_dir <- cli[["output_dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
mode <- match.arg(cli[["mode"]], c("pilot", "full"))
min_treated <- as_count(cli[["min_treated"]], "--min-treated", minimum = 5L)
if (!min_treated %in% c(5L, 10L)) {
  stop("--min-treated must be 10 (primary) or 5 (sensitivity)")
}
if (mode == "pilot" && min_treated != 10L) {
  stop("The frozen pilot is defined only for --min-treated 10")
}
B <- as_count(cli[["permutations"]], "--permutations", minimum = 19L)
seed <- as_count(cli[["seed"]], "--seed")
overwrite <- isTRUE(cli[["overwrite"]])

mode_pattern <- if (mode == "pilot") {
  "-pilot5000-"
} else {
  "-r[0-9]{4}-[0-9]{4}-"
}
all_rds <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
candidate <- all_rds[
  startsWith(basename(all_rds), "frangieh-exploratory-trans-e") &
    grepl(mode_pattern, basename(all_rds)) &
    grepl(paste0("-min", min_treated, "-B", B, "\\.rds$"),
          basename(all_rds))
]
candidate <- sort(candidate)
if (!length(candidate)) {
  stop("No matching shard RDS files found in ", input_dir)
}
cat("Strict Frangieh exploratory combiner\n")
cat("  mode/min/B/seed: ", mode, " / ", min_treated, " / ", B,
    " / ", seed, "\n", sep = "")
cat("  candidate files: ", length(candidate), "\n", sep = "")

expected_analysis <- "frangieh_label_invariant_exploratory_trans_pairs"
expected_scan_mode <- if (mode == "pilot") {
  "association_blind_pilot_5000"
} else {
  "declared_grid_or_range"
}
## The runner derives and publishes the frozen-pilot manifest whenever
## min_treated == 10, whether or not the shard itself is a --pilot shard, so
## the expected manifest set keys on min_treated rather than on the mode.
## Keying it on the mode made every full-grid min-treated-10 shard fail the
## provenance gate below.
expected_manifest_names <- c(
  "expression_audit", "response_universe", "selected_responses", "exposures"
)
if (min_treated == 10L) {
  expected_manifest_names <- c(
    expected_manifest_names, "association_blind_pilot_5000"
  )
}

baseline <- NULL
pair_chunks <- vector("list", length(candidate))
excluded_chunks <- vector("list", length(candidate))
shard_rows <- vector("list", length(candidate))

for (i in seq_along(candidate)) {
  file <- candidate[[i]]
  shard <- readRDS(file)
  needed <- c("config", "data_checks", "provenance", "exposure_manifest",
              "response_manifest", "pair_manifest", "pair_results")
  if (!all(needed %in% names(shard))) {
    stop("Shard lacks required objects: ", basename(file))
  }
  cfg <- shard$config
  checks <- shard$data_checks
  prov <- shard$provenance
  pairs <- shard$pair_results
  declared <- shard$pair_manifest
  required_provenance <- c(
    "counts_backend", "counts_source_file", "counts_source_size_bytes",
    "counts_source_mtime_utc", "counts_source_md5",
    "counts_identity_verification", "input_md5", "script_md5",
    "score_helper_md5", "package_versions"
  )
  if (!all(required_provenance %in% names(prov)) ||
      !identical(names(checks$manifest_md5), expected_manifest_names) ||
      anyNA(checks$manifest_md5) ||
      any(!grepl("^[0-9a-f]{32}$", checks$manifest_md5)) ||
      length(prov$counts_source_md5) != 1L || is.na(prov$counts_source_md5) ||
      !grepl("^[0-9a-f]{32}$", prov$counts_source_md5) ||
      !prov$counts_backend %in% c("ondisc_odm", "scperturb_h5ad") ||
      is.null(names(prov$input_md5)) ||
      is.null(names(prov$package_versions))) {
    stop("Manifest hashes or provenance are incomplete: ", basename(file))
  }
  if (!identical(cfg$analysis, expected_analysis) ||
      !identical(cfg$scan_mode, expected_scan_mode) ||
      !identical(as.integer(cfg$min_treated), min_treated) ||
      !identical(as.integer(cfg$permutations), B) ||
      !identical(as.integer(cfg$seed), seed) ||
      !identical(as.integer(cfg$response_panel_size), 5000L) ||
      !identical(cfg$primary_test,
                 "two-sided absolute-statistic global label permutation")) {
    stop("Configuration mismatch: ", basename(file))
  }
  if (mode == "pilot" &&
      (!identical(as.integer(cfg$declared_pilot_pairs), 5000L) ||
       is.null(shard$pilot_manifest))) {
    stop("Pilot declaration missing or invalid: ", basename(file))
  }
  if (!is.data.frame(pairs) || !is.data.frame(declared) ||
      nrow(pairs) != nrow(declared) || nrow(pairs) < 1L) {
    stop("Pair table/manifest size mismatch: ", basename(file))
  }
  required_pair_fields <- c(
    "pair_id", "exposure_index", "exposure_gene", "response_rank",
    "response_gene", "pilot_pair_index", "n_treated",
    "n_treated_batches", "max_treated_batch_share", "global_mean_count",
    "global_nonzero"
  )
  required_declared_fields <- c(
    "exposure_index", "exposure_gene", "response_rank", "response_gene",
    "pilot_pair_index"
  )
  if (!all(required_pair_fields %in% names(pairs)) ||
      !all(required_declared_fields %in% names(declared))) {
    stop("Pair identity fields are incomplete: ", basename(file))
  }
  declared_ids <- paste(declared$exposure_gene, declared$response_gene,
                        sep = "::")
  declared_match <- match(pairs$pair_id, declared_ids)
  declared_row <- declared[declared_match, , drop = FALSE]
  if (anyDuplicated(pairs$pair_id) || anyDuplicated(declared_ids) ||
      anyNA(declared_match) ||
      !identical(as.integer(pairs$exposure_index),
                 as.integer(declared_row$exposure_index)) ||
      !identical(as.character(pairs$exposure_gene),
                 as.character(declared_row$exposure_gene)) ||
      !identical(as.integer(pairs$response_rank),
                 as.integer(declared_row$response_rank)) ||
      !identical(as.character(pairs$response_gene),
                 as.character(declared_row$response_gene)) ||
      !identical(as.integer(pairs$pilot_pair_index),
                 as.integer(declared_row$pilot_pair_index)) ||
      any(pairs$exposure_gene == pairs$response_gene)) {
    stop("Pair identity audit failed: ", basename(file))
  }
  exposure_match <- match(
    as.integer(pairs$exposure_index),
    as.integer(shard$exposure_manifest$exposure_index)
  )
  response_match <- match(
    as.integer(pairs$response_rank),
    as.integer(shard$response_manifest$response_rank)
  )
  exposure_row <- shard$exposure_manifest[exposure_match, , drop = FALSE]
  response_row <- shard$response_manifest[response_match, , drop = FALSE]
  if (anyNA(exposure_match) || anyNA(response_match) ||
      !identical(as.character(pairs$exposure_gene),
                 as.character(exposure_row$exposure_gene)) ||
      !identical(as.character(pairs$response_gene),
                 as.character(response_row$response_gene)) ||
      !identical(as.integer(pairs$n_treated),
                 as.integer(exposure_row$n_treated)) ||
      !identical(as.integer(pairs$n_treated_batches),
                 as.integer(exposure_row$n_treated_batches)) ||
      !identical(as.numeric(pairs$max_treated_batch_share),
                 as.numeric(exposure_row$max_treated_batch_share)) ||
      !identical(as.numeric(pairs$global_mean_count),
                 as.numeric(response_row$global_mean_count)) ||
      !identical(as.integer(pairs$global_nonzero),
                 as.integer(response_row$global_nonzero))) {
    stop("Pair rows do not match embedded exposure/response manifests: ",
         basename(file))
  }
  if (mode == "pilot") {
    pilot_ids_in_shard <- paste(
      shard$pilot_manifest$exposure_gene,
      shard$pilot_manifest$response_gene, sep = "::"
    )
    pilot_match <- match(pairs$pair_id, pilot_ids_in_shard)
    pilot_row <- shard$pilot_manifest[pilot_match, , drop = FALSE]
    if (anyNA(pilot_match) ||
        !identical(as.integer(pairs$exposure_index),
                   as.integer(pilot_row$exposure_index)) ||
        !identical(as.integer(pairs$response_rank),
                   as.integer(pilot_row$response_rank)) ||
        !identical(as.integer(pairs$pilot_pair_index),
                   as.integer(pilot_row$pilot_pair_index))) {
      stop("Pair rows do not match the frozen pilot manifest: ",
           basename(file))
    }
  }
  ## Pairs whose NB nuisance fit degenerates are excluded by the pre-declared,
  ## label-invariant rule: that fit is a function of (Y, Z_design) alone and is
  ## held fixed across the permutation orbit, so excluding them cannot depend
  ## on the exposure assignment.  Every other non-ok status is fatal -- in
  ## particular a degenerate OBSERVED score depends on X and must never be
  ## dropped on that basis.
  excluded <- pairs$status == "degenerate_nuisance_fit"
  if (any(!pairs$status[!excluded] %in% "ok")) {
    stop("Shard contains a non-ok pair outside the declared exclusion rule: ",
         basename(file), " (",
         paste(unique(pairs$status[!excluded & pairs$status != "ok"]),
               collapse = ", "), ")")
  }
  shard_excluded <- pairs[excluded, c("pair_id", "pilot_pair_index",
                                      "exposure_index", "exposure_gene",
                                      "response_gene", "response_rank",
                                      "status", "error_message"),
                          drop = FALSE]
  pairs <- pairs[!excluded, , drop = FALSE]
  if (nrow(pairs) < 1L) {
    stop("Every pair in the shard was excluded: ", basename(file))
  }
  if (anyNA(pairs$p_two_sided_barry) || anyNA(pairs$p_two_sided_empirical)) {
    stop("Shard contains a missing-p pair: ", basename(file))
  }
  pair_index <- (as.double(pairs$exposure_index) - 1) * 5000 +
    as.double(pairs$response_rank)
  expected_pair_seed <- as.integer(
    (as.double(seed) + 10007 * pair_index) %% 2147483646 + 1
  )
  if (!identical(as.integer(pairs$pair_seed), expected_pair_seed) ||
      any(pairs$hits_abs_barry < 0 | pairs$hits_abs_barry > B) ||
      any(pairs$hits_abs_empirical < 0 | pairs$hits_abs_empirical > B) ||
      max(abs(pairs$p_two_sided_barry -
                (1 + pairs$hits_abs_barry) / (B + 1))) > 1e-15 ||
      max(abs(pairs$p_two_sided_empirical -
                (1 + pairs$hits_abs_empirical) / (B + 1))) > 1e-15 ||
      any(sign(pairs$information_stat) != sign(pairs$empirical_stat))) {
    stop("Seed/hit/p-value audit failed: ", basename(file))
  }

  signature <- list(
    manifest_md5 = checks$manifest_md5,
    exposure_manifest = shard$exposure_manifest,
    response_manifest = shard$response_manifest,
    pilot_manifest = if (mode == "pilot") shard$pilot_manifest else NULL,
    counts_backend = prov$counts_backend,
    counts_source_file = prov$counts_source_file,
    counts_source_size_bytes = prov$counts_source_size_bytes,
    counts_source_mtime_utc = prov$counts_source_mtime_utc,
    counts_source_md5 = prov$counts_source_md5,
    counts_identity_verification = prov$counts_identity_verification,
    input_md5 = prov$input_md5,
    runner_script_md5 = prov$script_md5,
    score_helper_md5 = prov$score_helper_md5,
    package_versions = prov$package_versions,
    rng_kind = cfg$rng_kind,
    permutation_batch_size = cfg$permutation_batch_size
  )
  if (is.null(baseline)) {
    baseline <- signature
  } else {
    if (!same_named_atomic(signature$manifest_md5,
                           baseline$manifest_md5) ||
        !identical(signature$exposure_manifest,
                   baseline$exposure_manifest) ||
        !identical(signature$response_manifest,
                   baseline$response_manifest) ||
        !identical(signature$pilot_manifest, baseline$pilot_manifest) ||
        !identical(signature$counts_backend, baseline$counts_backend) ||
        !identical(signature$counts_source_file,
                   baseline$counts_source_file) ||
        !identical(signature$counts_source_size_bytes,
                   baseline$counts_source_size_bytes) ||
        !identical(signature$counts_source_mtime_utc,
                   baseline$counts_source_mtime_utc) ||
        !identical(signature$counts_source_md5,
                   baseline$counts_source_md5) ||
        !identical(signature$counts_identity_verification,
                   baseline$counts_identity_verification) ||
        !same_named_atomic(signature$input_md5, baseline$input_md5) ||
        !identical(signature$runner_script_md5,
                   baseline$runner_script_md5) ||
        !identical(signature$score_helper_md5,
                   baseline$score_helper_md5) ||
        !same_named_atomic(signature$package_versions,
                           baseline$package_versions) ||
        !identical(signature$rng_kind, baseline$rng_kind) ||
        !identical(signature$permutation_batch_size,
                   baseline$permutation_batch_size)) {
      stop("Manifest or provenance mismatch: ", basename(file))
    }
  }
  pair_chunks[[i]] <- pairs
  excluded_chunks[[i]] <- shard_excluded
  shard_rows[[i]] <- data.frame(
    shard_file = basename(file), n_pairs = nrow(pairs),
    n_excluded = nrow(shard_excluded),
    exposure_index_min = min(pairs$exposure_index),
    exposure_index_max = max(pairs$exposure_index),
    response_rank_min = min(pairs$response_rank),
    response_rank_max = max(pairs$response_rank),
    stringsAsFactors = FALSE
  )
  rm(shard, pairs, declared, shard_excluded)
  if (i %% 25L == 0L || i == length(candidate)) {
    cat("  validated ", i, "/", length(candidate), " shards\n", sep = "")
  }
}

manifest_filenames <- c(
  expression_audit = "frangieh-exploratory-expression-audit.csv",
  response_universe = "frangieh-exploratory-response-universe.csv",
  selected_responses = "frangieh-exploratory-selected-responses.csv",
  exposures = sprintf(
    "frangieh-exploratory-exposures-min-treated-%02d.csv", min_treated
  ),
  association_blind_pilot_5000 =
    "frangieh-exploratory-association-blind-pilot-5000.csv"
)
manifest_files <- file.path(
  input_dir, "manifests", manifest_filenames[expected_manifest_names]
)
if (any(!file.exists(manifest_files))) {
  stop("Shared manifest file(s) are absent: ", paste(
    basename(manifest_files[!file.exists(manifest_files)]), collapse = ", "
  ))
}
disk_manifest_md5 <- unname(tools::md5sum(manifest_files))
names(disk_manifest_md5) <- expected_manifest_names
if (!same_named_atomic(disk_manifest_md5, baseline$manifest_md5)) {
  stop("On-disk manifest hashes do not match embedded shard provenance")
}
locked_expression_audit <- utils::read.csv(
  manifest_files[[which(expected_manifest_names == "expression_audit")]],
  stringsAsFactors = FALSE
)
if (nrow(locked_expression_audit) != 1L ||
    locked_expression_audit$schema_version != 3L ||
    !identical(locked_expression_audit$backend,
               baseline$counts_backend) ||
    !identical(locked_expression_audit$counts_source_path,
               baseline$counts_source_file) ||
    !identical(as.numeric(locked_expression_audit$counts_source_size_bytes),
               as.numeric(baseline$counts_source_size_bytes)) ||
    !identical(locked_expression_audit$counts_source_mtime_utc,
               baseline$counts_source_mtime_utc) ||
    !identical(locked_expression_audit$counts_source_md5,
               baseline$counts_source_md5)) {
  stop("Locked expression audit and shard counts-source identity disagree")
}

pair_results <- do.call(rbind, pair_chunks)
excluded_pairs <- do.call(rbind, excluded_chunks)
if (is.null(excluded_pairs)) {
  excluded_pairs <- data.frame(
    pair_id = character(), pilot_pair_index = integer(),
    exposure_index = integer(), exposure_gene = character(),
    response_gene = character(), response_rank = integer(),
    status = character(), error_message = character(),
    stringsAsFactors = FALSE
  )
}
rownames(excluded_pairs) <- NULL
shard_inventory <- do.call(rbind, shard_rows)
rownames(pair_results) <- NULL
if (anyDuplicated(pair_results$pair_id)) {
  stop("Duplicate pair IDs across shards; refusing to combine")
}

eligible_exposures <- baseline$exposure_manifest[
  baseline$exposure_manifest$eligible, , drop = FALSE
]
eligible_exposures <- eligible_exposures[order(
  eligible_exposures$exposure_index
), ]
responses <- baseline$response_manifest[order(
  baseline$response_manifest$response_rank
), ]
expected_exposures <- if (min_treated == 10L) 241L else 244L
if (nrow(eligible_exposures) != expected_exposures ||
    nrow(responses) != 5000L || any(responses$is_perturbation_target)) {
  stop("Baseline exposure/response manifest failed its size/control audit")
}

if (mode == "pilot") {
  pilot <- baseline$pilot_manifest[order(
    baseline$pilot_manifest$pilot_pair_index
  ), ]
  pilot_ids <- paste(pilot$exposure_gene, pilot$response_gene, sep = "::")
  exposure_balance <- table(pilot$exposure_index)
  covered_ids <- c(pair_results$pair_id, excluded_pairs$pair_id)
  covered_pilot_index <- c(as.integer(pair_results$pilot_pair_index),
                           as.integer(excluded_pairs$pilot_pair_index))
  if (length(covered_ids) != 5000L || nrow(pilot) != 5000L ||
      anyDuplicated(pilot$pilot_pair_index) ||
      !identical(sort(as.integer(pilot$pilot_pair_index)), seq_len(5000L)) ||
      anyDuplicated(pilot$response_rank) ||
      !identical(sort(as.integer(pilot$response_rank)), seq_len(5000L)) ||
      length(exposure_balance) != 241L || min(exposure_balance) != 20L ||
      max(exposure_balance) != 21L ||
      !setequal(covered_ids, pilot_ids) ||
      anyNA(covered_pilot_index) ||
      !identical(sort(covered_pilot_index), seq_len(5000L))) {
    stop("Frozen-pilot coverage/balance is incomplete; refusing to combine")
  }
  pair_results <- pair_results[order(pair_results$pilot_pair_index), ]
  pilot_row <- pilot[match(pair_results$pair_id, pilot_ids), , drop = FALSE]
  if (!identical(as.integer(pair_results$pilot_pair_index),
                 as.integer(pilot_row$pilot_pair_index)) ||
      !identical(as.integer(pair_results$exposure_index),
                 as.integer(pilot_row$exposure_index)) ||
      !identical(as.integer(pair_results$response_rank),
                 as.integer(pilot_row$response_rank))) {
    stop("Frozen-pilot row mapping failed after union; refusing to combine")
  }
} else {
  expected_n <- expected_exposures * 5000L
  grid_index <- c(as.integer(pair_results$exposure_index),
                  as.integer(excluded_pairs$exposure_index))
  exposure_counts <- table(grid_index)
  response_counts <- table(c(as.integer(pair_results$response_rank),
                             as.integer(excluded_pairs$response_rank)))
  if (nrow(pair_results) + nrow(excluded_pairs) != expected_n ||
      length(exposure_counts) != expected_exposures ||
      any(exposure_counts != 5000L) ||
      !identical(sort(as.integer(names(exposure_counts))),
                 eligible_exposures$exposure_index) ||
      length(response_counts) != 5000L ||
      any(response_counts != expected_exposures) ||
      !identical(sort(as.integer(names(response_counts))), seq_len(5000L))) {
    stop("Full-grid coverage is incomplete; refusing to combine")
  }
  pair_results <- pair_results[order(pair_results$exposure_index,
                                     pair_results$response_rank), ]
}
rownames(pair_results) <- NULL

alphas <- c(0.001, 0.005, 0.01, 0.05, 0.10)
paired_summary <- do.call(rbind, lapply(alphas, function(alpha) {
  reject_barry <- pair_results$p_two_sided_barry <= alpha
  reject_empirical <- pair_results$p_two_sided_empirical <= alpha
  barry_only <- sum(reject_barry & !reject_empirical)
  empirical_only <- sum(!reject_barry & reject_empirical)
  discordant <- barry_only + empirical_only
  exact_p <- if (discordant == 0L) {
    1
  } else {
    min(1, 2 * stats::pbinom(
      min(barry_only, empirical_only), discordant, prob = 0.5
    ))
  }
  data.frame(
    alpha = alpha,
    n_tests = nrow(pair_results),
    n_barry = sum(reject_barry),
    n_empirical = sum(reject_empirical),
    n_both = sum(reject_barry & reject_empirical),
    n_barry_only = barry_only,
    n_empirical_only = empirical_only,
    n_neither = sum(!reject_barry & !reject_empirical),
    empirical_minus_barry = sum(reject_empirical) - sum(reject_barry),
    exact_paired_discordance_p = exact_p,
    stringsAsFactors = FALSE
  )
}))

plot_ready <- pair_results[, c(
  "pair_id", "pilot_pair_index", "exposure_index", "exposure_gene",
  "response_rank", "response_gene", "n_treated", "global_mean_count",
  "global_nonzero", "information_stat", "empirical_stat",
  "p_two_sided_barry", "p_two_sided_empirical"
)]
plot_ready$score_direction <- ifelse(
  plot_ready$information_stat < 0, "negative",
  ifelse(plot_ready$information_stat > 0, "positive", "zero")
)
plot_ready$neglog10_p_barry <- -log10(plot_ready$p_two_sided_barry)
plot_ready$neglog10_p_empirical <- -log10(plot_ready$p_two_sided_empirical)
plot_ready$empirical_minus_barry_neglog10 <-
  plot_ready$neglog10_p_empirical - plot_ready$neglog10_p_barry
plot_ready$smaller_raw_p <- ifelse(
  plot_ready$p_two_sided_barry < plot_ready$p_two_sided_empirical, "barry",
  ifelse(plot_ready$p_two_sided_empirical < plot_ready$p_two_sided_barry,
         "empirical", "tie")
)

min_raw_p <- 1 / (B + 1)
resolution_note <- paste0(
  "Raw permutation p-values only. No BH/q-values computed. B=", B,
  " implies minimum attainable p=", format(min_raw_p, digits = 8),
  "; this resolution is not treated as sufficient evidence for ",
  "multiplicity-controlled discovery over ", nrow(pair_results), " tests."
)

tag <- paste("frangieh-exploratory-trans-combined", mode,
             paste0("min", min_treated), paste0("B", B), sep = "-")
rds_path <- file.path(output_dir, paste0(tag, ".rds"))
csv_path <- file.path(output_dir, paste0(tag, ".csv"))
summary_path <- file.path(output_dir,
                          paste0(tag, "-raw-p-comparison-summary.csv"))
plot_path <- file.path(output_dir, paste0(tag, "-plot-ready.csv"))
destinations <- c(rds_path, csv_path, summary_path, plot_path)
if (any(file.exists(destinations)) && !overwrite) {
  stop("Combined output exists; use --overwrite: ",
       paste(basename(destinations[file.exists(destinations)]), collapse = ", "))
}

combined <- list(
  config = list(
    analysis = "frangieh_exploratory_trans_strict_combination",
    mode = mode, min_treated = min_treated, permutations = B, seed = seed,
    counts_backend = baseline$counts_backend,
    n_tests = nrow(pair_results), n_shards = length(candidate),
    n_declared_pairs = nrow(pair_results) + nrow(excluded_pairs),
    n_excluded_degenerate_nuisance_fit = nrow(excluded_pairs),
    exclusion_rule = paste(
      "pairs whose negative-binomial nuisance fit on (Y, Z_design)",
      "degenerates are excluded; that fit is held fixed across the",
      "permutation orbit, so the criterion is invariant to every",
      "guide-label permutation and cannot depend on the exposure"
    ),
    p_value_resolution = min_raw_p,
    multiplicity_correction = "not_computed",
    multiplicity_note = resolution_note
  ),
  manifest_md5 = baseline$manifest_md5,
  exposure_manifest = baseline$exposure_manifest,
  response_manifest = baseline$response_manifest,
  pilot_manifest = baseline$pilot_manifest,
  shard_inventory = shard_inventory,
  excluded_pairs = excluded_pairs,
  shard_files = basename(candidate),
  shard_md5 = unname(tools::md5sum(candidate)),
  pair_results = pair_results,
  raw_p_comparison_summary = paired_summary,
  plot_ready = plot_ready,
  provenance = list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    command = commandArgs(),
    combiner_script = script_path,
    combiner_script_md5 = unname(tools::md5sum(script_path)),
    runner_script_md5 = baseline$runner_script_md5,
    score_helper_md5 = baseline$score_helper_md5,
    counts_backend = baseline$counts_backend,
    counts_source_file = baseline$counts_source_file,
    counts_source_size_bytes = baseline$counts_source_size_bytes,
    counts_source_mtime_utc = baseline$counts_source_mtime_utc,
    counts_source_md5 = baseline$counts_source_md5,
    counts_identity_verification = baseline$counts_identity_verification,
    input_md5 = baseline$input_md5,
    package_versions = baseline$package_versions,
    session_info = utils::sessionInfo()
  )
)
names(combined$shard_md5) <- basename(candidate)

tmp_rds <- tempfile(pattern = paste0(tag, "."), tmpdir = output_dir,
                    fileext = ".rds.tmp")
tmp_csv <- tempfile(pattern = paste0(tag, "."), tmpdir = output_dir,
                    fileext = ".csv.tmp")
tmp_summary <- tempfile(pattern = paste0(tag, "."), tmpdir = output_dir,
                        fileext = ".summary.csv.tmp")
tmp_plot <- tempfile(pattern = paste0(tag, "."), tmpdir = output_dir,
                     fileext = ".plot.csv.tmp")
saveRDS(combined, tmp_rds, compress = "gzip")
utils::write.csv(pair_results, tmp_csv, row.names = FALSE, na = "")
utils::write.csv(paired_summary, tmp_summary, row.names = FALSE, na = "")
utils::write.csv(plot_ready, tmp_plot, row.names = FALSE, na = "")
atomic_publish(tmp_rds, rds_path, overwrite)
atomic_publish(tmp_csv, csv_path, overwrite)
atomic_publish(tmp_summary, summary_path, overwrite)
atomic_publish(tmp_plot, plot_path, overwrite)

cat("Coverage and provenance audit passed.\n")
cat("  declared pairs: ", nrow(pair_results) + nrow(excluded_pairs), "\n",
    sep = "")
cat("  combined tests: ", nrow(pair_results), "\n", sep = "")
cat("  excluded (degenerate nuisance fit): ", nrow(excluded_pairs), "\n",
    sep = "")
if (nrow(excluded_pairs)) {
  print(excluded_pairs[, c("pair_id", "response_rank", "error_message")],
        row.names = FALSE)
}
cat("  RDS: ", rds_path, "\n", sep = "")
cat("  CSV: ", csv_path, "\n", sep = "")
cat("  summary: ", summary_path, "\n", sep = "")
cat("  plot-ready: ", plot_path, "\n", sep = "")
cat("  NOTE: ", resolution_note, "\n", sep = "")
