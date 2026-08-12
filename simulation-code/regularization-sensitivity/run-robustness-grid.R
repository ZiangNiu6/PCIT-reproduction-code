#!/usr/bin/env Rscript

# Extend the Section 4.2 locked-tuning experiment over a prespecified grid.
#
# The baseline is q = 178 and V = 7.5.  We cross three penalized dimensions
# with half, baseline, and twice the baseline signal variance.  In every cell,
# X- and Y-side penalties are selected separately by minimum LOOCV (the +1SE
# diagnostics are retained), and the selected X-side penalty is held fixed
# over complete raw-X permutations.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("Run this file with Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", script_argument))
script_directory <- dirname(script_file)
project_root <- normalizePath(file.path(script_directory, "..", "..", ".."))
results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
regularization_output_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)
simulation_directory <- script_directory
output_directory <- Sys.getenv(
  "PCIT_ROBUSTNESS_OUTPUT",
  unset = file.path(regularization_output_directory, "robustness-grid")
)

parse_positive_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(name, " must be a positive integer.")
  }
  value
}

parse_numeric_grid <- function(name, default) {
  value <- as.numeric(strsplit(
    Sys.getenv(name, unset = default),
    ",",
    fixed = TRUE
  )[[1L]])
  if (length(value) < 1L || any(!is.finite(value))) {
    stop(name, " must be a comma-separated numeric grid.")
  }
  value
}

parse_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(
    name,
    unset = if (default) "true" else "false"
  ))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop(name, " must be true or false.")
  }
  value %in% c("true", "1", "yes")
}

R <- parse_positive_integer("PCIT_ROBUSTNESS_R", "1000")
B <- parse_positive_integer("PCIT_ROBUSTNESS_B", "499")
signal_rank <- parse_positive_integer("PCIT_ROBUSTNESS_SIGNAL_RANK", "16")
q_values <- as.integer(parse_numeric_grid(
  "PCIT_ROBUSTNESS_Q",
  "98,178,238"
))
signal_variances <- parse_numeric_grid(
  "PCIT_ROBUSTNESS_SIGNAL_VARIANCE",
  "3.75,7.5,15"
)
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(detected_cores) != 1L ||
    !is.finite(detected_cores) ||
    detected_cores < 1L) {
  detected_cores <- 1L
}
workers <- parse_positive_integer(
  "PCIT_ROBUSTNESS_WORKERS",
  as.character(min(4L, detected_cores))
)
cell_cores <- parse_positive_integer("PCIT_ROBUSTNESS_CELL_CORES", "1")
force <- parse_flag("PCIT_ROBUSTNESS_FORCE", FALSE)

if (any(q_values < signal_rank) ||
    any(q_values %% 2L != 0L) ||
    any(q_values >= 298L)) {
  stop("Each q must be even, at least signal_rank, and smaller than 298.")
}
if (signal_rank %% 2L != 0L) {
  stop("PCIT_ROBUSTNESS_SIGNAL_RANK must be even.")
}
if (any(signal_variances < 0)) {
  stop("Signal variances must be nonnegative.")
}
if (B < 19L) {
  stop("PCIT_ROBUSTNESS_B must be at least 19.")
}

format_value <- function(value) {
  gsub("\\.", "p", format(value, scientific = FALSE, trim = TRUE))
}

cell_grid <- expand.grid(
  learner = c("LRR", "KRR"),
  q_penalized = q_values,
  signal_variance = signal_variances,
  stringsAsFactors = FALSE
)
cell_grid$cell_id <- paste0(
  tolower(cell_grid$learner),
  "-q", cell_grid$q_penalized,
  "-v", vapply(cell_grid$signal_variance, format_value, character(1))
)

lrr_script <- file.path(
  simulation_directory,
  "lrr-selected-calibration.R"
)
krr_script <- file.path(
  simulation_directory,
  "krr-selected-calibration.R"
)
if (!file.exists(lrr_script) || !file.exists(krr_script)) {
  stop("Could not find the LRR/KRR confirmation scripts.")
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
cell_directory <- file.path(output_directory, "cells")
log_directory <- file.path(output_directory, "logs")
dir.create(cell_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(log_directory, recursive = TRUE, showWarnings = FALSE)

expected_files <- function(cell) {
  directory <- file.path(cell_directory, cell$cell_id)
  if (cell$learner == "LRR") {
    c(
      summary = file.path(directory, "lrr-calibration-summary.csv"),
      diagnostics = file.path(directory, "lrr-tuning-diagnostics.csv"),
      results = file.path(directory, "lrr-confirmation-results.rds")
    )
  } else {
    prefix <- file.path(directory, "krr")
    c(
      summary = paste0(prefix, "-calibration-summary.csv"),
      diagnostics = paste0(prefix, "-tuning-diagnostics.csv"),
      results = paste0(prefix, "-results.rds")
    )
  }
}

cell_is_complete <- function(cell) {
  paths <- expected_files(cell)
  if (!all(file.exists(paths))) {
    return(FALSE)
  }
  result <- tryCatch(readRDS(paths[["results"]]), error = function(e) NULL)
  if (is.null(result)) {
    return(FALSE)
  }
  design <- result$design
  result_R <- if (!is.null(result$R)) result$R else design$R
  result_B <- if (!is.null(result$B)) result$B else design$B
  isTRUE(design$q_penalized == cell$q_penalized) &&
    isTRUE(abs(design$signal_variance - cell$signal_variance) < 1e-12) &&
    isTRUE(design$signal_rank == signal_rank) &&
    isTRUE(result_R == R) &&
    isTRUE(result_B == B)
}

run_cell <- function(index) {
  cell <- cell_grid[index, , drop = FALSE]
  paths <- expected_files(cell)
  directory <- dirname(paths[["results"]])
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_directory, paste0(cell$cell_id, ".log"))

  if (!force && cell_is_complete(cell)) {
    return(data.frame(
      cell_id = cell$cell_id,
      status = "reused",
      log_file = log_file,
      stringsAsFactors = FALSE
    ))
  }

  if (cell$learner == "LRR") {
    program <- lrr_script
    environment <- c(
      paste0("PCIT_LRR_Q_PENALIZED=", cell$q_penalized),
      paste0("PCIT_LRR_SIGNAL_VARIANCE=", cell$signal_variance),
      paste0("PCIT_LRR_SIGNAL_RANK=", signal_rank),
      "PCIT_LRR_SIGNAL_PATTERN=leading",
      paste0("PCIT_LRR_R=", R),
      paste0("PCIT_LRR_B=", B),
      paste0("PCIT_LRR_CORES=", cell_cores),
      paste0("PCIT_LRR_OUTPUT=", directory)
    )
  } else {
    program <- krr_script
    prefix <- file.path(directory, "krr")
    environment <- c(
      paste0("PCIT_KRR_CONFIRM_Q_PENALIZED=", cell$q_penalized),
      paste0(
        "PCIT_KRR_CONFIRM_SIGNAL_VARIANCE=",
        cell$signal_variance
      ),
      paste0("PCIT_KRR_CONFIRM_SIGNAL_RANK=", signal_rank),
      "PCIT_KRR_CONFIRM_SIGNAL_PATTERN=leading",
      paste0("PCIT_KRR_CONFIRM_R=", R),
      paste0("PCIT_KRR_CONFIRM_B=", B),
      paste0("PCIT_KRR_CONFIRM_CORES=", cell_cores),
      paste0("PCIT_KRR_CONFIRM_OUTPUT_PREFIX=", prefix)
    )
  }

  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", program),
    stdout = log_file,
    stderr = log_file,
    env = environment
  )
  if (status != 0L || !cell_is_complete(cell)) {
    stop(
      "Simulation cell ", cell$cell_id,
      " failed; inspect ", log_file, "."
    )
  }
  data.frame(
    cell_id = cell$cell_id,
    status = "completed",
    log_file = log_file,
    stringsAsFactors = FALSE
  )
}

message(
  "Running ", nrow(cell_grid), " robustness cells with R=", R,
  ", B=", B, ", workers=", workers,
  ", and cell_cores=", cell_cores, "."
)
indices <- seq_len(nrow(cell_grid))
if (.Platform$OS.type == "unix" && workers > 1L) {
  run_results <- parallel::mclapply(
    indices,
    run_cell,
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  run_results <- lapply(indices, run_cell)
}
failed <- vapply(run_results, inherits, logical(1), "try-error")
if (any(failed)) {
  stop(
    "At least one robustness cell failed: ",
    paste(cell_grid$cell_id[failed], collapse = ", ")
  )
}
run_manifest <- do.call(rbind, run_results)
write.csv(
  merge(cell_grid, run_manifest, by = "cell_id", sort = FALSE),
  file.path(output_directory, "run-manifest.csv"),
  row.names = FALSE
)

normalize_summary <- function(cell) {
  paths <- expected_files(cell)
  source <- read.csv(paths[["summary"]], check.names = FALSE)
  if (cell$learner == "LRR") {
    data.frame(
      learner = "LRR",
      q_penalized = cell$q_penalized,
      q_total = cell$q_penalized + 2L,
      hard_residual_rank = 298L - cell$q_penalized,
      signal_variance = cell$signal_variance,
      signal_rank = signal_rank,
      selection_rule = source$selection_rule,
      calibration = source$calibration,
      rejections = source$rejections,
      R = source$R,
      B = source$B,
      estimate = source$estimate,
      lower_95 = source$lower_95_exact,
      upper_95 = source$upper_95_exact,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      learner = "KRR",
      q_penalized = cell$q_penalized,
      q_total = cell$q_penalized + 2L,
      hard_residual_rank = 298L - cell$q_penalized,
      signal_variance = cell$signal_variance,
      signal_rank = signal_rank,
      selection_rule = source$cv_rule,
      calibration = source$method,
      rejections = source$rejections,
      R = source$R,
      B = source$B,
      estimate = source$estimate,
      lower_95 = source$lower_95,
      upper_95 = source$upper_95,
      stringsAsFactors = FALSE
    )
  }
}

normalize_diagnostics <- function(cell) {
  paths <- expected_files(cell)
  source <- read.csv(paths[["diagnostics"]], check.names = FALSE)
  result <- readRDS(paths[["results"]])
  if (cell$learner == "LRR") {
    rule <- source$selection_rule
    mean_x <- source$mean_observed_x_effective_dimension
    median_x <- source$median_observed_x_effective_dimension
    mean_y <- source$mean_observed_y_effective_dimension
    median_y <- source$median_observed_y_effective_dimension
    mean_gamma <- source$mean_gamma_plugin
    mean_alpha <- source$mean_alpha_normal_plugin
    raw <- result$raw
    raw_rule <- raw$selection_rule
    x_candidate <- raw$observed_x_candidate
    y_candidate <- raw$observed_y_candidate
    penalty <- result$candidate_grid$penalty_absolute
    eigenvalues <- exp(seq(
      log(4 * 300),
      log(0.25 * 300),
      length.out = cell$q_penalized
    ))
  } else {
    rule <- source$cv_rule
    mean_x <- source$mean_observed_x_effective_dimension
    median_x <- source$median_observed_x_effective_dimension
    mean_y <- source$mean_observed_y_effective_dimension
    median_y <- source$median_observed_y_effective_dimension
    mean_gamma <- source$mean_gamma_plugin
    mean_alpha <- source$mean_alpha_normal_plugin
    raw <- result$raw_diagnostics
    raw_rule <- raw$cv_rule
    x_candidate <- raw$observed_x_candidate
    y_candidate <- raw$observed_y_candidate
    penalty <- result$candidate_grid$penalty_ratio *
      stats::median(exp(seq(
        log(1),
        log(0.01),
        length.out = cell$q_penalized
      )))
    eigenvalues <- exp(seq(
      log(1),
      log(0.01),
      length.out = cell$q_penalized
    ))
  }

  weights <- vapply(penalty, function(lambda) {
    if (lambda == 0) {
      rep(0, cell$q_penalized)
    } else if (is.infinite(lambda)) {
      rep(1, cell$q_penalized)
    } else {
      lambda / (eigenvalues + lambda)
    }
  }, numeric(cell$q_penalized))
  active_x <- seq.int(1L, signal_rank, by = 2L)
  active_y <- seq.int(2L, signal_rank, by = 2L)
  candidate_x_rms <- sqrt(
    cell$signal_variance / length(active_x) *
      colSums(weights[active_x, , drop = FALSE]^2)
  )
  candidate_y_rms <- sqrt(
    cell$signal_variance / length(active_y) *
      colSums(weights[active_y, , drop = FALSE]^2)
  )
  mean_x_signal_rms <- vapply(rule, function(rule_name) {
    mean(candidate_x_rms[x_candidate[raw_rule == rule_name]])
  }, numeric(1))
  median_x_signal_rms <- vapply(rule, function(rule_name) {
    stats::median(candidate_x_rms[x_candidate[raw_rule == rule_name]])
  }, numeric(1))
  mean_y_signal_rms <- vapply(rule, function(rule_name) {
    mean(candidate_y_rms[y_candidate[raw_rule == rule_name]])
  }, numeric(1))
  median_y_signal_rms <- vapply(rule, function(rule_name) {
    stats::median(candidate_y_rms[y_candidate[raw_rule == rule_name]])
  }, numeric(1))
  data.frame(
    learner = cell$learner,
    q_penalized = cell$q_penalized,
    q_total = cell$q_penalized + 2L,
    hard_residual_rank = 298L - cell$q_penalized,
    signal_variance = cell$signal_variance,
    signal_rank = signal_rank,
    selection_rule = rule,
    mean_x_effective_dimension = mean_x,
    median_x_effective_dimension = median_x,
    mean_y_effective_dimension = mean_y,
    median_y_effective_dimension = median_y,
    mean_x_residual_signal_rms = mean_x_signal_rms,
    median_x_residual_signal_rms = median_x_signal_rms,
    mean_y_residual_signal_rms = mean_y_signal_rms,
    median_y_residual_signal_rms = median_y_signal_rms,
    mean_gamma_plugin = mean_gamma,
    mean_alpha_normal_plugin = mean_alpha,
    stringsAsFactors = FALSE
  )
}

normalize_paired_difference <- function(cell) {
  paths <- expected_files(cell)
  result <- readRDS(paths[["results"]])
  raw <- if (cell$learner == "LRR") {
    result$raw
  } else {
    result$raw_diagnostics
  }
  rule_column <- if (cell$learner == "LRR") {
    "selection_rule"
  } else {
    "cv_rule"
  }
  rules <- c("Minimum CV", "+1SE")
  do.call(rbind, lapply(rules, function(rule_name) {
    rule_data <- raw[raw[[rule_column]] == rule_name, , drop = FALSE]
    difference <- as.integer(rule_data$normal_reject) -
      as.integer(rule_data$locked_reject)
    estimate <- mean(difference)
    standard_error <- stats::sd(difference) / sqrt(length(difference))
    data.frame(
      learner = cell$learner,
      q_penalized = cell$q_penalized,
      signal_variance = cell$signal_variance,
      selection_rule = rule_name,
      normal_minus_permutation = estimate,
      paired_se = standard_error,
      paired_lower_95 = max(-1, estimate - 1.96 * standard_error),
      paired_upper_95 = min(1, estimate + 1.96 * standard_error),
      normal_only_rejections = sum(difference == 1L),
      permutation_only_rejections = sum(difference == -1L),
      R = nrow(rule_data),
      stringsAsFactors = FALSE
    )
  }))
}

summaries <- do.call(rbind, lapply(seq_len(nrow(cell_grid)), function(index) {
  normalize_summary(cell_grid[index, , drop = FALSE])
}))
diagnostics <- do.call(rbind, lapply(
  seq_len(nrow(cell_grid)),
  function(index) normalize_diagnostics(cell_grid[index, , drop = FALSE])
))
paired_differences <- do.call(rbind, lapply(
  seq_len(nrow(cell_grid)),
  function(index) {
    normalize_paired_difference(cell_grid[index, , drop = FALSE])
  }
))

summaries <- summaries[order(
  summaries$learner,
  summaries$q_penalized,
  summaries$signal_variance,
  match(summaries$selection_rule, c("Minimum CV", "+1SE")),
  summaries$calibration
), ]
rownames(summaries) <- NULL
diagnostics <- diagnostics[order(
  diagnostics$learner,
  diagnostics$q_penalized,
  diagnostics$signal_variance,
  match(diagnostics$selection_rule, c("Minimum CV", "+1SE"))
), ]
rownames(diagnostics) <- NULL

write.csv(
  summaries,
  file.path(output_directory, "robustness-grid-summary.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(output_directory, "robustness-grid-diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  paired_differences,
  file.path(output_directory, "paired-calibration-differences.csv"),
  row.names = FALSE
)

minimum <- summaries[summaries$selection_rule == "Minimum CV", ]
normal <- minimum[minimum$calibration == "Normal calibration", ]
permutation <- minimum[
  minimum$calibration == "Locked-tuning refit",
]
keys <- c("learner", "q_penalized", "signal_variance")
comparison <- merge(
  normal,
  permutation,
  by = keys,
  suffixes = c("_normal", "_permutation")
)
comparison$normal_absolute_error <- abs(comparison$estimate_normal - 0.05)
comparison$permutation_absolute_error <- abs(
  comparison$estimate_permutation - 0.05
)
comparison$absolute_error_reduction <-
  comparison$normal_absolute_error - comparison$permutation_absolute_error
comparison$permutation_ci_contains_nominal <-
  comparison$lower_95_permutation <= 0.05 &
  comparison$upper_95_permutation >= 0.05
comparison$permutation_not_materially_liberal <-
  comparison$estimate_permutation <= 0.07
comparison <- merge(
  comparison,
  paired_differences[
    paired_differences$selection_rule == "Minimum CV",
    c(
      "learner",
      "q_penalized",
      "signal_variance",
      "normal_minus_permutation",
      "paired_se",
      "paired_lower_95",
      "paired_upper_95",
      "normal_only_rejections",
      "permutation_only_rejections"
    )
  ],
  by = keys,
  sort = FALSE
)
comparison <- comparison[order(
  comparison$learner,
  comparison$q_penalized,
  comparison$signal_variance
), ]
write.csv(
  comparison,
  file.path(output_directory, "minimum-cv-comparison.csv"),
  row.names = FALSE
)

aggregate_metrics <- do.call(rbind, lapply(
  split(comparison, comparison$learner),
  function(data) {
    data.frame(
      learner = data$learner[[1L]],
      cells = nrow(data),
      normal_min = min(data$estimate_normal),
      normal_max = max(data$estimate_normal),
      permutation_min = min(data$estimate_permutation),
      permutation_max = max(data$estimate_permutation),
      normal_mean_absolute_error = mean(data$normal_absolute_error),
      permutation_mean_absolute_error = mean(
        data$permutation_absolute_error
      ),
      normal_max_absolute_error = max(data$normal_absolute_error),
      permutation_max_absolute_error = max(
        data$permutation_absolute_error
      ),
      cells_permutation_closer = sum(data$absolute_error_reduction > 0),
      cells_permutation_ci_contains_nominal = sum(
        data$permutation_ci_contains_nominal
      ),
      cells_permutation_not_materially_liberal = sum(
        data$permutation_not_materially_liberal
      ),
      stringsAsFactors = FALSE
    )
  }
))
write.csv(
  aggregate_metrics,
  file.path(output_directory, "aggregate-metrics.csv"),
  row.names = FALSE
)

print(comparison[, c(
  "learner",
  "q_penalized",
  "signal_variance",
  "estimate_normal",
  "estimate_permutation",
  "absolute_error_reduction",
  "permutation_ci_contains_nominal"
)], row.names = FALSE, digits = 4)
print(aggregate_metrics, row.names = FALSE, digits = 4)
message("Saved robustness-grid outputs under ", normalizePath(output_directory))
