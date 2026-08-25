#!/usr/bin/env Rscript

###############################################################################
# Summarize and plot the sparse-exposure and variance-misspecification
# comparison of parametric, permutation, and sign-flip score calibrations.
#
# Usage from anywhere:
#   Rscript --vanilla reproduction-code/simulation-code/plotting-code/\
#     plot-nb-permutation-signflip.R
#
# Optional positional arguments:
#   1. simulation results RDS
#   2. parameter-grid R source (or a legacy RDS)
#
# By default, the script reads
#   results/simulation-results/
#     nb-permutation-signflip-continuous-dispersion-estimated-theta/
#     nb-permutation-signflip_results.rds
# and writes CSV tables and diagnostic figures beside that file. The two
# manuscript-focused figures are written directly to figures-and-tables.
# Figure conventions adapt the simulation summaries in
# https://github.com/Katsevich-Lab/spacrt-manuscript.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

sim_name <- "nb-permutation-signflip"
results_namespace <-
  "nb-permutation-signflip-continuous-dispersion-estimated-theta"
expected_schema <- "continuous_nb_dispersion_estimated_theta"

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
manuscript_directory <- file.path(
  project_root,
  "manuscript",
  "figures-and-tables"
)

results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}

resolve_path <- function(path) {
  if (grepl("^/", path)) {
    path
  } else {
    file.path(project_root, path)
  }
}

arguments <- commandArgs(trailingOnly = TRUE)
results_file <- if (length(arguments) >= 1L) {
  resolve_path(arguments[[1L]])
} else {
  file.path(
    results_root, "simulation-results", results_namespace,
    paste0(sim_name, "_results.rds")
  )
}
parameter_grid_file <- if (length(arguments) >= 2L) {
  resolve_path(arguments[[2L]])
} else {
  file.path(
    project_root,
    "reproduction-code", "simulation-code", "sim-spec-formal",
    paste0(sim_name, "-grid.R")
  )
}

if (!file.exists(results_file)) {
  stop(
    "Missing simulation results: ", results_file,
    "\nRun the ", sim_name, " simulation before plotting."
  )
}
if (!file.exists(parameter_grid_file)) {
  stop(
    "Missing parameter grid: ", parameter_grid_file,
    "\nRestore the simulation specifier sources before plotting."
  )
}

results_directory <- dirname(normalizePath(results_file))

required_grid_columns <- c(
  "grid_id", "n", "null_design", "scenario_id", "nuisance_mode",
  "dgp_schema",
  "z_distribution",
  "x_prevalence", "x_intercept", "x_z_slope",
  "mu0", "beta_z", "nb_dispersion", "nb_size",
  "variance_model", "working_family",
  "working_dispersion_mode", "working_nb_size",
  "nb_size_intercept", "nb_size_z_slope",
  "gamma_signal", "fold_change", "B_resample",
  "alpha", "ground_truth"
)
required_output_columns <- c(
  "p.value", "statistic", "success", "failure_reason", "degenerate",
  "theta_hat", "n_exposed", "n_exposed_nonzero",
  "score_eff_n", "max_score_share", "implementation_schema"
)
method_levels <- c(
  "parametric_score",
  "information_studentized_permutation",
  "studentized_permutation",
  "hgf_signflip_mc",
  "standardized_signflip"
)
method_labels <- c(
  "parametric_score" = "Parametric score",
  "information_studentized_permutation" =
    "Barry et al.",
  "studentized_permutation" = "Empirical-studentized permutation",
  "hgf_signflip_mc" = "Hemerik et al.",
  "standardized_signflip" = "De Santis et al."
)
# Preserve the colours of the three shared methods in combine.pdf. The two
# sign-flip procedures use the remaining green/cyan manuscript colours.
method_colors <- c(
  "Parametric score" = "#F8766D",
  "Barry et al." = "#C77CFF",
  "Empirical-studentized permutation" = "#619CFF",
  "Hemerik et al." = "#00BA38",
  "De Santis et al." = "#00BFC4"
)

scenario_levels <- c(
  "common_nb",
  "rare_nb",
  "severe_variance"
)
scenario_labels <- c(
  "common_nb" = "Common X\nCorrect NB variance",
  "rare_nb" = "Rare X\nCorrect NB variance",
  "severe_variance" = "Moderately sparse X\nCovariate-dependent NB dispersion"
)
focused_scenario_labels <- c(
  "common_nb" = "High SNR",
  "rare_nb" = "Low SNR",
  "severe_variance" = "Misspecified NB variance"
)

theme_spacrt <- function(base_size = 11, remove_major_grid = FALSE) {
  output <- theme_bw(base_size = base_size) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = base_size),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      strip.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = base_size + 1),
      plot.tag.position = c(0, 1)
    )
  if (remove_major_grid) {
    output <- output + theme(
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  }
  output
}

as_logical_strict <- function(x, column_name) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x)) {
    output <- rep(NA, length(x))
    output[!is.na(x)] <- x[!is.na(x)] != 0
    return(output)
  }
  value <- tolower(trimws(as.character(x)))
  output <- rep(NA, length(value))
  output[value %in% c("true", "t", "1", "yes", "y", "alternative", "alt")] <- TRUE
  output[value %in% c("false", "f", "0", "no", "n", "null", "h0")] <- FALSE
  bad <- !is.na(value) & is.na(output)
  if (any(bad)) {
    stop(
      "Could not interpret ", column_name, " as logical: ",
      paste(unique(value[bad]), collapse = ", ")
    )
  }
  output
}

wilson_interval <- function(successes, trials, conf_level = 0.95) {
  successes <- as.numeric(successes)
  trials <- as.numeric(trials)
  lower <- upper <- rep(NA_real_, length(trials))
  valid <- is.finite(successes) & is.finite(trials) &
    trials > 0 & successes >= 0 & successes <= trials
  if (any(valid)) {
    z <- qnorm(1 - (1 - conf_level) / 2)
    p <- successes[valid] / trials[valid]
    denominator <- 1 + z^2 / trials[valid]
    center <- (p + z^2 / (2 * trials[valid])) / denominator
    half_width <- z * sqrt(
      p * (1 - p) / trials[valid] + z^2 / (4 * trials[valid]^2)
    ) / denominator
    lower[valid] <- pmax(0, center - half_width)
    upper[valid] <- pmin(1, center + half_width)
  }
  data.frame(lower = lower, upper = upper)
}

add_wilson_columns <- function(data, count_column, total_column, prefix) {
  interval <- wilson_interval(
    data[[count_column]],
    data[[total_column]]
  )
  data[[paste0(prefix, "_lower")]] <- interval$lower
  data[[paste0(prefix, "_upper")]] <- interval$upper
  data
}

expand_nominal_alpha <- function(data, levels) {
  bind_rows(lapply(levels, function(level) {
    mutate(data, nominal_alpha = level)
  }))
}

format_number <- function(x) {
  trimws(formatC(x, format = "fg", digits = 4))
}

percent_label <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f%%", 100 * x))
}

signal_label <- function(x) {
  output <- sprintf("%.5f", x)
  output <- sub("0+$", "", output)
  output <- sub("\\.$", "", output)
  output[abs(x) < 1e-12] <- "0"
  sub("^0\\.", ".", output)
}

null_design_labels <- c(
  "MI" = "MI",
  "CI" = "CI"
)

add_grid_plot_labels <- function(data) {
  prevalence_values <- sort(unique(data$x_prevalence))
  prevalence_levels <- paste0("pX = ", format_number(prevalence_values))
  dispersion_pairs <- unique(
    data[is.finite(data$nb_dispersion) & is.finite(data$nb_size),
      c("nb_dispersion", "nb_size"), drop = FALSE]
  )
  dispersion_pairs <- dispersion_pairs[
    order(dispersion_pairs$nb_dispersion),
    , drop = FALSE
  ]
  dispersion_levels <- paste0(
    "NB dispersion alpha = ", format_number(dispersion_pairs$nb_dispersion),
    " (theta = ", format_number(dispersion_pairs$nb_size), ")"
  )
  if (any(data$variance_model == "covariate_dependent_nb")) {
    dispersion_levels <- c(
      dispersion_levels,
      "Covariate-dependent truth; estimated common theta"
    )
  }

  data |>
    mutate(
      signal_strength = abs(.data$gamma_signal),
      null_design_label = factor(
        unname(null_design_labels[.data$null_design]),
        levels = unname(null_design_labels[c("MI", "CI")])
      ),
      prevalence_label = factor(
        paste0("pX = ", format_number(.data$x_prevalence)),
        levels = prevalence_levels
      ),
      dispersion_label = factor(
        if_else(
          is.finite(.data$nb_dispersion) & is.finite(.data$nb_size),
          paste0(
            "NB dispersion alpha = ", format_number(.data$nb_dispersion),
            " (theta = ", format_number(.data$nb_size), ")"
          ),
          "Covariate-dependent truth; estimated common theta"
        ),
        levels = dispersion_levels
      ),
      scenario_label = factor(
        unname(scenario_labels[.data$scenario_id]),
        levels = unname(scenario_labels[scenario_levels])
      )
    )
}

save_plot_pair <- function(
    plot, stem, width, height, output_directory = results_directory,
    save_png = TRUE) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_directory, paste0(sim_name, "-", stem, ".pdf"))
  png_file <- file.path(output_directory, paste0(sim_name, "-", stem, ".png"))
  # Use a compressed PDF 1.5 device so that figures embed cleanly in the
  # manuscript's pdfLaTeX output without a PDF-version warning or conversion.
  pdf_device <- function(filename, ...) {
    grDevices::pdf(
      file = filename,
      version = "1.5",
      useDingbats = FALSE,
      compress = TRUE,
      ...
    )
  }
  ggsave(
    pdf_file, plot,
    width = width, height = height, units = "in",
    device = pdf_device, limitsize = FALSE
  )
  if (save_png) {
    ggsave(
      png_file, plot,
      width = width, height = height, units = "in",
      dpi = 300, bg = "white", limitsize = FALSE
    )
  }
  message("Saved ", normalizePath(pdf_file))
  if (save_png) message("Saved ", normalizePath(png_file))
}

write_summary <- function(data, stem) {
  output_file <- file.path(
    results_directory,
    paste0(sim_name, "-", stem, ".csv")
  )
  write.csv(data, output_file, row.names = FALSE, na = "")
  message("Saved ", normalizePath(output_file))
}

###############################################################################
# Load and normalize the simulatr output.
###############################################################################

results_object <- readRDS(results_file)
results <- if (
  is.list(results_object) &&
    !is.data.frame(results_object) &&
    !is.null(results_object$results)
) {
  results_object$results
} else {
  results_object
}
if (!is.data.frame(results)) {
  stop("The results RDS must be a data frame or contain a $results data frame.")
}
if (!all(c("grid_id", "run_id", "method", "output") %in% names(results))) {
  stop("Results must contain grid_id, run_id, method, and output columns.")
}
if (!is.list(results$output)) {
  stop("The results output column must be a list-column.")
}

if (grepl("[.]rds$", parameter_grid_file, ignore.case = TRUE)) {
  parameter_grid <- as.data.frame(readRDS(parameter_grid_file))
} else {
  grid_environment <- new.env(parent = globalenv())
  sys.source(parameter_grid_file, envir = grid_environment)
  if (!exists("parameter_grid", envir = grid_environment, inherits = FALSE)) {
    stop("Parameter-grid source did not create parameter_grid: ", parameter_grid_file)
  }
  parameter_grid <- as.data.frame(grid_environment$parameter_grid)
}
missing_grid_columns <- setdiff(required_grid_columns, names(parameter_grid))
if (length(missing_grid_columns) > 0L) {
  stop(
    "Parameter grid is missing required columns: ",
    paste(missing_grid_columns, collapse = ", ")
  )
}
if (anyDuplicated(parameter_grid$grid_id)) {
  stop("Parameter-grid grid_id values must be unique.")
}
if (!identical(
  unique(as.character(parameter_grid$dgp_schema)),
  expected_schema
)) {
  stop(
    "Parameter grid must use dgp_schema = ", expected_schema, "; ",
    "legacy grids are not accepted."
  )
}
parameter_grid$scenario_id <- as.character(parameter_grid$scenario_id)
if (anyNA(parameter_grid$scenario_id) ||
    any(!nzchar(parameter_grid$scenario_id))) {
  stop("Parameter-grid scenario_id values may not be missing or empty.")
}
unexpected_scenarios <- setdiff(
  unique(parameter_grid$scenario_id),
  scenario_levels
)
missing_scenarios <- setdiff(
  scenario_levels,
  unique(parameter_grid$scenario_id)
)
if (length(unexpected_scenarios) > 0L || length(missing_scenarios) > 0L) {
  stop(
    "Expected scenario_id values ", paste(scenario_levels, collapse = ", "),
    "; unexpected: ",
    if (length(unexpected_scenarios) > 0L) {
      paste(unexpected_scenarios, collapse = ", ")
    } else {
      "none"
    },
    "; missing: ",
    if (length(missing_scenarios) > 0L) {
      paste(missing_scenarios, collapse = ", ")
    } else {
      "none"
    },
    "."
  )
}
if (!setequal(unique(parameter_grid$null_design), c("MI", "CI"))) {
  stop("Parameter grid must contain both MI and CI null designs.")
}
parameter_grid$ground_truth <- as_logical_strict(
  parameter_grid$ground_truth,
  "ground_truth"
)
if (anyNA(parameter_grid$ground_truth) ||
    any(!is.finite(parameter_grid$gamma_signal)) ||
    any(parameter_grid$gamma_signal > 0)) {
  stop("The focused grid requires nonmissing truth and nonpositive signals.")
}
grid_signal_contract <- parameter_grid |>
  group_by(.data$null_design, .data$scenario_id) |>
  summarise(
    n_signal_rows = n(),
    n_null = sum(!.data$ground_truth),
    n_alternative = sum(.data$ground_truth),
    n_signal_values = n_distinct(.data$gamma_signal),
    includes_zero = any(abs(.data$gamma_signal) < 1e-12),
    equal_spacing_error = {
      magnitude <- sort(unique(abs(.data$gamma_signal)))
      if (length(magnitude) == 5L) {
        max(abs(diff(magnitude) - mean(diff(magnitude))))
      } else {
        Inf
      }
    },
    .groups = "drop"
  )
bad_signal_contract <- grid_signal_contract |>
  filter(
    .data$n_signal_rows != 5L |
      .data$n_null != 1L |
      .data$n_alternative != 4L |
      .data$n_signal_values != 5L |
      !.data$includes_zero |
      .data$equal_spacing_error > 1e-10
  )
if (nrow(grid_signal_contract) != 6L || nrow(bad_signal_contract) > 0L) {
  stop(
    "Every MI/CI scenario must have five equally spaced signal values: ",
    "one null at zero and four alternatives."
  )
}

results <- results |>
  mutate(
    grid_id = suppressWarnings(as.integer(as.character(.data$grid_id))),
    run_id = suppressWarnings(as.integer(as.character(.data$run_id))),
    method = as.character(.data$method)
  )
if (anyNA(results$grid_id) || anyNA(results$run_id) || anyNA(results$method)) {
  stop("grid_id, run_id, and method may not be missing.")
}
unknown_grid_ids <- setdiff(unique(results$grid_id), parameter_grid$grid_id)
if (length(unknown_grid_ids) > 0L) {
  stop(
    "Results contain grid_id values absent from the parameter grid: ",
    paste(unknown_grid_ids, collapse = ", ")
  )
}
unexpected_methods <- setdiff(unique(results$method), method_levels)
if (length(unexpected_methods) > 0L) {
  stop("Unexpected methods: ", paste(unexpected_methods, collapse = ", "))
}
if (anyDuplicated(results[c("grid_id", "run_id", "method")])) {
  stop("Results contain duplicate grid_id/run_id/method rows.")
}

expected_groups <- tidyr::expand_grid(
  grid_id = parameter_grid$grid_id,
  method = method_levels
)
completion_summary <- results |>
  count(.data$grid_id, .data$method, name = "n_repetitions") |>
  right_join(expected_groups, by = c("grid_id", "method")) |>
  arrange(.data$grid_id, factor(.data$method, levels = method_levels))
missing_groups <- completion_summary |>
  filter(is.na(.data$n_repetitions))
if (nrow(missing_groups) > 0L) {
  stop(
    "Results are incomplete: missing ", nrow(missing_groups),
    " grid/method combinations."
  )
}

run_id_summary <- results |>
  group_by(.data$grid_id, .data$method) |>
  summarise(
    n_repetitions = n(),
    minimum_run_id = min(.data$run_id),
    maximum_run_id = max(.data$run_id),
    n_distinct_run_id = n_distinct(.data$run_id),
    .groups = "drop"
  )
bad_run_sequences <- run_id_summary |>
  filter(
    .data$minimum_run_id != 1L |
      .data$n_distinct_run_id != .data$n_repetitions |
      .data$maximum_run_id != .data$n_repetitions
  )
if (nrow(bad_run_sequences) > 0L) {
  stop(
    "Results are incomplete: run_id must be exactly 1, ..., B within every ",
    "grid/method combination."
  )
}
repetition_counts <- sort(unique(run_id_summary$n_repetitions))
if (length(repetition_counts) != 1L) {
  stop(
    "Results are incomplete: repetition counts differ across grid/method ",
    "combinations (", paste(repetition_counts, collapse = ", "), ")."
  )
}
if (repetition_counts != 2000L) {
  stop(
    "Canonical manuscript figures require exactly 2,000 repetitions per ",
    "grid/method combination; found ", repetition_counts, "."
  )
}
message(
  "Validated ", nrow(parameter_grid), " grid cells x ",
  length(method_levels), " methods x ", repetition_counts,
  " repetitions (", nrow(results), " rows)."
)

columns_to_join <- setdiff(required_grid_columns[-1L], names(results))
results <- results |>
  left_join(
    parameter_grid[, c("grid_id", columns_to_join), drop = FALSE],
    by = "grid_id"
  ) |>
  mutate(
    method = as.character(.data$method),
    null_design = as.character(.data$null_design),
    scenario_id = as.character(.data$scenario_id),
    nuisance_mode = as.character(.data$nuisance_mode),
    dgp_schema = as.character(.data$dgp_schema),
    z_distribution = as.character(.data$z_distribution),
    variance_model = as.character(.data$variance_model),
    working_family = as.character(.data$working_family),
    working_dispersion_mode = as.character(.data$working_dispersion_mode)
  ) |>
  unnest_wider("output", strict = FALSE, names_repair = "check_unique")

missing_result_columns <- setdiff(
  c(required_grid_columns, required_output_columns),
  names(results)
)
if (length(missing_result_columns) > 0L) {
  stop(
    "Joined results are missing required columns: ",
    paste(missing_result_columns, collapse = ", ")
  )
}

missing_methods <- setdiff(method_levels, unique(results$method))
if (length(missing_methods) > 0L) {
  warning("Results do not yet contain: ", paste(missing_methods, collapse = ", "))
}

results <- results |>
  mutate(
    implementation_schema = as.character(.data$implementation_schema),
    ground_truth = as_logical_strict(.data$ground_truth, "ground_truth"),
    success_flag = as_logical_strict(.data$success, "success"),
    degenerate_flag = as_logical_strict(.data$degenerate, "degenerate"),
    p_value = suppressWarnings(as.numeric(.data[["p.value"]])),
    statistic = suppressWarnings(as.numeric(.data$statistic)),
    theta_hat = suppressWarnings(as.numeric(.data$theta_hat)),
    n_exposed = suppressWarnings(as.numeric(.data$n_exposed)),
    n_exposed_nonzero = suppressWarnings(as.numeric(.data$n_exposed_nonzero)),
    score_eff_n = suppressWarnings(as.numeric(.data$score_eff_n)),
    max_score_share = suppressWarnings(as.numeric(.data$max_score_share)),
    failure_reason = as.character(.data$failure_reason),
    valid_p = is.finite(.data$p_value) &
      .data$p_value >= 0 & .data$p_value <= 1,
    reported_failure = is.na(.data$success_flag) | !.data$success_flag,
    invalid_p = !.data$valid_p,
    analysis_valid = !.data$reported_failure & .data$valid_p,
    analysis_failure = !.data$analysis_valid,
    degenerate_flag = !is.na(.data$degenerate_flag) & .data$degenerate_flag,
    method_label = factor(
      unname(method_labels[.data$method]),
      levels = unname(method_labels[method_levels])
    )
  )

if (any(is.na(results$ground_truth))) {
  stop("ground_truth may not be missing.")
}
if (anyNA(results$dgp_schema) || anyNA(results$implementation_schema) ||
    any(results$dgp_schema != expected_schema) ||
    any(results$implementation_schema != expected_schema) ||
    any(results$implementation_schema != results$dgp_schema)) {
  stop(
    "Results must carry matching ", expected_schema, " DGP and ",
    "implementation schemas; legacy outputs are not accepted."
  )
}
unexpected_nuisance_modes <- setdiff(unique(results$nuisance_mode), "fitted")
if (length(unexpected_nuisance_modes) > 0L) {
  stop(
    "This plotter expects the expanded fitted-nuisance grid; found: ",
    paste(unexpected_nuisance_modes, collapse = ", ")
  )
}
unexpected_null_designs <- setdiff(unique(results$null_design), c("MI", "CI"))
if (length(unexpected_null_designs) > 0L) {
  stop("Unexpected null designs: ", paste(unexpected_null_designs, collapse = ", "))
}
if (any(!is.finite(results$x_prevalence) | results$x_prevalence <= 0)) {
  stop("x_prevalence must be finite and positive.")
}
severe_variance_row <- results$scenario_id == "severe_variance"
ordinary_nb_row <- !severe_variance_row
if (any(results$z_distribution != "normal")) {
  stop("Every focused scenario must use Z distributed as N(0, 1).")
}
if (any(
  ordinary_nb_row &
    (!is.finite(results$nb_dispersion) | results$nb_dispersion <= 0 |
      !is.finite(results$nb_size) | results$nb_size <= 0)
)) {
  stop("Ordinary NB rows require finite positive nb_dispersion and nb_size.")
}
if (any(
  ordinary_nb_row &
    abs(results$nb_dispersion * results$nb_size - 1) > 1e-10
)) {
  stop("nb_size must equal 1 / nb_dispersion in ordinary NB rows.")
}
if (any(
  ordinary_nb_row &
    (results$variance_model != "constant_nb" |
      results$working_family != "negative_binomial" |
      results$working_dispersion_mode != "estimated" |
      !is.na(results$working_nb_size) |
      !is.na(results$nb_size_intercept) |
      !is.na(results$nb_size_z_slope))
)) {
  stop(
    "Ordinary NB rows require a constant-NB outcome and an estimated common ",
    "NB working dispersion."
  )
}
estimated_valid_row <-
  results$analysis_valid & results$working_dispersion_mode == "estimated"
if (any(
  estimated_valid_row &
    (!is.finite(results$theta_hat) | results$theta_hat <= 0)
)) {
  stop("Analysis-valid rows require a finite positive estimated NB size.")
}
if (any(
  severe_variance_row &
    (!is.na(results$nb_dispersion) | !is.na(results$nb_size))
)) {
  stop("nb_dispersion and nb_size must be NA in severe-variance rows.")
}
if (any(
  severe_variance_row &
    (!is.finite(results$x_prevalence) |
      abs(results$x_prevalence - 0.14) > 1e-12 |
      !is.finite(results$x_z_slope) |
      abs(
        results$x_z_slope -
          ifelse(results$null_design == "MI", 0, 0.5)
      ) > 1e-12 |
      !is.finite(results$x_intercept) |
      abs(
        results$x_intercept -
          ifelse(
            results$null_design == "MI",
            -1.815289966638249,
            -1.902689139395545
          )
      ) > 1e-10 |
      !is.finite(results$mu0) |
      abs(results$mu0 - 10) > 1e-12 |
      !is.finite(results$beta_z) |
      abs(results$beta_z - 0.5) > 1e-12 |
      results$variance_model != "covariate_dependent_nb" |
      results$working_family != "negative_binomial" |
      results$working_dispersion_mode != "estimated" |
      !is.na(results$working_nb_size) |
      !is.finite(results$nb_size_intercept) |
      abs(results$nb_size_intercept - 46) > 1e-12 |
      !is.finite(results$nb_size_z_slope) |
      abs(results$nb_size_z_slope + 2.1) > 1e-12)
)) {
  stop(
    "Severe-variance rows must use the canonical continuous-Z exposure and ",
    "mean models, theta(Z) = 46 exp(-2.1 Z), and a common NB working size ",
    "estimated by glm.nb."
  )
}
if (any(results$analysis_valid & !is.finite(results$statistic))) {
  warning("Some successful rows have non-finite statistics; see diagnostics.")
}

scenario_variables <- c(
  "method", "n", "null_design", "scenario_id", "nuisance_mode",
  "dgp_schema", "implementation_schema",
  "z_distribution",
  "x_prevalence", "x_intercept", "x_z_slope",
  "mu0", "beta_z", "nb_dispersion", "nb_size",
  "variance_model", "working_family",
  "working_dispersion_mode", "working_nb_size",
  "nb_size_intercept", "nb_size_z_slope",
  "gamma_signal", "fold_change", "B_resample",
  "alpha", "ground_truth"
)

summarize_rejections <- function(data, cutoff_column) {
  data |>
    mutate(
      .cutoff = .data[[cutoff_column]],
      .eligible = .data$analysis_valid & is.finite(.data$.cutoff),
      .reject = .data$.eligible & .data$p_value <= .data$.cutoff
    ) |>
    group_by(across(all_of(c(scenario_variables, "nominal_alpha")))) |>
    summarise(
      n_attempted = n(),
      n_reported_success = sum(!.data$reported_failure),
      n_valid = sum(.data$.eligible),
      rejections = sum(.data$.reject),
      critical_value = if (any(is.finite(.data$.cutoff))) {
        .data$.cutoff[which(is.finite(.data$.cutoff))[1L]]
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    mutate(
      estimate = if_else(
        .data$n_valid > 0,
        .data$rejections / .data$n_valid,
        NA_real_
      )
    ) |>
    add_wilson_columns("rejections", "n_valid", "mc")
}

###############################################################################
# Type-I error at 0.05 and 0.01, with 95% Wilson Monte Carlo intervals.
###############################################################################

typeI_alpha_levels <- c(0.05, 0.01)
typeI_summary <- results |>
  filter(!.data$ground_truth) |>
  expand_nominal_alpha(typeI_alpha_levels) |>
  summarize_rejections("nominal_alpha") |>
  mutate(
    method_label = unname(method_labels[.data$method]),
    mc_interval = "Wilson 95%"
  ) |>
  arrange(
    .data$nominal_alpha, .data$null_design, .data$scenario_id,
    .data$nuisance_mode,
    .data$n, .data$x_prevalence, .data$mu0,
    .data$nb_dispersion, .data$nb_size,
    factor(.data$method, levels = method_levels)
  )

###############################################################################
# Nominal-cutoff and size-adjusted rejection rates at the configured alpha.
# Under alternatives these rates are power; at gamma = 0 they are null
# rejection rates. Including the null makes the line plots start at signal 0.
# The size-adjusted critical value is the type-1 empirical alpha quantile of
# successful null p-values for the same method and DGP design.
###############################################################################

power_alpha_levels <- sort(unique(results$alpha[is.finite(results$alpha)]))
if (length(power_alpha_levels) == 0L) {
  power_alpha_levels <- 0.05
}

raw_power_summary <- results |>
  expand_nominal_alpha(power_alpha_levels) |>
  summarize_rejections("nominal_alpha") |>
  mutate(
    power_type = "raw_power",
    null_n_valid = NA_integer_,
    null_rejections_at_cutoff = NA_integer_,
    null_size_at_cutoff = NA_real_,
    null_size_lower = NA_real_,
    null_size_upper = NA_real_
  )

size_adjustment_keys <- c(
  "method", "n", "null_design", "scenario_id", "nuisance_mode",
  "dgp_schema", "implementation_schema",
  "z_distribution",
  "x_prevalence", "x_intercept", "x_z_slope",
  "mu0", "beta_z", "nb_dispersion", "nb_size",
  "variance_model", "working_family",
  "working_dispersion_mode", "working_nb_size",
  "nb_size_intercept", "nb_size_z_slope",
  "B_resample", "alpha",
  "nominal_alpha"
)

null_for_adjustment <- results |>
  filter(!.data$ground_truth) |>
  expand_nominal_alpha(power_alpha_levels)

size_adjustment_cutoffs <- null_for_adjustment |>
  group_by(across(all_of(size_adjustment_keys))) |>
  summarise(
    null_n_attempted = n(),
    null_n_valid = sum(.data$analysis_valid),
    size_adjusted_cutoff = if (sum(.data$analysis_valid) > 0L) {
      unname(quantile(
        .data$p_value[.data$analysis_valid],
        probs = first(.data$nominal_alpha),
        type = 1,
        na.rm = TRUE
      ))
    } else {
      NA_real_
    },
    .groups = "drop"
  )

null_cutoff_attainment <- null_for_adjustment |>
  left_join(size_adjustment_cutoffs, by = size_adjustment_keys) |>
  mutate(
    eligible_at_cutoff = .data$analysis_valid &
      is.finite(.data$size_adjusted_cutoff),
    reject_at_cutoff = .data$eligible_at_cutoff &
      .data$p_value <= .data$size_adjusted_cutoff
  ) |>
  group_by(across(all_of(c(size_adjustment_keys, "size_adjusted_cutoff")))) |>
  summarise(
    null_n_attempted = n(),
    null_n_valid = sum(.data$eligible_at_cutoff),
    null_rejections_at_cutoff = sum(.data$reject_at_cutoff),
    .groups = "drop"
  ) |>
  mutate(
    null_size_at_cutoff = if_else(
      .data$null_n_valid > 0,
      .data$null_rejections_at_cutoff / .data$null_n_valid,
      NA_real_
    )
  ) |>
  add_wilson_columns(
    "null_rejections_at_cutoff", "null_n_valid", "null_size"
  )

size_adjustment_cutoffs <- size_adjustment_cutoffs |>
  select(-all_of(c("null_n_attempted", "null_n_valid"))) |>
  left_join(
    null_cutoff_attainment,
    by = c(size_adjustment_keys, "size_adjusted_cutoff")
  ) |>
  mutate(
    method_label = unname(method_labels[.data$method]),
    quantile_type = 1L,
    mc_interval = "Wilson 95%"
  ) |>
  arrange(
    .data$nominal_alpha, .data$null_design, .data$scenario_id,
    .data$nuisance_mode,
    factor(.data$method, levels = method_levels)
  )

adjusted_power_data <- results |>
  expand_nominal_alpha(power_alpha_levels) |>
  left_join(
    select(
      size_adjustment_cutoffs,
      all_of(size_adjustment_keys),
      "size_adjusted_cutoff", "null_n_valid",
      "null_rejections_at_cutoff", "null_size_at_cutoff",
      "null_size_lower", "null_size_upper"
    ),
    by = size_adjustment_keys
  )

if (nrow(adjusted_power_data) > 0L &&
    any(!is.finite(adjusted_power_data$size_adjusted_cutoff))) {
  warning(
    "Some cells have no matching valid null cutoff; their size-adjusted ",
    "rejection rate is reported as NA."
  )
}

adjusted_power_summary <- adjusted_power_data |>
  summarize_rejections("size_adjusted_cutoff") |>
  left_join(
    distinct(
      adjusted_power_data,
      across(all_of(c(scenario_variables, "nominal_alpha"))),
      .data$null_n_valid, .data$null_rejections_at_cutoff,
      .data$null_size_at_cutoff, .data$null_size_lower,
      .data$null_size_upper
    ),
    by = c(scenario_variables, "nominal_alpha")
  ) |>
  mutate(power_type = "size_adjusted_power")

power_summary <- bind_rows(raw_power_summary, adjusted_power_summary) |>
  mutate(
    method_label = unname(method_labels[.data$method]),
    mc_interval = "Wilson 95%"
  ) |>
  arrange(
    .data$power_type, .data$nominal_alpha,
    .data$null_design, .data$scenario_id, .data$nuisance_mode,
    .data$n, .data$x_prevalence, .data$mu0,
    .data$nb_dispersion, .data$nb_size,
    .data$fold_change, factor(.data$method, levels = method_levels)
  )

###############################################################################
# Failures and diagnostics.
###############################################################################

failure_summary <- results |>
  group_by(across(all_of(scenario_variables))) |>
  summarise(
    n_attempted = n(),
    n_reported_success = sum(!.data$reported_failure),
    n_reported_failure = sum(.data$reported_failure),
    n_analysis_valid = sum(.data$analysis_valid),
    n_analysis_failure = sum(.data$analysis_failure),
    n_degenerate = sum(.data$degenerate_flag),
    n_invalid_p = sum(.data$invalid_p),
    .groups = "drop"
  ) |>
  mutate(
    reported_failure_rate = .data$n_reported_failure / .data$n_attempted,
    analysis_failure_rate = .data$n_analysis_failure / .data$n_attempted,
    degenerate_rate = .data$n_degenerate / .data$n_attempted,
    invalid_p_rate = .data$n_invalid_p / .data$n_attempted
  ) |>
  add_wilson_columns(
    "n_reported_failure", "n_attempted", "reported_failure"
  ) |>
  add_wilson_columns(
    "n_analysis_failure", "n_attempted", "analysis_failure"
  ) |>
  add_wilson_columns("n_degenerate", "n_attempted", "degenerate") |>
  add_wilson_columns("n_invalid_p", "n_attempted", "invalid_p") |>
  mutate(
    method_label = unname(method_labels[.data$method]),
    mc_interval = "Wilson 95%"
  ) |>
  arrange(
    .data$ground_truth, .data$null_design, .data$scenario_id,
    .data$nuisance_mode,
    .data$fold_change, factor(.data$method, levels = method_levels)
  )

failure_denominators <- failure_summary |>
  select(all_of(scenario_variables), "n_attempted")

failure_reason_summary <- results |>
  filter(.data$analysis_failure) |>
  mutate(
    failure_reason = case_when(
      .data$invalid_p & !.data$reported_failure ~ "invalid_or_missing_p_value",
      is.na(.data$failure_reason) | trimws(.data$failure_reason) == "" ~
        "unspecified_failure",
      TRUE ~ .data$failure_reason
    )
  ) |>
  count(
    across(all_of(scenario_variables)),
    .data$failure_reason,
    name = "failures_with_reason"
  ) |>
  left_join(failure_denominators, by = scenario_variables) |>
  mutate(
    fraction_of_attempts = .data$failures_with_reason / .data$n_attempted,
    method_label = unname(method_labels[.data$method])
  ) |>
  arrange(
    .data$null_design, .data$scenario_id, .data$nuisance_mode,
    .data$fold_change,
    factor(.data$method, levels = method_levels),
    desc(.data$failures_with_reason)
  )

diagnostic_columns <- c(
  "statistic", "theta_hat", "n_exposed", "n_exposed_nonzero",
  "score_eff_n", "max_score_share"
)
diagnostic_summary <- results |>
  select(
    all_of(scenario_variables), "analysis_valid",
    all_of(diagnostic_columns)
  ) |>
  pivot_longer(
    cols = all_of(diagnostic_columns),
    names_to = "diagnostic",
    values_to = "value"
  ) |>
  group_by(across(all_of(c(scenario_variables, "diagnostic")))) |>
  summarise(
    n_attempted = n(),
    n_available = sum(is.finite(.data$value)),
    n_available_success = sum(.data$analysis_valid & is.finite(.data$value)),
    mean = if (any(is.finite(.data$value))) {
      mean(.data$value[is.finite(.data$value)])
    } else {
      NA_real_
    },
    sd = if (sum(is.finite(.data$value)) > 1L) {
      sd(.data$value[is.finite(.data$value)])
    } else {
      NA_real_
    },
    minimum = if (any(is.finite(.data$value))) {
      min(.data$value[is.finite(.data$value)])
    } else {
      NA_real_
    },
    q05 = if (any(is.finite(.data$value))) {
      unname(quantile(.data$value, 0.05, na.rm = TRUE))
    } else {
      NA_real_
    },
    q25 = if (any(is.finite(.data$value))) {
      unname(quantile(.data$value, 0.25, na.rm = TRUE))
    } else {
      NA_real_
    },
    median = if (any(is.finite(.data$value))) {
      median(.data$value, na.rm = TRUE)
    } else {
      NA_real_
    },
    q75 = if (any(is.finite(.data$value))) {
      unname(quantile(.data$value, 0.75, na.rm = TRUE))
    } else {
      NA_real_
    },
    q95 = if (any(is.finite(.data$value))) {
      unname(quantile(.data$value, 0.95, na.rm = TRUE))
    } else {
      NA_real_
    },
    maximum = if (any(is.finite(.data$value))) {
      max(.data$value[is.finite(.data$value)])
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  mutate(method_label = unname(method_labels[.data$method])) |>
  arrange(
    .data$diagnostic, .data$null_design, .data$scenario_id,
    .data$nuisance_mode,
    .data$fold_change, factor(.data$method, levels = method_levels)
  )

###############################################################################
# Save tables.
###############################################################################

write_summary(typeI_summary, "typeI-summary")
write_summary(power_summary, "power-summary")
write_summary(size_adjustment_cutoffs, "size-adjustment-cutoffs")
write_summary(failure_summary, "failure-summary")
write_summary(failure_reason_summary, "failure-reasons")
write_summary(diagnostic_summary, "diagnostic-summary")
write_summary(completion_summary, "completion-summary")

###############################################################################
# Type-I figure.
###############################################################################

if (nrow(typeI_summary) > 0L) {
  typeI_plot_data <- typeI_summary |>
    mutate(
      method_label = factor(
        .data$method_label,
        levels = unname(method_labels[method_levels])
      ),
      alpha_label = factor(
        paste0("alpha = ", format_number(.data$nominal_alpha)),
        levels = paste0(
          "alpha = ", format_number(sort(unique(.data$nominal_alpha), decreasing = TRUE))
        )
      )
    ) |>
    add_grid_plot_labels()

  nominal_lines <- typeI_plot_data |>
    distinct(.data$alpha_label, .data$nominal_alpha)

  typeI_line_data <- typeI_plot_data |>
    group_by(
      .data$alpha_label, .data$null_design_label,
      .data$scenario_label, .data$method_label
    ) |>
    filter(n_distinct(.data$x_prevalence) > 1L) |>
    ungroup()

  typeI_plot <- ggplot(
    typeI_plot_data,
    aes(
      x = .data$x_prevalence,
      y = .data$estimate,
      colour = .data$method_label,
      group = .data$method_label
    )
  ) +
    geom_hline(
      data = nominal_lines,
      aes(yintercept = .data$nominal_alpha),
      inherit.aes = FALSE,
      linetype = "dashed", colour = "black", linewidth = 0.55
    ) +
    geom_errorbar(
      aes(ymin = .data$mc_lower, ymax = .data$mc_upper),
      width = 0, linewidth = 0.55
    ) +
    geom_line(data = typeI_line_data, linewidth = 0.75) +
    geom_point(size = 2.7) +
    facet_grid(
      rows = vars(.data$alpha_label, .data$null_design_label),
      cols = vars(.data$scenario_label),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_x_log10(
      breaks = sort(unique(typeI_plot_data$x_prevalence)),
      labels = format_number
    ) +
    scale_colour_manual(values = method_colors, name = NULL) +
    scale_y_continuous(
      labels = percent_label,
      expand = expansion(mult = c(0.03, 0.12))
    ) +
    labs(
      x = "Marginal exposure prevalence",
      y = "Empirical Type-I error (95% MC interval)",
      title = "Score-test calibration under the null",
      caption = paste0(
        "Dashed lines are the nominal levels; intervals are Wilson 95% ",
        "Monte Carlo intervals. MI = marginal independence; ",
        "CI = conditional independence."
      )
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_spacrt(base_size = 11) +
    theme(
      plot.caption = element_text(hjust = 0)
    )

  save_plot_pair(
    typeI_plot, "typeI",
    width = 13,
    height = 8
  )
}

###############################################################################
# Nominal-cutoff and size-adjusted rejection-rate figure.
###############################################################################

if (nrow(power_summary) > 0L && any(is.finite(power_summary$estimate))) {
  power_plot_data <- power_summary |>
    mutate(
      method_label = factor(
        .data$method_label,
        levels = unname(method_labels[method_levels])
      ),
      power_label = factor(
        recode(
          .data$power_type,
          raw_power = "Nominal-cutoff rejection",
          size_adjusted_power = "Size-adjusted rejection"
        ),
        levels = c("Nominal-cutoff rejection", "Size-adjusted rejection")
      ),
      alpha_label = paste0("alpha = ", format_number(.data$nominal_alpha))
    ) |>
    add_grid_plot_labels()

  power_nominal_lines <- power_plot_data |>
    distinct(
      .data$power_label, .data$null_design_label, .data$alpha_label,
      .data$scenario_label, .data$nominal_alpha
    )

  power_plot <- ggplot(
    power_plot_data,
    aes(
      x = .data$signal_strength,
      y = .data$estimate,
      colour = .data$method_label,
      group = .data$method_label
    )
  ) +
    geom_hline(
      data = power_nominal_lines,
      aes(yintercept = .data$nominal_alpha),
      inherit.aes = FALSE,
      linetype = "dashed", colour = "black", linewidth = 0.55
    ) +
    geom_errorbar(
      aes(ymin = .data$mc_lower, ymax = .data$mc_upper),
      width = 0, linewidth = 0.5
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2.5) +
    facet_grid(
      rows = vars(
        .data$power_label, .data$null_design_label, .data$alpha_label
      ),
      cols = vars(.data$scenario_label),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_colour_manual(values = method_colors, name = "Calibration") +
    scale_x_continuous(
      breaks = sort(unique(power_plot_data$signal_strength)),
      labels = signal_label
    ) +
    scale_y_continuous(
      limits = c(0, 1), labels = percent_label,
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    labs(
      x = "Signal magnitude |gamma| (log-fold decrease)",
      y = "Rejection rate (95% MC interval)",
      title = "Nominal-cutoff and size-adjusted rejection rates",
      caption = paste0(
        "Size adjustment uses the type-1 empirical null alpha-quantile ",
        "for the same method and design. MI = marginal independence; ",
        "CI = conditional independence."
      )
    ) +
    guides(
      colour = guide_legend(order = 1, nrow = 1, byrow = TRUE)
    ) +
    theme_spacrt(base_size = 11) +
    theme(plot.caption = element_text(hjust = 0))

  n_power_rows <- n_distinct(interaction(
    power_plot_data$power_label,
    power_plot_data$null_design_label,
    power_plot_data$alpha_label,
    drop = TRUE
  ))
  save_plot_pair(
    power_plot, "power",
    width = 13,
    height = max(8, 2 * n_power_rows)
  )
}

###############################################################################
# Manuscript-style summary: Type-I error above rejection-rate curves.
###############################################################################

if (exists("typeI_plot", inherits = FALSE) &&
    exists("power_plot", inherits = FALSE)) {
  manuscript_summary_plot <- (
    (typeI_plot +
      labs(
        tag = "a",
        title = "Type-I error",
        caption = NULL
      ) +
      theme(legend.position = "none")) /
      (power_plot +
        labs(
          tag = "b",
          title = "Nominal-cutoff and size-adjusted rejection rates",
          caption = NULL
        ) +
        theme(legend.position = "bottom", legend.box = "horizontal")) +
      plot_layout(heights = c(1, 1.05)) +
      plot_annotation(
        caption = paste0(
          "MI = marginal independence; CI = conditional independence. ",
          "Intervals are Wilson 95% Monte Carlo intervals."
        )
      )
  )

  save_plot_pair(
    manuscript_summary_plot, "manuscript-summary",
    width = 13,
    height = 14
  )
}

###############################################################################
# Focused 3 x 3 manuscript figures, one each for MI and CI.
#
# Rows contrast common X, rare X, and severe variance misspecification. Columns
# show null-p-value calibration, Type-I error, and nominal-cutoff rejection
# rates. The rejection-rate plots include gamma = 0, so their first point is
# the corresponding null rejection rate.
###############################################################################

focused_scenarios <- scenario_levels
focused_alpha <- 0.05
focused_power_types <- c("raw_power", "size_adjusted_power")

focused_typeI <- typeI_summary |>
  filter(
    .data$nominal_alpha == focused_alpha,
    .data$scenario_id %in% focused_scenarios
  )
focused_power <- power_summary |>
  filter(
    .data$power_type %in% focused_power_types,
    .data$nominal_alpha == focused_alpha,
    .data$scenario_id %in% focused_scenarios
  ) |>
  mutate(signal_strength = abs(.data$gamma_signal))
focused_null_pvalues <- results |>
  filter(
    !.data$ground_truth,
    .data$analysis_valid,
    .data$scenario_id %in% focused_scenarios
  )

if (nrow(focused_typeI) > 0L &&
    nrow(focused_power) > 0L &&
    nrow(focused_null_pvalues) > 0L) {
  focused_method_levels <- unname(method_labels[method_levels])
  qq_resampling_floor <- 1 / (max(focused_null_pvalues$B_resample) + 1)
  # The analytic parametric p-value can be far below the Monte Carlo floor.
  # Keep the scientifically relevant tail legible and clamp only its plotted
  # position, never the value used in a rejection-rate calculation.
  qq_panel_floor <- min(qq_resampling_floor / 10^0.12, 1e-4)
  # Leave a small unlabeled buffer beyond the smallest displayed QQ tick.
  # Points clamped to the plotting floor are then fully visible rather than
  # being bisected by the panel boundary.
  qq_display_floor <- qq_panel_floor / 10^0.15
  qq_exponents <- seq(0, floor(log10(qq_panel_floor)), by = -1)
  qq_breaks <- 10^qq_exponents
  qq_labels <- as.expression(lapply(qq_exponents, function(exponent) {
    if (exponent == 0) 1 else bquote(10^.(exponent))
  }))
  typeI_upper <- max(
    0.15,
    ceiling(100 * max(focused_typeI$mc_upper, na.rm = TRUE) * 1.05) / 100
  )

  focused_text_panel <- function(
      label, font_size = 11.5, font_face = "bold", rotation = 0) {
    patchwork::wrap_elements(full = grid::textGrob(
      label,
      rot = rotation,
      gp = grid::gpar(fontsize = font_size, fontface = font_face)
    ))
  }

  make_focused_qq <- function(null_design_value, scenario_value, show_x) {
    qq_data <- focused_null_pvalues |>
      filter(
        .data$null_design == null_design_value,
        .data$scenario_id == scenario_value
      ) |>
      mutate(
        method_label = factor(
          unname(method_labels[.data$method]),
          levels = focused_method_levels
        )
      ) |>
      group_by(.data$method, .data$method_label) |>
      arrange(.data$p_value, .by_group = TRUE) |>
      mutate(
        qq_rank = row_number(),
        qq_n = n(),
        expected_p = (.data$qq_rank - 0.5) / .data$qq_n,
        observed_p_plot = pmax(.data$p_value, qq_panel_floor)
      ) |>
      ungroup()

    qq_group_sizes <- qq_data |>
      distinct(.data$method, .data$qq_n)
    if (nrow(qq_group_sizes) == 0L) {
      stop("Focused QQ panel has no valid null p-values.")
    }
    # If method-specific numerical failures produce unequal valid counts, use
    # the smallest group for one conservative common pointwise envelope.
    qq_n <- min(qq_group_sizes$qq_n)
    qq_band_data <- tibble(
      qq_rank = seq_len(qq_n),
      expected_p = (.data$qq_rank - 0.5) / qq_n,
      p_lower = qbeta(0.025, .data$qq_rank, qq_n + 1 - .data$qq_rank),
      p_upper = qbeta(0.975, .data$qq_rank, qq_n + 1 - .data$qq_rank)
    )

    output <- ggplot(
      qq_data,
      aes(
        x = .data$expected_p,
        y = .data$observed_p_plot,
        colour = .data$method_label,
        group = .data$method_label
      )
    ) +
      geom_ribbon(
        data = qq_band_data,
        mapping = aes(
          x = .data$expected_p,
          ymin = .data$p_lower,
          ymax = .data$p_upper
        ),
        inherit.aes = FALSE,
        fill = "grey75",
        colour = NA,
        alpha = 0.5
      ) +
      katlabutils::stat_qq_points(
        ymin = qq_panel_floor,
        size = 1.25,
        max_pts_to_plot = 1000
      ) +
      geom_abline(
        intercept = 0,
        slope = 1,
        colour = "black",
        linewidth = 0.55
      ) +
      scale_colour_manual(values = method_colors, name = NULL) +
      scale_x_continuous(
        trans = katlabutils::revlog_trans(10),
        breaks = qq_breaks,
        labels = qq_labels
      ) +
      scale_y_continuous(
        trans = katlabutils::revlog_trans(10),
        breaks = qq_breaks,
        labels = qq_labels
      ) +
      coord_cartesian(
        xlim = c(qq_display_floor, 1),
        ylim = c(qq_display_floor, 1),
        expand = FALSE
      ) +
      labs(
        x = if (show_x) "Expected null p-value" else NULL,
        y = "Observed p-value"
      ) +
      theme_spacrt(base_size = 11, remove_major_grid = TRUE) +
      theme(
        legend.position = "none",
        plot.margin = margin(3, 3, 3, 3)
      )
    if (!show_x) {
      output <- output + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    output
  }

  make_focused_typeI <- function(null_design_value, scenario_value) {
    plot_data <- focused_typeI |>
      filter(
        .data$null_design == null_design_value,
        .data$scenario_id == scenario_value
      ) |>
      mutate(
        method_label = factor(
          .data$method_label,
          levels = focused_method_levels
        )
      )

    ggplot(
      plot_data,
      aes(x = .data$method_label, y = .data$estimate, fill = .data$method_label)
    ) +
      geom_col(width = 0.62) +
      geom_errorbar(
        aes(ymin = .data$mc_lower, ymax = .data$mc_upper),
        width = 0.14,
        linewidth = 0.5
      ) +
      geom_hline(
        yintercept = focused_alpha,
        linetype = "dashed",
        colour = "black",
        linewidth = 0.5
      ) +
      geom_text(
        aes(y = .data$mc_upper, label = sprintf("%.3f", .data$estimate)),
        vjust = -0.45,
        size = 3
      ) +
      scale_fill_manual(values = method_colors, name = NULL) +
      scale_y_continuous(
        limits = c(0, typeI_upper),
        breaks = seq(0, typeI_upper, by = 0.05),
        labels = percent_label,
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        x = NULL,
        y = "Rejection rate"
      ) +
      theme_spacrt(base_size = 11) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(3, 3, 3, 3)
      )
  }

  make_focused_power <- function(
      null_design_value, scenario_value, power_type_value, show_x,
      show_legend = FALSE) {
    plot_data <- focused_power |>
      filter(
        .data$null_design == null_design_value,
        .data$scenario_id == scenario_value,
        .data$power_type == power_type_value
      ) |>
      mutate(
        method_label = factor(
          .data$method_label,
          levels = focused_method_levels
        )
      )

    output <- ggplot(
      plot_data,
      aes(
        x = .data$signal_strength,
        y = .data$estimate,
        colour = .data$method_label,
        group = .data$method_label
      )
    ) +
      geom_hline(
        yintercept = focused_alpha,
        linetype = "dashed",
        colour = "black",
        linewidth = 0.5
      ) +
      geom_errorbar(
        aes(ymin = .data$mc_lower, ymax = .data$mc_upper),
        width = 0,
        linewidth = 0.45
      ) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 2.1) +
      scale_colour_manual(values = method_colors, name = NULL) +
      scale_x_continuous(
        breaks = sort(unique(plot_data$signal_strength)),
        labels = signal_label
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        labels = percent_label,
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        x = if (show_x) {
          "Signal magnitude |gamma|\n(log-fold decrease)"
        } else {
          NULL
        },
        y = "Rejection rate"
      ) +
      guides(
        colour = guide_legend(nrow = 1, byrow = TRUE)
      ) +
      theme_spacrt(base_size = 11) +
      theme(
        legend.position = if (show_legend) "bottom" else "none",
        legend.key.width = grid::unit(0.14, "in"),
        legend.text = element_text(size = 10.5),
        legend.spacing.x = grid::unit(0.025, "in"),
        plot.margin = margin(3, 3, 3, 3)
      )
    if (!show_x) {
      output <- output + theme(
        axis.title.x = element_blank()
      )
    }
    output
  }

  make_focused_figure <- function(null_design_value) {
    qq_common <- make_focused_qq(
      null_design_value, "common_nb", show_x = FALSE
    ) +
      labs(tag = "a")
    typeI_common <- make_focused_typeI(
      null_design_value, "common_nb"
    ) + labs(tag = "b")
    power_legend_source <- make_focused_power(
      null_design_value, "common_nb", "raw_power",
      show_x = FALSE, show_legend = TRUE
    )
    raw_power_common <- power_legend_source +
      labs(tag = "c") +
      theme(legend.position = "none")

    qq_rare <- make_focused_qq(
      null_design_value, "rare_nb", show_x = FALSE
    ) +
      labs(tag = "d")
    typeI_rare <- make_focused_typeI(
      null_design_value, "rare_nb"
    ) + labs(tag = "e")
    raw_power_rare <- make_focused_power(
      null_design_value, "rare_nb", "raw_power", show_x = FALSE
    ) + labs(tag = "f")

    qq_severe <- make_focused_qq(
      null_design_value, "severe_variance", show_x = TRUE
    ) + labs(tag = "g")
    typeI_severe <- make_focused_typeI(
      null_design_value, "severe_variance"
    ) + labs(tag = "h")
    raw_power_severe <- make_focused_power(
      null_design_value, "severe_variance", "raw_power", show_x = TRUE
    ) + labs(tag = "i")

    shared_legend <- cowplot::get_plot_component(
      power_legend_source,
      "guide-box-bottom"
    )

    patchwork::wrap_plots(
      A = patchwork::plot_spacer(),
      B = focused_text_panel("Null p-value QQ"),
      C = focused_text_panel("Type-I error (alpha = 0.05)"),
      D = focused_text_panel("Rejection rate (p <= 0.05)"),
      E = focused_text_panel(
        unname(focused_scenario_labels["common_nb"]),
        font_size = 10.5,
        rotation = 90
      ),
      F = qq_common,
      G = typeI_common,
      H = raw_power_common,
      I = focused_text_panel(
        unname(focused_scenario_labels["rare_nb"]),
        font_size = 10.5,
        rotation = 90
      ),
      J = qq_rare,
      K = typeI_rare,
      L = raw_power_rare,
      M = focused_text_panel(
        unname(focused_scenario_labels["severe_variance"]),
        font_size = 10.5,
        rotation = 90
      ),
      N = qq_severe,
      O = typeI_severe,
      P = raw_power_severe,
      Q = patchwork::wrap_elements(full = shared_legend),
      design = "
        ABCD
        EFGH
        ####
        IJKL
        ####
        MNOP
        QQQQ
      ",
      widths = c(0.20, 1, 1, 1),
      heights = c(0.14, 1, 0.13, 1, 0.13, 1, 0.18)
    )
  }

  for (design_value in c("MI", "CI")) {
    focused_figure <- make_focused_figure(design_value)
    save_plot_pair(
      focused_figure,
      paste0(tolower(design_value), "-focused-3x3"),
      width = 10.5,
      height = 10.3,
      output_directory = manuscript_directory,
      save_png = FALSE
    )
  }
}

###############################################################################
# Failure-rate figure.
###############################################################################

failure_plot_data <- failure_summary |>
  select(
    all_of(scenario_variables), "method_label",
    starts_with("analysis_failure_"),
    starts_with("degenerate_"),
    starts_with("invalid_p_")
  ) |>
  select(
    -any_of(c(
      "analysis_failure_rate_lower", "analysis_failure_rate_upper",
      "degenerate_rate_lower", "degenerate_rate_upper",
      "invalid_p_rate_lower", "invalid_p_rate_upper"
    ))
  ) |>
  pivot_longer(
    cols = c(
      "analysis_failure_rate", "degenerate_rate", "invalid_p_rate"
    ),
    names_to = "failure_metric",
    values_to = "rate"
  ) |>
  mutate(
    lower = case_when(
      .data$failure_metric == "analysis_failure_rate" ~ .data$analysis_failure_lower,
      .data$failure_metric == "degenerate_rate" ~ .data$degenerate_lower,
      TRUE ~ .data$invalid_p_lower
    ),
    upper = case_when(
      .data$failure_metric == "analysis_failure_rate" ~ .data$analysis_failure_upper,
      .data$failure_metric == "degenerate_rate" ~ .data$degenerate_upper,
      TRUE ~ .data$invalid_p_upper
    ),
    failure_label = factor(
      recode(
        .data$failure_metric,
        analysis_failure_rate = "Analysis failure",
        degenerate_rate = "Degenerate score",
        invalid_p_rate = "Invalid/missing p-value"
      ),
      levels = c(
        "Analysis failure", "Degenerate score", "Invalid/missing p-value"
      )
    ),
    method_label = factor(
      .data$method_label,
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_grid_plot_labels()

if (nrow(failure_plot_data) > 0L) {
  failure_line_data <- failure_plot_data |>
    group_by(
      .data$failure_label, .data$null_design_label,
      .data$scenario_label, .data$method_label
    ) |>
    filter(n_distinct(.data$signal_strength) > 1L) |>
    ungroup()

  failure_plot <- ggplot(
    failure_plot_data,
    aes(
      x = .data$signal_strength,
      y = .data$rate,
      colour = .data$method_label,
      group = .data$method_label
    )
  ) +
    geom_errorbar(
      aes(ymin = .data$lower, ymax = .data$upper),
      width = 0, linewidth = 0.45
    ) +
    geom_line(
      data = failure_line_data,
      linewidth = 0.7
    ) +
    geom_point(size = 2.3) +
    facet_grid(
      rows = vars(.data$failure_label, .data$null_design_label),
      cols = vars(.data$scenario_label),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_colour_manual(values = method_colors, name = "Calibration") +
    scale_x_continuous(
      breaks = sort(unique(failure_plot_data$signal_strength))
    ) +
    scale_y_continuous(
      labels = percent_label,
      expand = expansion(mult = c(0.03, 0.15))
    ) +
    labs(
      x = "Signal magnitude |gamma| (log-fold decrease)",
      y = "Rate (95% MC interval)",
      title = "Failures and degenerate statistics"
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_spacrt(base_size = 10)

  save_plot_pair(
    failure_plot, "failure-rates",
    width = 13,
    height = 11
  )
}

###############################################################################
# Diagnostic-distribution figure (median and interquartile range).
###############################################################################

diagnostic_labels <- c(
  "statistic" = "Test statistic",
  "theta_hat" = "Fitted/working NB theta",
  "n_exposed" = "Number exposed",
  "n_exposed_nonzero" = "Exposed with nonzero count",
  "score_eff_n" = "Effective score sample size",
  "max_score_share" = "Largest squared-score share"
)

diagnostic_plot_data <- diagnostic_summary |>
  filter(is.finite(.data$median)) |>
  mutate(
    diagnostic_label = factor(
      unname(diagnostic_labels[.data$diagnostic]),
      levels = unname(diagnostic_labels[diagnostic_columns])
    ),
    method_label = factor(
      .data$method_label,
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_grid_plot_labels()

if (nrow(diagnostic_plot_data) > 0L) {
  diagnostic_line_data <- diagnostic_plot_data |>
    group_by(
      .data$diagnostic_label, .data$null_design_label,
      .data$scenario_label, .data$method_label
    ) |>
    filter(n_distinct(.data$signal_strength) > 1L) |>
    ungroup()

  diagnostic_plot <- ggplot(
    diagnostic_plot_data,
    aes(
      x = .data$signal_strength,
      y = .data$median,
      colour = .data$method_label,
      group = .data$method_label
    )
  ) +
    geom_errorbar(
      aes(ymin = .data$q25, ymax = .data$q75),
      width = 0, linewidth = 0.45
    ) +
    geom_line(
      data = diagnostic_line_data,
      linewidth = 0.7
    ) +
    geom_point(size = 2.2) +
    facet_grid(
      rows = vars(.data$diagnostic_label, .data$null_design_label),
      cols = vars(.data$scenario_label),
      scales = "free",
      space = "free_x"
    ) +
    scale_colour_manual(values = method_colors, name = "Calibration") +
    scale_x_continuous(
      breaks = sort(unique(diagnostic_plot_data$signal_strength))
    ) +
    labs(
      x = "Signal magnitude |gamma| (log-fold decrease)",
      y = "Median (interquartile range)",
      title = "Simulation diagnostics"
    ) +
    guides(
      colour = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_spacrt(base_size = 9, remove_major_grid = TRUE)

  save_plot_pair(
    diagnostic_plot, "diagnostics",
    width = 13,
    height = 18
  )
}

message(
  "Completed ", sim_name, " summaries in ",
  normalizePath(results_directory)
)
