#!/usr/bin/env Rscript

###############################################################################
# Plotting code for pcm-additive simulation
#
# Single plot: 3 rows x 3 columns
#   Row 1: X independent of Z (Additive Gaussian)
#   Row 2: Additive Gaussian (X dep Z)
#   Row 3: Heteroscedastic (X dep Z)
#   Column 1: Null Type-I error at the row's power-analysis sample size
#   Column 2: Rejection rate (raw), varying signal strength
#   Column 3: Empirical-cutoff rejection rate, varying signal strength
#
#   Resizing: for each (method, dgp), the empirical alpha-quantile of p-values
#   at gamma=0 is used as the rejection cutoff for all gamma values.
#
#   colour/fill = procedure (Permutation, Asymptotic)
#   line/border type = test (GCM solid, PCM dashed)
#
# The panel layout and typography follow manuscript Figure 3: a header row,
# three scenario rows, panel tags a--i, left-side row labels, and one shared
# bottom legend.
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

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
reproduction_root <- file.path(project_root, "reproduction-code")
source(file.path(reproduction_root, "plotting-style.R"))

cli_arguments <- commandArgs(trailingOnly = TRUE)
supported_arguments <- "--validate-only"
unknown_arguments <- setdiff(cli_arguments, supported_arguments)
if (length(unknown_arguments) > 0L) {
  stop(
    "Unknown argument(s): ", paste(unknown_arguments, collapse = ", "),
    ". Supported argument: --validate-only"
  )
}
validate_only <- "--validate-only" %in% cli_arguments

manuscript_directory <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables"
)

############################# Preprocess the data ##############################
sim_name <- "pcm-additive"
result_namespace <- "pcm-additive-matched"
result_schema <- "pcm_additive_matched"
combined_method <- "matched_gcm_pcm_rf"
display_methods <- c(
  "perm_gcm_rf", "asym_gcm_rf", "perm_pcm_rf", "asym_pcm_rf"
)
results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
results_dir <- file.path(
  results_root,
  "simulation-results",
  result_namespace,
  paste0(sim_name, "_results.rds")
)

if (!file.exists(results_dir)) {
  stop(
    "Matched pcm-additive result not found at ", results_dir, ". ",
    "Run `bash run_all_simulation.sh full`; the legacy pcm-additive result ",
    "is intentionally not accepted."
  )
}
results_raw <- readRDS(results_dir)

grid_file <- file.path(
  dirname(script_directory),
  "sim-spec-formal",
  "pcm-additive-grid.R"
)
grid_environment <- new.env(parent = globalenv())
sys.source(grid_file, envir = grid_environment)
current_parameter_grid <- as.data.frame(grid_environment$parameter_grid)
expected_dgps <- c("setting1_xindep", "setting1", "setting2")
current_grid_ok <-
  nrow(current_parameter_grid) == 27L &&
  !anyDuplicated(
    current_parameter_grid[c("n", "dgp", "gamma_signal")]
  ) &&
  setequal(current_parameter_grid$dgp, expected_dgps)
if (!current_grid_ok ||
    !all(c("grid_id", "B_perm", "alpha") %in% names(current_parameter_grid)) ||
    any(current_parameter_grid$B_perm != 499L) ||
    any(current_parameter_grid$alpha != 0.05)) {
  stop(
    "pcm-additive parameter grid is stale or incompatible; rerun ",
    "reproduction-code/simulation-code/sim-spec-formal/pcm-additive-grid.R"
  )
}

if (!is.list(results_raw) || !is.data.frame(results_raw$results) ||
    !all(c("grid_id", "run_id", "method", "output") %in%
         names(results_raw$results))) {
  stop("The saved matched pcm-additive result has an incompatible structure.")
}

results_table <- results_raw$results
if (!is.list(results_table$output) || any(vapply(
  results_table$output, is.null, logical(1)
))) {
  stop("Every matched pcm-additive row must contain a non-null output list.")
}
if (!is.numeric(results_table$grid_id) || !is.numeric(results_table$run_id) ||
    anyNA(results_table$grid_id) || anyNA(results_table$run_id) ||
    any(!is.finite(results_table$grid_id)) ||
    any(!is.finite(results_table$run_id)) ||
    any(results_table$grid_id != as.integer(results_table$grid_id)) ||
    any(results_table$run_id != as.integer(results_table$run_id))) {
  stop("grid_id and run_id must be finite, non-missing integers.")
}
results_table$grid_id <- as.integer(results_table$grid_id)
results_table$run_id <- as.integer(results_table$run_id)
results_table$method <- as.character(results_table$method)

if (!identical(unique(results_table$method), combined_method)) {
  stop(
    "The matched pcm-additive result must contain only method `",
    combined_method, "`."
  )
}

expected_grid_ids <- sort(as.integer(current_parameter_grid$grid_id))
result_grid_ids <- sort(unique(results_table$grid_id))
if (!identical(result_grid_ids, expected_grid_ids)) {
  stop(
    "The matched pcm-additive results do not match the current 27-cell grid."
  )
}

expected_replications <- 1000L
expected_rows <- length(expected_grid_ids) * expected_replications
if (nrow(results_table) != expected_rows) {
  stop(
    "The matched pcm-additive result must contain exactly ", expected_rows,
    " rows (27 grid cells x 1,000 replications); found ",
    nrow(results_table), "."
  )
}

result_keys <- paste(results_table$grid_id, results_table$run_id, sep = ":")
if (anyDuplicated(result_keys)) {
  duplicate_index <- which(duplicated(result_keys))[1L]
  stop(
    "Duplicate matched pcm-additive grid/run key: ",
    result_keys[[duplicate_index]], "."
  )
}
expected_run_ids <- seq_len(expected_replications)
run_ids_by_grid <- split(results_table$run_id, results_table$grid_id)
complete_run_ids <- vapply(
  run_ids_by_grid,
  function(run_ids) identical(sort(run_ids), expected_run_ids),
  logical(1)
)
if (length(complete_run_ids) != length(expected_grid_ids) ||
    !all(complete_run_ids)) {
  bad_grid <- names(complete_run_ids)[which(!complete_run_ids)[1L]]
  stop(
    "Every grid cell must contain run_id 1:1000 exactly once; first ",
    "incomplete grid_id: ", bad_grid, "."
  )
}

parameter_grid <- current_parameter_grid

stop_bad_output <- function(row_index, message) {
  stop(
    "Invalid matched output at grid_id=", results_table$grid_id[[row_index]],
    ", run_id=", results_table$run_id[[row_index]], ": ", message,
    call. = FALSE
  )
}

is_whole_number <- function(value) {
  is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value == as.integer(value)
}

p_value_matrix <- matrix(
  NA_real_,
  nrow = nrow(results_table),
  ncol = length(display_methods),
  dimnames = list(NULL, display_methods)
)
statistic_matrix <- matrix(
  NA_real_,
  nrow = nrow(results_table),
  ncol = 2L,
  dimnames = list(NULL, c("gcm", "pcm"))
)
reference_configuration <- NULL
n_by_grid <- stats::setNames(parameter_grid$n, parameter_grid$grid_id)

for (row_index in seq_len(nrow(results_table))) {
  output <- results_table$output[[row_index]]
  if (!is.list(output)) {
    stop_bad_output(row_index, "output is not a list")
  }
  if (!identical(output$schema_version, result_schema)) {
    stop_bad_output(
      row_index,
      paste0("schema_version must be `", result_schema, "`")
    )
  }
  if ("error" %in% names(output) &&
      !is.null(output$error) && length(output$error) > 0L) {
    stop_bad_output(row_index, "output contains an error")
  }

  p_values <- output$p_values
  if (!is.numeric(p_values) ||
      !identical(names(p_values), display_methods) ||
      anyNA(p_values) || any(!is.finite(p_values)) ||
      any(p_values < 0 | p_values > 1)) {
    stop_bad_output(
      row_index,
      paste0(
        "p_values must be finite probabilities named ",
        paste(display_methods, collapse = ", ")
      )
    )
  }

  statistics <- output$statistics
  if (!is.numeric(statistics) ||
      !identical(names(statistics), c("gcm", "pcm")) ||
      anyNA(statistics) || any(!is.finite(statistics))) {
    stop_bad_output(
      row_index,
      "statistics must contain finite named values gcm and pcm"
    )
  }

  alternatives <- output$alternatives
  expected_alternatives <- c(gcm = "two.sided", pcm = "greater")
  if (!is.character(alternatives) ||
      !identical(alternatives, expected_alternatives)) {
    stop_bad_output(
      row_index,
      "alternatives must be c(gcm='two.sided', pcm='greater')"
    )
  }

  matched <- output$matched
  if (!is.logical(matched) ||
      !identical(names(matched), c("gcm", "pcm", "all")) ||
      anyNA(matched) || !all(matched)) {
    stop_bad_output(row_index, "all matched-statistic flags must be TRUE")
  }
  flags <- output$flags
  expected_flags <- c(
    "matched_statistics", "explicit_pcm_split",
    "fixed_rf_randomness", "single_thread_rf"
  )
  if (!is.logical(flags) || !identical(names(flags), expected_flags) ||
      anyNA(flags) || !all(flags)) {
    stop_bad_output(row_index, "all implementation-integrity flags must be TRUE")
  }

  expected_asymptotic <- c(
    asym_gcm_rf = 2 * stats::pnorm(abs(statistics[["gcm"]]),
                                  lower.tail = FALSE),
    asym_pcm_rf = stats::pnorm(statistics[["pcm"]], lower.tail = FALSE)
  )
  observed_asymptotic <- p_values[names(expected_asymptotic)]
  if (!isTRUE(all.equal(
    unname(observed_asymptotic),
    unname(expected_asymptotic),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop_bad_output(
      row_index,
      "asymptotic p-values do not reconstruct from the stored statistics"
    )
  }

  configuration <- output$configuration
  configuration_fields <- c(
    "B_max", "h", "alpha", "rf_num_trees", "rf_num_threads", "n_cores"
  )
  if (!is.list(configuration) ||
      !all(configuration_fields %in% names(configuration))) {
    stop_bad_output(row_index, "configuration is missing required fields")
  }
  configuration_values <- vapply(
    configuration_fields,
    function(field) {
      value <- configuration[[field]]
      if (!is.numeric(value) || length(value) != 1L ||
          is.na(value) || !is.finite(value)) {
        stop_bad_output(
          row_index,
          paste0("configuration$", field, " must be a finite scalar")
        )
      }
      as.numeric(value)
    },
    numeric(1)
  )
  if (configuration_values[["B_max"]] != 499 ||
      configuration_values[["h"]] != 15 ||
      configuration_values[["alpha"]] != 0.05 ||
      any(configuration_values[c(
        "rf_num_trees", "rf_num_threads", "n_cores"
      )] <= 0)) {
    stop_bad_output(row_index, "configuration is not the canonical study")
  }
  if (is.null(reference_configuration)) {
    reference_configuration <- configuration_values
  } else if (!identical(configuration_values, reference_configuration)) {
    stop_bad_output(row_index, "configuration changes across replications")
  }

  seeds <- output$seeds
  if (!is.numeric(seeds) || length(seeds) == 0L || is.null(names(seeds)) ||
      any(!nzchar(names(seeds))) || anyDuplicated(names(seeds)) ||
      anyNA(seeds) || any(!is.finite(seeds)) ||
      any(seeds != as.integer(seeds))) {
    stop_bad_output(row_index, "seeds must be finite named integers")
  }

  pcm_split <- output$pcm_split
  if (!is.list(pcm_split) ||
      !all(c("D1_n", "D2_n", "D2_checksum") %in% names(pcm_split)) ||
      !is_whole_number(pcm_split$D1_n) ||
      !is_whole_number(pcm_split$D2_n) ||
      pcm_split$D1_n <= 0 || pcm_split$D2_n <= 0 ||
      pcm_split$D1_n + pcm_split$D2_n !=
        n_by_grid[[as.character(results_table$grid_id[[row_index]])]] ||
      length(pcm_split$D2_checksum) != 1L ||
      is.na(pcm_split$D2_checksum) ||
      !nzchar(as.character(pcm_split$D2_checksum))) {
    stop_bad_output(row_index, "pcm_split is missing or incompatible with n")
  }

  adaptive <- output$adaptive
  if (!is.list(adaptive) ||
      !identical(names(adaptive), c("gcm", "pcm"))) {
    stop_bad_output(row_index, "adaptive must contain gcm and pcm diagnostics")
  }
  for (test_name in c("gcm", "pcm")) {
    diagnostic <- adaptive[[test_name]]
    required_diagnostics <- c(
      "actual_B", "K", "t_star", "stopped_early", "stop_reason"
    )
    if (!is.list(diagnostic) ||
        !all(required_diagnostics %in% names(diagnostic)) ||
        !is_whole_number(diagnostic$actual_B) ||
        !is_whole_number(diagnostic$K) ||
        diagnostic$actual_B < 1L ||
        diagnostic$actual_B > configuration_values[["B_max"]] ||
        diagnostic$K < 0L || diagnostic$K > diagnostic$actual_B ||
        !is.numeric(diagnostic$t_star) ||
        length(diagnostic$t_star) != 1L || is.na(diagnostic$t_star) ||
        !is.logical(diagnostic$stopped_early) ||
        length(diagnostic$stopped_early) != 1L ||
        is.na(diagnostic$stopped_early) ||
        !identical(
          diagnostic$stopped_early,
          diagnostic$actual_B < configuration_values[["B_max"]]
        ) ||
        !is.character(diagnostic$stop_reason) ||
        length(diagnostic$stop_reason) != 1L ||
        is.na(diagnostic$stop_reason) || !nzchar(diagnostic$stop_reason)) {
      stop_bad_output(
        row_index,
        paste0("invalid adaptive diagnostics for ", test_name)
      )
    }

    if (diagnostic$K < configuration_values[["h"]]) {
      if (!is.infinite(diagnostic$t_star) || diagnostic$t_star < 0) {
        stop_bad_output(
          row_index,
          paste0(test_name, " t_star must be Inf before the h-th exceedance")
        )
      }
      expected_permutation <- configuration_values[["h"]] /
        (diagnostic$actual_B + configuration_values[["h"]] - diagnostic$K)
    } else {
      if (!is.finite(diagnostic$t_star) || diagnostic$t_star < 1 ||
          diagnostic$t_star > diagnostic$actual_B) {
        stop_bad_output(
          row_index,
          paste0(test_name, " t_star is incompatible with adaptive stopping")
        )
      }
      expected_permutation <- configuration_values[["h"]] /
        diagnostic$t_star
    }
    permutation_name <- paste0("perm_", test_name, "_rf")
    if (!isTRUE(all.equal(
      unname(p_values[[permutation_name]]),
      unname(expected_permutation),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      stop_bad_output(
        row_index,
        paste0(
          permutation_name,
          " does not reconstruct from the adaptive diagnostics"
        )
      )
    }
  }

  p_value_matrix[row_index, ] <- p_values[display_methods]
  statistic_matrix[row_index, ] <- statistics[c("gcm", "pcm")]
}

########################## Expand matched p-values ############################
pvals_df <- dplyr::bind_cols(
  results_table[c("grid_id", "run_id")],
  as.data.frame(p_value_matrix, check.names = FALSE)
) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(display_methods),
    names_to = "method",
    values_to = "p.value"
  ) |>
  dplyr::left_join(parameter_grid, by = "grid_id")

if (nrow(pvals_df) != expected_rows * length(display_methods) ||
    anyNA(pvals_df$p.value) || any(!is.finite(pvals_df$p.value))) {
  stop("Expansion of the matched p-values failed its integrity check.")
}

if (validate_only) {
  cat(
    "Validated matched pcm-additive result: ", expected_rows,
    " complete grid/run outputs and ", nrow(pvals_df),
    " finite displayed p-values.\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
}

alpha <- 0.05

add_wilson_bounds <- function(data, confidence = 0.95) {
  z_value <- stats::qnorm(1 - (1 - confidence) / 2)
  estimate <- data$n_reject / data$n_valid
  denominator <- 1 + z_value^2 / data$n_valid
  center <- (estimate + z_value^2 / (2 * data$n_valid)) / denominator
  half_width <- z_value * sqrt(
    estimate * (1 - estimate) / data$n_valid +
      z_value^2 / (4 * data$n_valid^2)
  ) / denominator

  data$mean_reject <- estimate
  data$reject_lower <- pmax(0, center - half_width)
  data$reject_upper <- pmin(1, center + half_width)
  data
}

########################## Compute standard rejection rates ####################
summary_df <- pvals_df |>
  dplyr::mutate(reject = as.integer(p.value <= alpha)) |>
  dplyr::group_by(method, dgp, n, gamma_signal) |>
  dplyr::summarise(
    n_reject = sum(reject, na.rm = TRUE),
    n_valid = sum(!is.na(reject)),
    .groups = "drop"
  ) |>
  add_wilson_bounds()

########################## Compute resized rejection rates #####################
# For each (method, dgp), find empirical alpha-quantile from gamma=0 p-values
# at the fixed n used in power analysis (n=600 for dep-Z, n=100 for xindep)

# Step 1: get null p-values for each (method, dgp)
null_pvals <- pvals_df |>
  dplyr::filter(gamma_signal == 0) |>
  dplyr::select(method, dgp, n, p.value)

# For power DGPs, we want the cutoff at the fixed n:
#   setting1/2 at n=600, setting1_xindep at n=100.  Each DGP also has its
#   own gamma=0 sample-size grid for the Type-I column.
power_dgps <- c("setting1_xindep", "setting1", "setting2")
power_n <- c(setting1_xindep = 100, setting1 = 600, setting2 = 600)

empirical_cutoffs <- null_pvals |>
  dplyr::filter(dgp %in% power_dgps) |>
  dplyr::filter(n == power_n[dgp]) |>
  dplyr::group_by(method, dgp) |>
  dplyr::summarise(
    cutoff = quantile(p.value, probs = alpha, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: join cutoffs back to all p-values and compute resized rejection
resized_df <- pvals_df |>
  dplyr::filter(dgp %in% power_dgps) |>
  dplyr::left_join(empirical_cutoffs, by = c("method", "dgp")) |>
  dplyr::mutate(reject_resized = as.integer(p.value <= cutoff)) |>
  dplyr::group_by(method, dgp, n, gamma_signal) |>
  dplyr::summarise(
    n_reject = sum(reject_resized, na.rm = TRUE),
    n_valid = sum(!is.na(reject_resized)),
    .groups = "drop"
  ) |>
  add_wilson_bounds()

########################## Add calibration / test columns ######################
# Keep the two calibrations adjacent within each test statistic so the
# permutation-versus-asymptotic change is read as a direct pairwise comparison.
method_order <- c(
  "perm_gcm_rf", "asym_gcm_rf", "perm_pcm_rf", "asym_pcm_rf"
)
method_display <- c(
  "perm_gcm_rf" = "Permutation GCM",
  "asym_gcm_rf" = "Asymptotic GCM",
  "perm_pcm_rf" = "Permutation PCM",
  "asym_pcm_rf" = "Asymptotic PCM"
)

add_method_cols <- function(data) {
  data |>
    add_calibration_test_aesthetics() |>
    dplyr::mutate(
      method = factor(.data$method, levels = method_order),
      method_label = factor(
        unname(method_display[as.character(.data$method)]),
        levels = unname(method_display[method_order])
      )
    )
}

summary_df <- add_method_cols(summary_df)
resized_df <- add_method_cols(resized_df)

###############################################################################
# Row facet labels (DGP settings)
###############################################################################
setting_levels <- c("X independent of Z",
                    "Additive Gaussian",
                    "Heteroscedastic")

###############################################################################
# Assemble data for each column, with a shared `setting` facet variable
###############################################################################

# --- Column 1: Type-I error at the row's power-analysis sample size ---
df_col1 <- bind_rows(
  summary_df |>
    dplyr::filter(
      dgp == "setting1_xindep", gamma_signal == 0,
      n == unname(power_n["setting1_xindep"])
    ) |>
    dplyr::mutate(setting = "X independent of Z"),
  summary_df |>
    dplyr::filter(
      dgp == "setting1", gamma_signal == 0,
      n == unname(power_n["setting1"])
    ) |>
    dplyr::mutate(setting = "Additive Gaussian"),
  summary_df |>
    dplyr::filter(
      dgp == "setting2", gamma_signal == 0,
      n == unname(power_n["setting2"])
    ) |>
    dplyr::mutate(setting = "Heteroscedastic")
) |>
  dplyr::mutate(setting = factor(setting, levels = setting_levels))

# --- Column 2: Raw rejection rate, x = signal strength ---
df_col2 <- bind_rows(
  summary_df |>
    dplyr::filter(dgp == "setting1_xindep", n == 100) |>
    dplyr::mutate(setting = "X independent of Z"),
  summary_df |>
    dplyr::filter(dgp == "setting1", n == 600) |>
    dplyr::mutate(setting = "Additive Gaussian"),
  summary_df |>
    dplyr::filter(dgp == "setting2", n == 600) |>
    dplyr::mutate(setting = "Heteroscedastic")
) |>
  dplyr::mutate(setting = factor(setting, levels = setting_levels))

# --- Column 3: Empirical-cutoff rejection rate, x = signal strength ---
df_col3 <- bind_rows(
  resized_df |>
    dplyr::filter(dgp == "setting1_xindep", n == 100) |>
    dplyr::mutate(setting = "X independent of Z"),
  resized_df |>
    dplyr::filter(dgp == "setting1", n == 600) |>
    dplyr::mutate(setting = "Additive Gaussian"),
  resized_df |>
    dplyr::filter(dgp == "setting2", n == 600) |>
    dplyr::mutate(setting = "Heteroscedastic")
) |>
  dplyr::mutate(setting = factor(setting, levels = setting_levels))

###############################################################################
# Build the Figure-3-style 3 x 3 layout
###############################################################################

if (nrow(df_col1) != length(setting_levels) * length(method_order)) {
  stop("Expected one Type-I bar per method and scenario at the anchor n")
}

typeI_upper <- max(
  0.15,
  ceiling(max(df_col1$reject_upper, na.rm = TRUE) * 1.15 / 0.05) * 0.05
)
percent_label <- scales::label_percent(accuracy = 1)
signal_label <- scales::label_number(accuracy = 0.001, trim = TRUE)

text_panel <- function(
    label, font_size = 11.5, font_face = "bold", rotation = 0) {
  patchwork::wrap_elements(full = grid::textGrob(
    label,
    rot = rotation,
    gp = grid::gpar(fontsize = font_size, fontface = font_face)
  ))
}

make_typeI_panel <- function(setting_value) {
  plot_data <- df_col1 |>
    dplyr::filter(.data$setting == setting_value)

  ggplot(
    plot_data,
    aes(
      x = .data$method_label,
      y = .data$mean_reject,
      fill = .data$calibration,
      linetype = .data$test
    )
  ) +
    geom_col(
      width = 0.62,
      colour = "black",
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = .data$reject_lower, ymax = .data$reject_upper),
      width = 0.14,
      colour = "black",
      linetype = "solid",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = alpha,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.5
    ) +
    geom_text(
      aes(y = .data$reject_upper, label = sprintf("%.3f", .data$mean_reject)),
      vjust = -0.45,
      size = 3
    ) +
    scale_calibration_fill(guide = "none") +
    scale_test_linetype(guide = "none") +
    scale_y_continuous(
      limits = c(0, typeI_upper),
      breaks = seq(0, typeI_upper, by = 0.05),
      labels = percent_label,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = "Rejection rate") +
    theme_figure3() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
}

make_rejection_panel <- function(
    data, setting_value, show_x = FALSE, show_legend = FALSE) {
  plot_data <- data |>
    dplyr::filter(.data$setting == setting_value)

  output <- ggplot(
    plot_data,
    aes(
      x = .data$gamma_signal,
      y = .data$mean_reject,
      colour = .data$calibration,
      linetype = .data$test,
      group = interaction(.data$calibration, .data$test)
    )
  ) +
    geom_hline(
      yintercept = alpha,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.5
    ) +
    geom_errorbar(
      aes(ymin = .data$reject_lower, ymax = .data$reject_upper),
      width = 0,
      linetype = "solid",
      linewidth = 0.45
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.1) +
    scale_calibration_colour(
      name = "Procedure",
      guide = guide_legend(
        nrow = 1,
        order = 1,
        override.aes = list(linetype = "solid")
      )
    ) +
    scale_test_linetype(
      guide = guide_legend(
        nrow = 1,
        order = 2,
        override.aes = list(colour = "black")
      )
    ) +
    scale_x_continuous(
      breaks = sort(unique(plot_data$gamma_signal)),
      labels = signal_label
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = percent_label,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = if (show_x) "Signal strength gamma" else NULL,
      y = "Rejection rate"
    ) +
    theme_figure3() +
    theme(
      legend.position = if (show_legend) "bottom" else "none"
    )

  if (!show_x) {
    output <- output + theme(axis.title.x = element_blank())
  }
  output
}

typeI_xindep <- make_typeI_panel("X independent of Z") + labs(tag = "a")
raw_xindep_legend <- make_rejection_panel(
  df_col2, "X independent of Z", show_legend = TRUE
)
raw_xindep <- raw_xindep_legend +
  labs(tag = "b") +
  theme(legend.position = "none")
resized_xindep <- make_rejection_panel(
  df_col3, "X independent of Z"
) + labs(tag = "c")

typeI_additive <- make_typeI_panel("Additive Gaussian") + labs(tag = "d")
raw_additive <- make_rejection_panel(
  df_col2, "Additive Gaussian"
) + labs(tag = "e")
resized_additive <- make_rejection_panel(
  df_col3, "Additive Gaussian"
) + labs(tag = "f")

typeI_heteroscedastic <- make_typeI_panel("Heteroscedastic") + labs(tag = "g")
raw_heteroscedastic <- make_rejection_panel(
  df_col2, "Heteroscedastic", show_x = TRUE
) + labs(tag = "h")
resized_heteroscedastic <- make_rejection_panel(
  df_col3, "Heteroscedastic", show_x = TRUE
) + labs(tag = "i")

shared_legend <- cowplot::get_plot_component(
  raw_xindep_legend +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    ),
  "guide-box-bottom"
)

joint_plot <- patchwork::wrap_plots(
  A = patchwork::plot_spacer(),
  B = text_panel("Type-I error (alpha = 0.05)"),
  C = text_panel("Rejection rate (p <= 0.05)"),
  D = text_panel("Empirical-cutoff rejection rate"),
  E = text_panel("X independent of Z (n = 100)", 10.5, rotation = 90),
  F = typeI_xindep,
  G = raw_xindep,
  H = resized_xindep,
  I = text_panel("Additive Gaussian (n = 600)", 10.5, rotation = 90),
  J = typeI_additive,
  K = raw_additive,
  L = resized_additive,
  M = text_panel("Heteroscedastic (n = 600)", 10.5, rotation = 90),
  N = typeI_heteroscedastic,
  O = raw_heteroscedastic,
  P = resized_heteroscedastic,
  Q = patchwork::wrap_elements(full = shared_legend),
  design = "
    ABCD
    EFGH
    IJKL
    MNOP
    QQQQ
  ",
  widths = c(0.20, 1, 1, 1),
  heights = c(0.14, 1, 1, 1, 0.18)
)

########################## Save and synchronize ###############################
dir.create(manuscript_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(manuscript_directory, "pcm_additive.pdf")
ggsave(
  output_path,
  joint_plot,
  height = 9.5,
  width = 10.5,
  device = cairo_pdf
)
cat("Plot saved to ", output_path, "\n", sep = "")
