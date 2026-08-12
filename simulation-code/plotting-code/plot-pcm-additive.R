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

manuscript_directory <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables"
)
dir.create(manuscript_directory, recursive = TRUE, showWarnings = FALSE)

############################# Preprocess the data ##############################
sim_name <- "pcm-additive"
results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
results_dir <- file.path(
  results_root,
  "simulation-results",
  sim_name,
  paste0(sim_name, "_results.rds")
)

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
if (!current_grid_ok) {
  stop(
    "pcm-additive parameter grid is stale or incompatible; rerun ",
    "reproduction-code/simulation-code/sim-spec-formal/pcm-additive-grid.R"
  )
}
expected_methods <- c(
  "perm_gcm_rf", "perm_pcm_rf", "asym_gcm_rf", "asym_pcm_rf"
)
if (!is.list(results_raw) || !is.data.frame(results_raw$results) ||
    !all(c("grid_id", "method", "output") %in% names(results_raw$results)) ||
    !setequal(as.character(unique(results_raw$results$method)), expected_methods)) {
  stop("The saved pcm-additive result has an incompatible structure.")
}

result_grid_ids <- sort(unique(results_raw$results$grid_id))
result_counts <- results_raw$results |>
  dplyr::count(.data$grid_id, .data$method, name = "n_replications")
expected_pair_count <- length(result_grid_ids) * length(expected_methods)
if (nrow(result_counts) != expected_pair_count ||
    any(result_counts$n_replications != 1000L)) {
  stop(
    "The pcm-additive result must contain exactly 1,000 replications for ",
    "every grid-by-method combination."
  )
}
if (setequal(result_grid_ids, current_parameter_grid$grid_id)) {
  parameter_grid <- current_parameter_grid
} else {
  stop(
    "The saved pcm-additive results do not match the redesigned 27-cell ",
    "three-row simulation. Archive the legacy result file and rerun the ",
    "simulation before plotting."
  )
}

if (!setequal(unique(results_raw$results$grid_id), parameter_grid$grid_id)) {
  stop("The selected pcm-additive result cells do not match the plotting grid.")
}

# join results with parameter grid
full_results <- results_raw$results |>
  dplyr::left_join(parameter_grid, by = "grid_id")

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

########################## Extract raw p-values ################################
pvals_df <- full_results |>
  dplyr::mutate(method = as.character(method)) |>
  tidyr::unnest_wider(output)

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
output_path <- file.path(manuscript_directory, "pcm_additive.pdf")
ggsave(
  output_path,
  joint_plot,
  height = 9.5,
  width = 10.5,
  device = cairo_pdf
)
cat("Plot saved to ", output_path, "\n", sep = "")
