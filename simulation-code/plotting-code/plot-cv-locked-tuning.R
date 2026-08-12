#!/usr/bin/env Rscript

# Reproduce the LRR and KRR rows of the deleted Figure 1 with the controlled
# leading-rank-16, V = 7.5 simulation:
#
#   LRR: mean LOOCV curve | CV minimum | fixed-lambda error path;
#   KRR: mean LOOCV curve | CV minimum | fixed-lambda error path.
#
# The blue permutation bars in the middle panels use locked tuning: select the
# X-side tuning parameter once on the observed data, then hold that choice fixed
# while refitting every raw-X permutation. In the right panels, each candidate
# lambda is fixed in advance and reapplied to every raw-X permutation. No full
# re-CV arm is used.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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
results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
data_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)

lrr_file <- file.path(data_directory, "lrr-calibration-summary.csv")
krr_file <- file.path(
  data_directory,
  "krr-selected-calibration-summary.csv"
)
curve_file <- file.path(data_directory, "loocv-curves.csv")
lrr_fixed_path_file <- file.path(
  data_directory,
  "lrr-fixed-path-summary.csv"
)
krr_fixed_path_file <- file.path(
  data_directory,
  "krr-fixed-path-summary.csv"
)
required_files <- c(
  lrr_file,
  krr_file,
  curve_file,
  lrr_fixed_path_file,
  krr_fixed_path_file
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing figure input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

lrr_source <- read.csv(lrr_file, check.names = FALSE)
krr_source <- read.csv(krr_file, check.names = FALSE)
curves <- read.csv(curve_file, check.names = FALSE)
lrr_fixed_path_source <- read.csv(
  lrr_fixed_path_file,
  check.names = FALSE
)
krr_fixed_path_source <- read.csv(
  krr_fixed_path_file,
  check.names = FALSE
)

normal_color <- "#E69F00"
permutation_color <- "#0072B2"
calibration_colors <- c(
  "Normal" = normal_color,
  "Permutation" = permutation_color
)
curve_design <- unique(curves[, c(
  "q_penalized",
  "signal_variance",
  "signal_pattern",
  "signal_rank"
)])
if (nrow(curve_design) != 1L ||
    curve_design$q_penalized != 178L ||
    abs(curve_design$signal_variance - 7.5) > 1e-12 ||
    curve_design$signal_pattern != "leading" ||
    curve_design$signal_rank != 16L) {
  stop(
    "The LOOCV curves are not from the selected q=178, rank-16, V=7.5 setup."
  )
}

normalize_lrr <- function(data) {
  data |>
    transmute(
      learner = .data$learner,
      cv_rule = .data$selection_rule,
      method = .data$calibration,
      rejections = .data$rejections,
      R = .data$R,
      B = .data$B,
      estimate = .data$estimate,
      lower_95 = .data$lower_95_exact,
      upper_95 = .data$upper_95_exact
    )
}

normalize_krr <- function(data) {
  data |>
    transmute(
      learner = .data$learner,
      cv_rule = .data$cv_rule,
      method = .data$method,
      rejections = .data$rejections,
      R = .data$R,
      B = if_else(
        .data$method == "Normal calibration",
        NA_integer_,
        as.integer(.data$B)
      ),
      estimate = .data$estimate,
      lower_95 = .data$lower_95,
      upper_95 = .data$upper_95
    )
}

cv_results <- bind_rows(
  normalize_lrr(lrr_source),
  normalize_krr(krr_source)
) |>
  mutate(
    cv_rule = factor(
      .data$cv_rule,
      levels = c("Minimum CV", "+1SE")
    )
  ) |>
  arrange(.data$learner, .data$cv_rule, .data$method)

expected_cv_cells <- expand.grid(
  learner = c("LRR", "KRR"),
  cv_rule = c("Minimum CV", "+1SE"),
  method = c(
    "Normal calibration",
    "Locked-tuning refit"
  ),
  stringsAsFactors = FALSE
)
observed_cv_cells <- cv_results |>
  transmute(
    learner = as.character(.data$learner),
    cv_rule = as.character(.data$cv_rule),
    method = as.character(.data$method)
  )
if (nrow(cv_results) != nrow(expected_cv_cells) ||
    nrow(anti_join(expected_cv_cells, observed_cv_cells,
      by = c("learner", "cv_rule", "method")
    )) > 0L ||
    any(cv_results$R != 1000L) ||
    any(cv_results$B[cv_results$method == "Locked-tuning refit"] != 499L) ||
    any(!is.finite(cv_results$estimate))) {
  stop("The paired LRR/KRR calibration summaries are incomplete.")
}

# Keep the two empirical calibrations in a common, plot-ready file.
summary_file <- file.path(data_directory, "calibration-summary.csv")
write.csv(cv_results, summary_file, row.names = FALSE)

fixed_path_method_levels <- c(
  "Normal approximation",
  "Refitted permutation"
)

same_numeric_path <- function(left, right, tolerance = 1e-8) {
  if (length(left) != length(right) ||
      any(is.infinite(left) != is.infinite(right))) {
    return(FALSE)
  }
  finite <- is.finite(left) & is.finite(right)
  all(abs(left[finite] - right[finite]) <= tolerance)
}

validate_fixed_path <- function(data, learner) {
  data <- data |>
    arrange(.data$candidate)
  curve_path <- curves |>
    filter(
      .data$learner == .env$learner,
      .data$sample_type == "Observed X"
    ) |>
    arrange(.data$effective_dimension)
  expected_dimensions <- c(seq.int(120L, 290L, by = 10L), 295L, 298L)

  if (nrow(data) != 20L ||
      !identical(as.integer(data$candidate), seq_len(20L)) ||
      !same_numeric_path(data$effective_dimension, expected_dimensions) ||
      nrow(curve_path) != 20L ||
      !same_numeric_path(
        data$effective_dimension,
        curve_path$effective_dimension
      ) ||
      !same_numeric_path(data$penalty_ratio, curve_path$penalty_ratio) ||
      any(data$learner != learner) ||
      any(data$tuning_regime != "fixed lambda") ||
      any(data$R != 1000L) ||
      any(data$B != 499L) ||
      any(data$q_penalized != 178L) ||
      any(data$signal_pattern != "leading") ||
      any(data$signal_rank != 16L) ||
      any(abs(data$signal_variance - 7.5) > 1e-12)) {
    stop(
      "The fixed-lambda ", learner,
      " results do not match the displayed LOOCV candidate path."
    )
  }
  invisible(TRUE)
}

validate_fixed_path(lrr_fixed_path_source, "LRR")
validate_fixed_path(krr_fixed_path_source, "KRR")

normalize_lrr_fixed_path <- function(data) {
  common <- data |>
    arrange(.data$candidate)
  bind_rows(
    common |>
      transmute(
        learner = .data$learner,
        candidate = .data$candidate,
        effective_dimension = .data$effective_dimension,
        penalty_ratio = .data$penalty_ratio,
        penalty_absolute = .data$penalty_absolute,
        method = "Normal approximation",
        rejections = .data$normal_rejections,
        R = .data$R,
        B = NA_integer_,
        estimate = .data$normal_empirical,
        lower_95 = .data$normal_lower_95_exact,
        upper_95 = .data$normal_upper_95_exact
      ),
    common |>
      transmute(
        learner = .data$learner,
        candidate = .data$candidate,
        effective_dimension = .data$effective_dimension,
        penalty_ratio = .data$penalty_ratio,
        penalty_absolute = .data$penalty_absolute,
        method = "Refitted permutation",
        rejections = .data$permutation_rejections,
        R = .data$R,
        B = .data$B,
        estimate = .data$permutation_empirical,
        lower_95 = .data$permutation_lower_95_exact,
        upper_95 = .data$permutation_upper_95_exact
      )
  )
}

normalize_krr_fixed_path <- function(data) {
  common <- data |>
    arrange(.data$candidate)
  bind_rows(
    common |>
      transmute(
        learner = .data$learner,
        candidate = .data$candidate,
        effective_dimension = .data$effective_dimension,
        penalty_ratio = .data$penalty_ratio,
        penalty_absolute = .data$penalty_absolute,
        method = "Normal approximation",
        rejections = .data$normal_rejections,
        R = .data$R,
        B = NA_integer_,
        estimate = .data$normal_empirical,
        lower_95 = .data$normal_lower_95,
        upper_95 = .data$normal_upper_95
      ),
    common |>
      transmute(
        learner = .data$learner,
        candidate = .data$candidate,
        effective_dimension = .data$effective_dimension,
        penalty_ratio = .data$penalty_ratio,
        penalty_absolute = .data$penalty_absolute,
        method = "Refitted permutation",
        rejections = .data$locked_permutation_rejections,
        R = .data$R,
        B = .data$B,
        estimate = .data$locked_permutation_empirical,
        lower_95 = .data$locked_permutation_lower_95,
        upper_95 = .data$locked_permutation_upper_95
      )
  )
}

fixed_path_metadata <- bind_rows(
  lrr_fixed_path_source,
  krr_fixed_path_source
) |>
  select(
    "learner",
    "candidate",
    "tuning_regime",
    "n",
    "q_penalized",
    "q_total",
    "hard_residual_rank",
    "signal_pattern",
    "signal_rank",
    "signal_variance"
  )

fixed_path_results <- bind_rows(
  normalize_lrr_fixed_path(lrr_fixed_path_source),
  normalize_krr_fixed_path(krr_fixed_path_source)
) |>
  left_join(
    fixed_path_metadata,
    by = c("learner", "candidate")
  ) |>
  mutate(
    method = factor(.data$method, levels = fixed_path_method_levels),
    calibration = factor(
      if_else(
        .data$method == "Normal approximation",
        "Normal",
        "Permutation"
      ),
      levels = c("Normal", "Permutation")
    )
  ) |>
  arrange(.data$learner, .data$method, .data$candidate) |>
  relocate(
    "learner",
    "candidate",
    "tuning_regime",
    "n",
    "q_penalized",
    "q_total",
    "hard_residual_rank",
    "signal_pattern",
    "signal_rank",
    "signal_variance"
  )

if (nrow(fixed_path_results) != 80L ||
    any(!is.finite(fixed_path_results$estimate))) {
  stop("The combined fixed-lambda calibration path is incomplete.")
}

fixed_path_summary_file <- file.path(
  data_directory,
  "fixed-tuning-calibration-path.csv"
)
write.csv(fixed_path_results, fixed_path_summary_file, row.names = FALSE)

regular_method_levels <- c(
  "Normal approximation",
  "Refitted permutation"
)

make_regular_bar_data <- function(learner, cv_rule) {
  empirical <- cv_results |>
    filter(
      .data$learner == .env$learner,
      as.character(.data$cv_rule) == .env$cv_rule,
      .data$method %in% c(
        "Normal calibration",
        "Locked-tuning refit"
      )
    )
  if (nrow(empirical) != 2L) {
    stop("Missing bar data for ", learner, " / ", cv_rule, ".")
  }

  normal <- empirical |>
    filter(.data$method == "Normal calibration") |>
    transmute(
      method = "Normal approximation",
      rejection_rate = .data$estimate,
      calibration = "Normal"
    )
  locked <- empirical |>
    filter(.data$method == "Locked-tuning refit") |>
    transmute(
      method = "Refitted permutation",
      rejection_rate = .data$estimate,
      calibration = "Permutation"
    )

  bind_rows(normal, locked) |>
    mutate(
      method = factor(.data$method, levels = regular_method_levels),
      calibration = factor(
        .data$calibration,
        levels = c("Normal", "Permutation")
      ),
      setting = if_else(cv_rule == "Minimum CV", "CV minimum", "+1 SE")
    )
}

make_regular_bar_plot <- function(learner, cv_rule) {
  data <- make_regular_bar_data(learner, cv_rule)

  ggplot(
    data,
    aes(
      x = .data$method,
      y = .data$rejection_rate,
      fill = .data$calibration
    )
  ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.6
    ) +
    geom_col(
      width = 0.70,
      colour = "grey15",
      linewidth = 0.65
    ) +
    geom_text(
      aes(label = sprintf("%.3f", .data$rejection_rate)),
      vjust = -0.45,
      size = 3.8,
      show.legend = FALSE
    ) +
    facet_wrap(~ setting, nrow = 1) +
    scale_fill_manual(values = calibration_colors, guide = "none") +
    scale_x_discrete(labels = c(
      "Normal approximation" = "Normal",
      "Refitted permutation" = "Perm."
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = NULL,
      x = NULL,
      y = expression("Rejection rate at " * alpha * " = 0.05")
    ) +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x = element_text(size = 10.5),
      axis.title.x = element_blank(),
      strip.text = element_text(face = "bold", size = 10.5),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(3, 8, 3, 5.5)
    )
}

curve_titles <- c(
  "LRR" = "Linear ridge",
  "KRR" = "Kernel ridge"
)
curve_labels <- c(
  "LRR" = "Cross-validation with LRR",
  "KRR" = "Cross-validation with KRR"
)
dimension_breaks <- c(120, 160, 190, 220, 250, 298)

make_loocv_plot <- function(learner, row_label) {
  data <- curves |>
    filter(
      .data$learner == .env$learner,
      .data$sample_type == "Observed X"
    ) |>
    arrange(.data$effective_dimension) |>
    mutate(strip_label = curve_labels[[.env$learner]])
  if (nrow(data) == 0L) {
    stop("No observed-X LOOCV curve for ", learner, ".")
  }

  minimum_row <- data[which.min(data$mean_loocv_mse), , drop = FALSE]

  ggplot(
    data,
    aes(x = .data$effective_dimension, y = .data$mean_loocv_mse)
  ) +
    geom_line(colour = permutation_color, linewidth = 0.9) +
    geom_point(colour = permutation_color, size = 2.4) +
    geom_vline(
      xintercept = minimum_row$effective_dimension,
      linetype = "dashed",
      colour = normal_color,
      linewidth = 0.7
    ) +
    geom_point(
      data = minimum_row,
      colour = normal_color,
      size = 3.6
    ) +
    facet_wrap(~ strip_label, nrow = 1) +
    scale_x_continuous(
      breaks = dimension_breaks,
      labels = scales::label_number(accuracy = 1)
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      title = paste0(row_label, " ", curve_titles[[learner]]),
      x = expression("Residual effective dimension " * d[eff]),
      y = "LOOCV MSE"
    ) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
      strip.text = element_text(face = "bold", size = 10.5),
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.margin = margin(3, 8, 3, 5.5)
    )
}

fixed_path_colors <- c(
  "Normal approximation" = normal_color,
  "Refitted permutation" = permutation_color
)
fixed_path_labels <- c(
  "Normal approximation" = "Normal",
  "Refitted permutation" = "Perm."
)

make_fixed_path_plot <- function(learner) {
  data <- fixed_path_results |>
    filter(.data$learner == .env$learner) |>
    mutate(strip_label = "Fixed lambda path")
  if (nrow(data) != 40L) {
    stop("Missing fixed-lambda path data for ", learner, ".")
  }

  ggplot(
    data,
    aes(
      x = .data$effective_dimension,
      y = .data$estimate,
      colour = .data$method,
      group = .data$method
    )
  ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.6
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(
      aes(shape = .data$method),
      size = 2.35,
      show.legend = FALSE
    ) +
    facet_wrap(~ strip_label, nrow = 1) +
    scale_colour_manual(
      values = fixed_path_colors,
      breaks = c("Normal approximation", "Refitted permutation"),
      labels = fixed_path_labels[
        c("Normal approximation", "Refitted permutation")
      ],
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Normal approximation" = 16,
        "Refitted permutation" = 17
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = dimension_breaks,
      labels = scales::label_number(accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.025, 0.08))
    ) +
    labs(
      title = NULL,
      x = expression("Residual effective dimension " * d[eff]),
      y = expression("Rejection rate at " * alpha * " = 0.05")
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
      strip.text = element_text(face = "bold", size = 10.5),
      legend.position = "bottom",
      legend.text = element_text(size = 12.5),
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.25, "cm"),
      plot.margin = margin(3, 8, 3, 5.5)
    )
}

lrr_plot <- make_loocv_plot("LRR", "(a)") +
  make_regular_bar_plot("LRR", "Minimum CV") +
  make_fixed_path_plot("LRR") +
  plot_layout(widths = c(1, 1, 1))
krr_plot <- make_loocv_plot("KRR", "(b)") +
  make_regular_bar_plot("KRR", "Minimum CV") +
  make_fixed_path_plot("KRR") +
  plot_layout(widths = c(1, 1, 1))

figure_body <- (
  lrr_plot / krr_plot +
    plot_layout(heights = c(1, 1))
)
figure <- wrap_plots(
  figure_body,
  guide_area(),
  ncol = 1,
  heights = c(1, 0.07),
  guides = "collect"
)

pdf_file <- file.path(data_directory, "cv-locked-tuning-simulation.pdf")
png_file <- file.path(data_directory, "cv-locked-tuning-simulation.png")
standalone_file <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables",
  "cv-locked-tuning-simulation.pdf"
)
dir.create(dirname(standalone_file), recursive = TRUE, showWarnings = FALSE)

ggsave(
  pdf_file,
  figure,
  width = 12,
  height = 6.6,
  useDingbats = FALSE
)
ggsave(png_file, figure, width = 12, height = 6.6, dpi = 220)
if (!file.copy(pdf_file, standalone_file, overwrite = TRUE)) {
  stop("Could not copy manuscript figure to ", standalone_file)
}

message("Saved: ", normalizePath(pdf_file))
message("Saved: ", normalizePath(png_file))
message("Saved: ", normalizePath(summary_file))
message("Saved: ", normalizePath(fixed_path_summary_file))
message("Copied: ", normalizePath(standalone_file, mustWork = FALSE))
