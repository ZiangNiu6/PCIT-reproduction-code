#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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
regularization_output_directory <- Sys.getenv(
  "PCIT_REGULARIZATION_OUTPUT",
  unset = file.path(
    results_root, "simulation-results", "regularization-sensitivity"
  )
)
data_directory <- Sys.getenv(
  "PCIT_ROBUSTNESS_OUTPUT",
  unset = file.path(regularization_output_directory, "robustness-grid")
)
summary_file <- file.path(data_directory, "robustness-grid-summary.csv")
if (!file.exists(summary_file)) {
  stop("Missing robustness summary: ", summary_file)
}

source <- read.csv(summary_file, check.names = FALSE)
data <- source |>
  filter(.data$selection_rule == "Minimum CV") |>
  mutate(
    calibration = recode(
      .data$calibration,
      "Normal calibration" = "Normal",
      "Locked-tuning refit" = "Locked permutation"
    ),
    calibration = factor(
      .data$calibration,
      levels = c("Normal", "Locked permutation")
    ),
    learner = factor(.data$learner, levels = c("LRR", "KRR")),
    q_label = factor(
      paste0(
        "q = ", .data$q_penalized,
        "; hard residual rank = ", .data$hard_residual_rank
      ),
      levels = paste0(
        "q = ", sort(unique(.data$q_penalized)),
        "; hard residual rank = ",
        298L - sort(unique(.data$q_penalized))
      )
    ),
    label = sprintf("%.3f", .data$estimate)
  )

expected_rows <- 2L * 3L * 3L * 2L
if (nrow(data) != expected_rows ||
    any(!is.finite(data$estimate)) ||
    any(data$R != 1000L) ||
    any(data$B[data$calibration == "Locked permutation"] != 499L)) {
  stop("The robustness grid is incomplete or does not use R=1000, B=499.")
}

calibration_colors <- c(
  "Normal" = "#E69F00",
  "Locked permutation" = "#0072B2"
)
calibration_shapes <- c(
  "Normal" = 16,
  "Locked permutation" = 17
)

figure <- ggplot(
  data,
  aes(
    x = .data$signal_variance,
    y = .data$estimate,
    colour = .data$calibration,
    shape = .data$calibration,
    group = .data$calibration
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "red",
    linewidth = 0.65
  ) +
  geom_line(linewidth = 0.85, position = position_dodge(width = 0.35)) +
  geom_errorbar(
    aes(ymin = .data$lower_95, ymax = .data$upper_95),
    width = 0.28,
    linewidth = 0.65,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(
    size = 3.1,
    stroke = 0.8,
    position = position_dodge(width = 0.35)
  ) +
  geom_text(
    aes(
      label = .data$label,
      vjust = if_else(
        .data$calibration == "Normal",
        -1.15,
        1.55
      )
    ),
    size = 3.5,
    show.legend = FALSE,
    position = position_dodge(width = 0.35)
  ) +
  facet_grid(
    rows = vars(.data$learner),
    cols = vars(.data$q_label),
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = sort(unique(data$signal_variance)),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_colour_manual(values = calibration_colors, name = NULL) +
  scale_shape_manual(values = calibration_shapes, name = NULL) +
  labs(
    x = expression("Conditional-mean signal variance " * V),
    y = expression("Type-I error at " * alpha * " = 0.05")
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key.width = grid::unit(1.25, "cm"),
    plot.margin = margin(8, 10, 4, 6)
  )

pdf_file <- file.path(data_directory, "regularization-robustness-grid.pdf")
png_file <- file.path(data_directory, "regularization-robustness-grid.png")
ggsave(
  pdf_file,
  figure,
  width = 13.2,
  height = 7.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  png_file,
  figure,
  width = 13.2,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)
manuscript_file <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables",
  "regularization-robustness-grid.pdf"
)
dir.create(dirname(manuscript_file), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(pdf_file, manuscript_file, overwrite = TRUE)) {
  stop("Could not copy manuscript figure to ", manuscript_file)
}
message("Saved ", normalizePath(pdf_file))
message("Saved ", normalizePath(png_file))
message("Copied ", normalizePath(manuscript_file))
