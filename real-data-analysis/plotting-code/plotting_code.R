#!/usr/bin/env Rscript
## Plot the five-method Frangieh score-calibration experiment.
##
## Reads corrected postprocessed RDS files and writes the corresponding plots.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
analysis_dir <- normalizePath(file.path(script_dir, ".."))
project_root <- normalizePath(file.path(analysis_dir, "..", ".."))
results_root <- Sys.getenv(
  "LOCAL_PDRT_DATA_DIR",
  unset = file.path(project_root, "results")
)
source(file.path(analysis_dir, "method_config.R"))
results_dir <- Sys.getenv(
  "FRANGIEH_SCORE_RESULTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_RESULTS_SUBDIR)
)
plots_dir <- Sys.getenv(
  "FRANGIEH_SCORE_PLOTS_DIR",
  unset = file.path(results_root, "frangieh-2021", "plots")
)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

aggregated_file <- file.path(results_dir, "score_fdr_aggregated.rds")
pair_file <- file.path(results_dir, "score_pair_long.rds")
if (!file.exists(aggregated_file) || !file.exists(pair_file)) {
  stop(
    "Missing five-method postprocessed results. Run data_postprocessing.R ",
    "after all score_fdr_repNNN.rds jobs finish."
  )
}

aggregated <- readRDS(aggregated_file)
pair_long <- readRDS(pair_file)

validate_postprocessed <- function(x, required_columns, label) {
  if (!is.data.frame(x) || !all(required_columns %in% names(x))) {
    stop(label, " is not a valid ", PERTURB_SCHEMA_VERSION,
         " data frame; missing: ",
         paste(setdiff(required_columns, names(x)), collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(x$schema_version) ||
      !all(x$schema_version == PERTURB_SCHEMA_VERSION)) {
    stop(label, " has an incompatible schema_version. Rerun ",
         "data_postprocessing.R; do not mix legacy combined_* outputs with ",
         PERTURB_SCHEMA_VERSION, ".", call. = FALSE)
  }
  found_methods <- unique(as.character(x$method))
  if (!setequal(found_methods, PERTURB_METHOD_LEVELS)) {
    stop(label, " does not contain exactly the five configured methods",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_postprocessed(
  aggregated,
  c("schema_version", "method", "q", "mean_fdr", "se_fdr",
    "mean_power", "se_power"),
  basename(aggregated_file)
)
validate_postprocessed(
  pair_long,
  c("schema_version", "rep", "pair_idx", "is_alt", "method", "q",
    "p_value"),
  basename(pair_file)
)

HEADLINE_Q <- 0.05
aggregated <- aggregated |>
  filter(dplyr::near(q, HEADLINE_Q))
pair_long <- pair_long |>
  filter(dplyr::near(q, HEADLINE_Q))
if (nrow(aggregated) != length(PERTURB_METHOD_LEVELS)) {
  stop("Expected one q=0.05 aggregate row for each configured method")
}
if (nrow(pair_long) == 0L) stop("No pair-level rows found at q=0.05")

aggregated$method <- factor(
  aggregated$method,
  levels = PERTURB_METHOD_LEVELS,
  labels = unname(PERTURB_METHOD_DISPLAY[PERTURB_METHOD_LEVELS])
)
pair_long$method <- factor(
  pair_long$method,
  levels = PERTURB_METHOD_LEVELS,
  labels = unname(PERTURB_METHOD_DISPLAY[PERTURB_METHOD_LEVELS])
)
display_colors <- setNames(
  unname(PERTURB_METHOD_COLORS[PERTURB_METHOD_LEVELS]),
  unname(PERTURB_METHOD_DISPLAY[PERTURB_METHOD_LEVELS])
)

cat(sprintf("Plotting q = %.2f with %d methods.\n",
            HEADLINE_Q, length(PERTURB_METHOD_LEVELS)))

base_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

## FDR and power panels ------------------------------------------------------
## A bar has a zero baseline, which is undefined under a strict logarithmic
## transformation.  Show FDR as point estimates with standard-error bars on a
## log10 axis, while retaining the original linear bar encoding for power.
fdr_data <- aggregated |>
  transmute(
    method,
    metric = factor("FDR", levels = "FDR"),
    value = mean_fdr,
    se = se_fdr,
    lower = mean_fdr - se_fdr,
    upper = mean_fdr + se_fdr
  )
if (any(!is.finite(fdr_data$value)) || any(fdr_data$value <= 0) ||
    any(!is.finite(fdr_data$lower)) || any(fdr_data$lower <= 0)) {
  stop(
    "The log-scaled FDR panel requires strictly positive finite aggregate ",
    "FDR estimates and lower standard-error bounds."
  )
}

power_data <- aggregated |>
  transmute(
    method,
    metric = factor("Power", levels = "Power"),
    value = mean_power,
    se = se_power
  )

panel_theme <- base_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key = element_rect(colour = NA, fill = NA)
  )

p_fdr <- ggplot(fdr_data, aes(x = method, y = value, fill = method)) +
  geom_hline(
    yintercept = HEADLINE_Q,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.16,
    linewidth = 0.75
  ) +
  geom_point(
    shape = 21,
    size = 5,
    stroke = 0.75,
    colour = "black"
  ) +
  facet_wrap(~metric) +
  scale_fill_manual(
    name = NULL,
    values = display_colors,
    guide = "none"
  ) +
  scale_y_log10(
    limits = c(0.04, 0.8),
    breaks = c(0.05, 0.10, 0.20, 0.50),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(x = NULL, y = "False discovery rate (log scale)") +
  panel_theme +
  theme(legend.position = "none")

p_power <- ggplot(power_data, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(0, value - se), ymax = value + se),
    width = 0.25
  ) +
  facet_wrap(~metric) +
  scale_fill_manual(
    name = NULL,
    values = display_colors,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_y_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  panel_theme

shared_legend <- cowplot::get_plot_component(
  p_power +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal"
    ),
  "guide-box-bottom"
)
p_summary <- patchwork::wrap_plots(
  A = p_fdr + theme(legend.position = "none"),
  B = p_power + theme(legend.position = "none"),
  C = patchwork::wrap_elements(full = shared_legend),
  design = "
    AB
    CC
  ",
  widths = c(1, 1),
  heights = c(1, 0.14)
)

summary_output <- file.path(plots_dir, "fdr_power_summary.pdf")
ggsave(
  summary_output,
  p_summary,
  width = 15,
  height = 6,
  device = grDevices::pdf,
  useDingbats = FALSE
)

manuscript_summary <- file.path(
  project_root,
  "project-progress", "combine", "figures-and-tables",
  "frangieh-score5-fdr-power.pdf"
)
if (!file.copy(summary_output, manuscript_summary, overwrite = TRUE)) {
  stop("Failed to synchronize the Frangieh FDR/power figure to the manuscript")
}

## Null and alternative p-value histograms ---------------------------------
first_rep <- min(pair_long$rep)
first_method <- levels(pair_long$method)[[1L]]
first_block <- pair_long |>
  filter(rep == first_rep, method == first_method)
n_alt_per_rep <- sum(first_block$is_alt)
n_null_per_rep <- sum(!first_block$is_alt)

histogram_data <- pair_long |>
  mutate(
    method_panel = factor(
      method,
      levels = levels(method),
      labels = c(
        "Parametric",
        "Information perm.",
        "Empirical perm.",
        "HGF flip",
        "Standardized flip"
      )
    ),
    group = ifelse(
      is_alt,
      sprintf("Alt (n=%d/rep)", n_alt_per_rep),
      sprintf("Null (n=%d/rep)", n_null_per_rep)
    )
  )

p_histograms <- ggplot(histogram_data, aes(x = p_value, fill = method)) +
  geom_histogram(bins = 30, boundary = 0, colour = NA) +
  facet_grid(group ~ method_panel, scales = "free_y") +
  scale_fill_manual(values = display_colors, guide = "none") +
  labs(x = "p-value", y = "count") +
  base_theme

pvalue_output <- file.path(plots_dir, "pvalue_histograms.pdf")
ggsave(
  pvalue_output,
  p_histograms,
  width = 13,
  height = 5.5,
  device = grDevices::pdf,
  useDingbats = FALSE
)

manuscript_pvalues <- file.path(
  project_root,
  "project-progress", "combine", "figures-and-tables",
  "frangieh-score5-pvalue-histograms.pdf"
)
if (!file.copy(pvalue_output, manuscript_pvalues, overwrite = TRUE)) {
  stop("Failed to synchronize the Frangieh p-value figure to the manuscript")
}

cat("Wrote:\n")
cat("  ", file.path(plots_dir, "fdr_power_summary.pdf"), "\n", sep = "")
cat("  ", file.path(plots_dir, "pvalue_histograms.pdf"), "\n", sep = "")
cat("  ", manuscript_summary, "\n", sep = "")
cat("  ", manuscript_pvalues, "\n", sep = "")
