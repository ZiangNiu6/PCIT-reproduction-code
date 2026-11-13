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
source(file.path(analysis_dir, "..", "reproduction-helpers.R"))
project_root <- reproduction_root(script_dir)
results_root <- reproduction_results_root(project_root)
source(file.path(analysis_dir, "method_config.R"))
if (!isTRUE(capabilities("cairo"))) {
  stop("Cairo PDF support is required to embed publication figure fonts.")
}
results_dir <- Sys.getenv(
  "FRANGIEH_SCORE_RESULTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_RESULTS_SUBDIR)
)
plots_dir <- Sys.getenv(
  "FRANGIEH_SCORE_PLOTS_DIR",
  unset = file.path(results_root, "frangieh-2021", PERTURB_PLOTS_SUBDIR)
)
figure_override <- Sys.getenv("PCIT_FIGURE_DIR", unset = "")
manuscript_directory <- if (nzchar(figure_override)) {
  reproduction_figure_dir(project_root)
} else file.path(results_root, PERTURB_FIGURES_SUBDIR)

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
{
  source(file.path(analysis_dir, "workflow-io.R"))
  path <- file.path(results_dir, "postprocessing-contract.rds")
  if (!file.exists(path)) stop("GB2024 postprocessing-contract.rds is missing; rerun data_postprocessing.R.")
  manifest <- readRDS(path)
  table_files <- file.path(results_dir, c(
    "score_pair_long.rds", "score_fdr_summary.rds",
    "score_fdr_aggregated.rds", "score_diagnostics.rds"
  ))
  if (!is.list(manifest) || !identical(manifest$schema, "frangieh_postprocessed_bundle_v1") ||
      !identical(manifest$result_schema, PERTURB_SCHEMA_VERSION) ||
      !identical(manifest$cohort_policy, PERTURB_COHORT_POLICY) ||
      !identical(manifest$tables, frangieh_file_manifest(table_files))) {
    stop("GB2024 postprocessed results do not match their schema/policy/hash contract.")
  }
  frangieh_validate_preprocessing_contract(manifest$preprocessing_contract)
  replicate_files <- sort(list.files(
    results_dir, pattern = "^score_fdr_rep.*[.]rds$", full.names = TRUE
  ))
  expected_replicates <- sprintf("score_fdr_rep%03d.rds", 1:20)
  if (!identical(basename(replicate_files), expected_replicates) ||
      any(file.info(replicate_files)$isdir) ||
      !identical(manifest$replicate_files, frangieh_file_manifest(replicate_files))) {
    stop("GB2024 source replicates changed or are incomplete; rerun data_postprocessing.R before plotting.")
  }
}

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
         "data_postprocessing.R before plotting ",
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
    "mean_power", "se_power", "n_reps"),
  basename(aggregated_file)
)
validate_postprocessed(
  pair_long,
  c("schema_version", "rep", "pair_idx", "is_alt", "method", "q",
    "p_value"),
  basename(pair_file)
)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manuscript_directory, showWarnings = FALSE, recursive = TRUE)

HEADLINE_Q <- PERTURB_HEADLINE_Q
if (anyNA(aggregated$q) || anyNA(pair_long$q) ||
    !all(dplyr::near(aggregated$q, HEADLINE_Q)) ||
    !all(dplyr::near(pair_long$q, HEADLINE_Q))) {
  stop(sprintf("Publication inputs must contain only the configured q=%.2f results.", HEADLINE_Q))
}
if (nrow(aggregated) != length(PERTURB_METHOD_LEVELS)) {
  stop(sprintf("Expected one q=%.2f aggregate row for each configured method", HEADLINE_Q))
}
if (!all(aggregated$n_reps == 20L) || nrow(pair_long) != 100000L ||
    anyDuplicated(pair_long[c("rep", "method", "pair_idx")]) ||
    !setequal(unique(pair_long$rep), 1:20) ||
    anyNA(pair_long$p_value) || any(!is.finite(pair_long$p_value)) ||
    any(pair_long$p_value < 0 | pair_long$p_value > 1)) {
  stop("Expected 20 complete runs of 1,000 pairs for each of the five methods.")
}
block_counts <- pair_long |>
  group_by(rep, method) |>
  summarise(n_alt = sum(is_alt), n_null = sum(!is_alt), .groups = "drop")
if (nrow(block_counts) != 100L ||
    any(block_counts$n_alt != 50L | block_counts$n_null != 950L)) {
  stop("Each method/run must contain exactly 50 positive and 950 synthetic-null pairs.")
}

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
power_data <- aggregated |>
  transmute(
    method,
    metric = factor("Power", levels = "Power"),
    value = mean_power,
    se = se_power
  )

{
  for (data in list(fdr_data, power_data)) {
    if (any(!is.finite(data$value)) || any(data$value < 0 | data$value > 1) ||
        any(!is.finite(data$se)) || any(data$se < 0)) {
      stop("GB2024 plots require finite probability estimates and nonnegative standard errors.")
    }
  }
  use_log_fdr <- all(fdr_data$value > 0 & fdr_data$lower > 0)
  fdr_data <- fdr_data |>
    mutate(lower = pmax(0, lower), upper = pmin(1, upper))
  if (use_log_fdr) {
    log_range <- range(log10(c(fdr_data$lower, fdr_data$upper, HEADLINE_Q)))
    padding <- max(diff(log_range), 0.25) * 0.08
    fdr_limits <- 10^(log_range + c(-padding, padding))
    fdr_limits <- c(max(.Machine$double.xmin, fdr_limits[[1L]]),
                    min(1, fdr_limits[[2L]]))
    fdr_scale <- scale_y_log10(
      limits = fdr_limits, labels = scales::label_percent(),
      expand = expansion(mult = c(0.03, 0.06))
    )
    fdr_axis_label <- "False discovery rate (log scale)"
  } else {
    fdr_scale <- scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    )
    fdr_axis_label <- "False discovery rate (linear scale)"
  }
  power_limits <- c(0, 1)
  power_breaks <- seq(0, 1, by = 0.2)
  power_data <- power_data |>
    mutate(lower = pmax(0, value - se), upper = pmin(1, value + se))
}

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
    colour = "black"
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
  fdr_scale +
  labs(x = NULL, y = fdr_axis_label) +
  panel_theme +
  theme(legend.position = "none")

p_power <- ggplot(power_data, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.25
  ) +
  facet_wrap(~metric) +
  scale_fill_manual(
    name = NULL,
    values = display_colors,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_y_continuous(
    limits = power_limits,
    breaks = power_breaks,
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
  device = grDevices::cairo_pdf
)

manuscript_summary <- file.path(
  manuscript_directory,
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
  scale_x_continuous(labels = function(x) sub("^0[.]", ".", as.character(x))) +
  scale_fill_manual(values = display_colors, guide = "none") +
  labs(
    x = sprintf("p-value (adaptive methods calibrated at q = %.2f)", HEADLINE_Q),
    y = "count"
  ) +
  base_theme

pvalue_output <- file.path(plots_dir, "pvalue_histograms.pdf")
ggsave(
  pvalue_output,
  p_histograms,
  width = 13,
  height = 5.5,
  device = grDevices::cairo_pdf
)

manuscript_pvalues <- file.path(
  manuscript_directory,
  "frangieh-score5-pvalue-histograms.pdf"
)
if (!file.copy(pvalue_output, manuscript_pvalues, overwrite = TRUE)) {
  stop("Failed to synchronize the Frangieh p-value figure to the manuscript")
}

## A compact report of the selected analysis only ----------------------------
html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}
summary_rows <- vapply(seq_len(nrow(aggregated)), function(i) {
  row <- aggregated[i, ]
  sprintf(
    "<tr><th scope='row'>%s</th><td>%.2f%% &plusmn; %.2f%%</td><td>%.2f%% &plusmn; %.2f%%</td><td>%.2f</td><td>%d</td></tr>",
    html_escape(row$method), 100 * row$mean_fdr, 100 * row$se_fdr,
    100 * row$mean_power, 100 * row$se_power, row$mean_rejected, row$n_reps
  )
}, character(1L))
report_output <- file.path(plots_dir, "report.html")
writeLines(c(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  sprintf("<title>Frangieh five-method benchmark: q = %.2f</title>", HEADLINE_Q),
  "<style>body{font:17px/1.55 system-ui,sans-serif;max-width:1100px;margin:3rem auto;padding:0 1.2rem;color:#18212b}table{border-collapse:collapse;width:100%;font-size:15px}th,td{padding:.65rem .8rem;border-bottom:1px solid #d5dce3;text-align:right}th:first-child{text-align:left}thead{background:#edf2f6}a{color:#1261a0}</style>",
  "</head><body>",
  sprintf("<h1>Frangieh five-method benchmark at q = %.2f</h1>", HEADLINE_Q),
  "<p>Each of 20 benchmark runs contains 50 positive self-target pairs and 950 synthetic-null pairs. The same biological cells and positive pairs are reused across runs; these are not biological replicates.</p>",
  "<table><thead><tr><th>Method</th><th>Estimated FDR &plusmn; SE</th><th>Power &plusmn; SE</th><th>Mean discoveries</th><th>Runs</th></tr></thead><tbody>",
  summary_rows, "</tbody></table>",
  "<p>Estimated FDR is the mean of the run-level false discovery proportions, with zero assigned when a run has no discoveries. Power is the mean fraction of the 50 designated positive pairs recovered. Standard errors describe variability across these benchmark runs and do not quantify biological replication.</p>",
  sprintf("<p>All adaptive procedures were run at the selected target q = %.2f. The summaries use their saved decisions. The p-value display pools the q-specific values across runs; it does not re-threshold stopped adaptive p-values.</p>", HEADLINE_Q),
  "<p><a href='fdr_power_summary.pdf'>FDR and power figure</a> &middot; <a href='pvalue_histograms.pdf'>Positive and null p-value distributions</a></p>",
  "</body></html>"
), report_output)

cat("Wrote:\n")
cat("  ", file.path(plots_dir, "fdr_power_summary.pdf"), "\n", sep = "")
cat("  ", file.path(plots_dir, "pvalue_histograms.pdf"), "\n", sep = "")
cat("  ", manuscript_summary, "\n", sep = "")
cat("  ", manuscript_pvalues, "\n", sep = "")
cat("  ", report_output, "\n", sep = "")
