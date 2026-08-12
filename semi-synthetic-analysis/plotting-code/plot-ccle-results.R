#!/usr/bin/env Rscript
###############################################################################
# Plot CCLE Synthetic-Y 30-targets x s in {5,10} sweep (60 setups).
#
# Reads the combined result under $LOCAL_PDRT_DATA_DIR/ccle and writes the QQ,
# standalone CHIC2 calibration preview, and combined Type-I figure directly to
# project-progress/combine/figures-and-tables. Test type is encoded by colour;
# calibration is encoded by line or border type. Pass --intro-only to
# regenerate only the standalone CHIC2 calibration preview used as Figure 1.
###############################################################################

script_args <- commandArgs(trailingOnly = TRUE)
unknown_args <- setdiff(script_args, "--intro-only")
if (length(unknown_args) > 0L) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "),
       call. = FALSE)
}
intro_only <- "--intro-only" %in% script_args

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
})

script_file_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_file_arg) != 1L) {
  stop("Could not determine the plotting script path.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_file_arg[[1L]]),
  mustWork = TRUE
)
script_dir <- dirname(script_path)
project_root <- normalizePath(
  file.path(script_dir, "..", "..", ".."),
  mustWork = TRUE
)

data_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(data_root)) {
  stop(
    "LOCAL_PDRT_DATA_DIR must point to the external data/results directory.",
    call. = FALSE
  )
}
data_root <- path.expand(data_root)
if (!dir.exists(data_root)) {
  stop("LOCAL_PDRT_DATA_DIR does not exist: ", data_root, call. = FALSE)
}
data_root <- normalizePath(data_root, mustWork = TRUE)

style_path <- file.path(project_root, "reproduction-code", "plotting-style.R")
if (!file.exists(style_path)) {
  stop("Shared plotting style not found: ", style_path, call. = FALSE)
}
source(style_path)

results_path <- file.path(
  data_root,
  "ccle",
  "ccle-typeI-synthY-30targets-results.rds"
)
if (!file.exists(results_path)) {
  stop(
    "Combined CCLE results not found: ", results_path,
    "\nRun combine-ccle-results.R after all 60 setup jobs finish.",
    call. = FALSE
  )
}
results <- readRDS(results_path)

storage_methods <- c("asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive")
required_record_fields <- c(
  "setup_idx", "target", "sparsity_tag", "freq", "s", "R", "pvals",
  "rejection_rate", "na_rate"
)
validate_record <- function(record, key) {
  if (!is.list(record) ||
      !all(required_record_fields %in% names(record))) {
    stop("Malformed CCLE record: ", key, call. = FALSE)
  }
  if (!identical(as.integer(record$R), 500L) ||
      !is.matrix(record$pvals) || nrow(record$pvals) != 500L ||
      !identical(colnames(record$pvals), storage_methods) ||
      !identical(names(record$rejection_rate), storage_methods) ||
      !identical(names(record$na_rate), storage_methods)) {
    stop("Noncanonical CCLE method/replication schema in: ", key,
         call. = FALSE)
  }
  if (length(record$setup_idx) != 1L || is.na(record$setup_idx) ||
      length(record$target) != 1L || is.na(record$target) ||
      length(record$sparsity_tag) != 1L || is.na(record$sparsity_tag) ||
      length(record$freq) != 1L || !is.finite(record$freq) ||
      length(record$s) != 1L || is.na(record$s) ||
      !identical(key, sprintf("%s_s%d", record$target, record$s))) {
    stop("Invalid CCLE setup metadata in: ", key, call. = FALSE)
  }
  if (any(!is.finite(record$pvals)) ||
      any(record$pvals < 0 | record$pvals > 1) ||
      any(!is.finite(record$rejection_rate)) ||
      any(!is.finite(record$na_rate)) || any(record$na_rate != 0)) {
    stop("Missing or invalid CCLE p-values in: ", key, call. = FALSE)
  }
  expected_rejection <- colMeans(record$pvals <= 0.05)
  if (!isTRUE(all.equal(
    unname(record$rejection_rate), unname(expected_rejection),
    tolerance = 1e-12
  ))) {
    stop("Stored CCLE rejection rates disagree with p-values in: ", key,
         call. = FALSE)
  }
  invisible(TRUE)
}

if (!is.list(results) || length(results) != 60L || is.null(names(results)) ||
    anyNA(names(results)) || any(names(results) == "") ||
    anyDuplicated(names(results))) {
  stop("The canonical CCLE result must contain 60 uniquely named setups.",
       call. = FALSE)
}
invisible(Map(validate_record, results, names(results)))
setup_ids <- vapply(results, function(x) as.integer(x$setup_idx), integer(1L))
targets <- vapply(results, function(x) as.character(x$target), character(1L))
s_values <- vapply(results, function(x) as.integer(x$s), integer(1L))
sparsity_tags <- vapply(
  results, function(x) as.character(x$sparsity_tag), character(1L)
)
if (!identical(unname(sort(setup_ids)), seq_len(60L)) ||
    length(unique(targets)) != 30L ||
    !identical(sort(unique(s_values)), c(5L, 10L)) ||
    any(table(targets) != 2L) ||
    any(table(s_values) != 30L) ||
    any(table(targets, s_values) != 1L) ||
    !setequal(
      unique(sparsity_tags),
      c("low_sparsity", "medium_sparsity", "high_sparsity")
    ) ||
    any(table(sparsity_tags) != 20L)) {
  stop("The CCLE result does not match the canonical 30-target x 2-s grid.",
       call. = FALSE)
}

output_dir <- file.path(
  project_root,
  "project-progress",
  "combine",
  "figures-and-tables"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Pair calibration methods within GCM and PCM, then show Pearson as baseline.
method_levels <- c("perm_gcm", "asym_gcm", "perm_pcm", "asym_pcm", "naive")
method_labels <- c(
  "asym_gcm" = "GCM (asym)",
  "perm_gcm" = "GCM (perm)",
  "asym_pcm" = "PCM (asym)",
  "perm_pcm" = "PCM (perm)",
  "naive" = "Pearson"
)

## Internal tags refer to mutation FREQUENCY tercile; display labels invert that
## to X-vector SPARSITY: low-frequency mutations correspond to high-sparsity X.
sparsity_levels <- c("low_sparsity", "medium_sparsity", "high_sparsity")
sparsity_labels <- c("low_sparsity"    = "High sparsity",
                     "medium_sparsity" = "Medium sparsity",
                     "high_sparsity"   = "Low sparsity")

# ---------------------------------------------------------------------------
# 1. Long-format rejection table for the 60-setup boxplot
# ---------------------------------------------------------------------------
rej_long <- bind_rows(lapply(results, function(r) {
  data.frame(target         = r$target,
             sparsity_tag   = r$sparsity_tag,
             freq           = r$freq,
             s              = r$s,
             method         = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             stringsAsFactors = FALSE)
})) |>
  mutate(
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    ),
    sparsity_tag = factor(
      sparsity_tag,
      levels = sparsity_levels,
      labels = sparsity_labels[sparsity_levels]
    )
  ) |>
  add_calibration_test_aesthetics()

# ---------------------------------------------------------------------------
# 2. Pick 3 representative setups (one per tercile) at s = 10.
# ---------------------------------------------------------------------------
chosen <- list(
  list(key = "CHIC2_MUT_s10",  tag = "low_sparsity"),
  list(key = "PLK1_MUT_s10",   tag = "medium_sparsity"),
  list(key = "ATP11A_MUT_s10", tag = "high_sparsity")
)
chosen_keys <- vapply(chosen, function(x) x$key, character(1))
stopifnot(all(chosen_keys %in% names(results)))

setup_label <- function(rec) {
  sprintf("%s   (freq = %.3f, s = %d)",
          gsub("_MUT$", "", rec$target), rec$freq, rec$s)
}

setup_order <- vapply(chosen_keys, function(k) setup_label(results[[k]]), character(1))

# ---------------------------------------------------------------------------
# 3. QQ plot for the three chosen setups (Figure-4-style panels: method x setup)
# ---------------------------------------------------------------------------
pv_long <- bind_rows(lapply(chosen_keys, function(k) {
  r  <- results[[k]]
  pv <- as.data.frame(r$pvals)
  pv$rep <- seq_len(nrow(pv))
  pv |>
    pivot_longer(-rep, names_to = "method", values_to = "p_value") |>
    mutate(setup = setup_label(r))
})) |>
  filter(!is.na(p_value)) |>
  mutate(
    setup = factor(setup, levels = setup_order),
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

pv_qq <- pv_long |>
  group_by(setup, method) |>
  arrange(p_value, .by_group = TRUE) |>
  mutate(theoretical = ppoints(n())) |>
  ungroup()

# ---------------------------------------------------------------------------
# 3a. Standalone CHIC2 calibration preview for the introduction. This uses
#     test colour and calibration border type without repeating CHIC2 metadata
#     in a facet strip.
# ---------------------------------------------------------------------------
wilson_interval <- function(successes, repetitions, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  rate <- successes / repetitions
  denominator <- 1 + z^2 / repetitions
  center <- (rate + z^2 / (2 * repetitions)) / denominator
  half_width <- z * sqrt(
    rate * (1 - rate) / repetitions + z^2 / (4 * repetitions^2)
  ) / denominator
  data.frame(
    ci_lower = pmax(0, center - half_width),
    ci_upper = pmin(1, center + half_width)
  )
}

chic2_record <- results[["CHIC2_MUT_s10"]]
if (!is.list(chic2_record) ||
    !all(c("target", "s", "R", "pvals", "rejection_rate") %in%
         names(chic2_record))) {
  stop("Malformed CHIC2 Figure 1 result", call. = FALSE)
}
if (!identical(chic2_record$target, "CHIC2_MUT") ||
    !identical(as.integer(chic2_record$s), 10L) ||
    !identical(as.integer(chic2_record$R), 500L)) {
  stop("The Figure 1 result must be CHIC2_MUT, s=10, R=500", call. = FALSE)
}
chic2_pvals <- as.data.frame(chic2_record$pvals)[method_levels]
if (nrow(chic2_pvals) != chic2_record$R ||
    any(!is.finite(as.matrix(chic2_pvals))) ||
    any(as.matrix(chic2_pvals) < 0 | as.matrix(chic2_pvals) > 1)) {
  stop("Figure 1 contains missing or invalid p-values", call. = FALSE)
}
chic2_repetitions <- colSums(!is.na(chic2_pvals))
chic2_rejections <- colSums(chic2_pvals <= 0.05, na.rm = TRUE)
chic2_summary <- data.frame(
  method = method_levels,
  repetitions = unname(chic2_repetitions),
  rejections = unname(chic2_rejections),
  rejection_rate = unname(chic2_rejections / chic2_repetitions),
  stringsAsFactors = FALSE
) |>
  bind_cols(wilson_interval(chic2_rejections, chic2_repetitions)) |>
  mutate(
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

stopifnot(isTRUE(all.equal(
  unname(chic2_summary$rejection_rate),
  unname(chic2_record$rejection_rate[method_levels])
)))

theme_chic2_intro <- theme_bw(base_size = 9) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, colour = "black"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(0.28, "in"),
    legend.spacing.x = grid::unit(0.01, "in"),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 7.2, margin = margin(r = 2, l = 4)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0.06, "in"),
    plot.margin = margin(3, 2, 8, 2)
  )

chic2_bar_upper <- min(
  1,
  1.20 * max(chic2_summary$ci_upper, 0.05)
)
p_chic2_bar <- ggplot(
  chic2_summary,
  aes(x = method_label, y = rejection_rate)
) +
  geom_col(
    aes(fill = test, linetype = calibration),
    colour = "black",
    linewidth = 0.55,
    width = 0.64
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.13,
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.55
  ) +
  geom_text(
    aes(y = ci_upper, label = sprintf("%.3f", rejection_rate)),
    vjust = -0.5,
    size = 2.5
  ) +
  scale_test_fill() +
  scale_calibration_linetype() +
  scale_y_continuous(
    limits = c(0, chic2_bar_upper),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Rejection rate") +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(
        colour = "black",
        linetype = "solid",
        linewidth = 0.55
      )
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(
        fill = "white",
        colour = "black",
        linewidth = 0.8
      )
    )
  ) +
  theme_chic2_intro +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
chic2_calibration_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-CHIC2-calibration.pdf"
)
ggsave(
  chic2_calibration_path,
  p_chic2_bar,
  width = 5.5,
  height = 3.4
)
cat("Saved: ", chic2_calibration_path, "\n", sep = "")

if (intro_only) {
  quit(save = "no", status = 0L)
}

p_qq <- ggplot(
  pv_qq,
  aes(
    x = theoretical,
    y = p_value,
    colour = test,
    linetype = calibration,
    group = method_label
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_grid(method_label ~ setup) +
  scale_test_colour(guide = "none") +
  scale_calibration_linetype(guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Uniform(0,1) quantile", y = "Observed p-value") +
  theme_figure3(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

qq_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-30targets-pvalqq.pdf"
)
ggsave(qq_path, p_qq, width = 14, height = 11)
cat("Saved: ", qq_path, "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Settings shared by the combined figure's all-setup boxplot.
# ---------------------------------------------------------------------------
n_setups_per_facet <- rej_long |>
  filter(method == method_levels[1]) |>
  count(sparsity_tag)

## Pseudo-log y-axis (asinh-style; handles the few asym-GCM rejection rates
## that are exactly zero) with breaks at log decades, which are evenly spaced
## visually. The 0.05 nominal level is flagged by the dotted black line.
typeI_y_scale <- scale_y_continuous(
  trans  = scales::pseudo_log_trans(sigma = 0.001, base = 10),
  breaks = c(0.001, 0.01, 0.1, 1),
  labels = c("0.001", "0.01", "0.1", "1"),
  limits = c(0, 1),
  expand = expansion(mult = c(0.02, 0.02))
)

# ---------------------------------------------------------------------------
# 5. Rejection rates for the combined figure's representative-setting bars.
# ---------------------------------------------------------------------------
rej_chosen <- bind_rows(lapply(chosen_keys, function(k) {
  r <- results[[k]]
  data.frame(setup          = setup_label(r),
             method         = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             stringsAsFactors = FALSE)
})) |>
  mutate(
    setup = factor(setup, levels = setup_order),
    method = factor(method, levels = method_levels),
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_levels])
    )
  ) |>
  add_calibration_test_aesthetics()

# ---------------------------------------------------------------------------
# 6. Combined two-row figure: box (60 setups) on top, bars (3 setups) below.
#    Shared test-colour and calibration-linetype legends appear at the bottom;
#    x-axis labels are omitted because the legends identify all procedures.
# ---------------------------------------------------------------------------
strip_x <- theme(axis.text.x  = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank())

p_box_c <- ggplot(
  rej_long,
  aes(x = method_label, y = rejection_rate, fill = test)
) +
  geom_boxplot(
    aes(linetype = calibration),
    colour = "black",
    width = 0.65,
    linewidth = 0.55,
    outlier.size = 0.8,
    alpha = 0.85
  ) +
  geom_jitter(aes(colour = test), width = 0.12, height = 0,
              size = 0.55, alpha = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             colour = "black", linewidth = 0.6) +
  facet_wrap(~ sparsity_tag, nrow = 1) +
  scale_test_fill() +
  scale_test_colour(guide = "none") +
  scale_calibration_linetype() +
  typeI_y_scale +
  labs(y = expression("Empirical Type-I rate at " * alpha * " = 0.05")) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(alpha = 1, colour = "black", linetype = "solid")
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(fill = "white", colour = "black", alpha = 1)
    )
  ) +
  theme_figure3(base_size = 13) +
  theme(strip.text       = element_text(face = "bold", size = 11),
        legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.text      = element_text(size = 11),
        legend.key       = element_rect(colour = NA, fill = NA)) +
  strip_x

p_bar_c <- ggplot(
  rej_chosen,
  aes(x = method_label, y = rejection_rate, fill = test)
) +
  geom_col(
    aes(linetype = calibration),
    colour = "black",
    linewidth = 0.55,
    width = 0.7
  ) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             colour = "black", linewidth = 0.6) +
  facet_wrap(~ setup, nrow = 1) +
  scale_test_fill() +
  scale_calibration_linetype() +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(y = expression("Rejection rate at " * alpha * " = 0.05")) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(alpha = 1, colour = "black", linetype = "solid")
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(fill = "white", colour = "black", alpha = 1)
    )
  ) +
  theme_figure3(base_size = 13) +
  theme(strip.text       = element_text(size = 11, face = "bold"),
        legend.box       = "horizontal") +
  strip_x

# Extract the explicitly bottom-positioned guide box.  ggplot2 4 creates
# placeholder guide boxes for several positions, so matching the generic
# "guide-box" can select an empty placeholder instead of the visible legend.
legend_plot <- p_bar_c +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.key.width = grid::unit(0.40, "in"),
    legend.key = element_rect(colour = NA, fill = NA)
  )
shared_legend <- cowplot::get_plot_component(
  legend_plot,
  "guide-box-bottom"
)

p_box_nl <- p_box_c + theme(legend.position = "none")
p_bar_nl <- p_bar_c + theme(legend.position = "none")

p_panels <- cowplot::plot_grid(
  p_box_nl, p_bar_nl,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("(a)", "(b)"),
  label_size = 12,
  label_fontface = "bold",
  rel_heights = c(1, 1)
)

p_combined <- cowplot::plot_grid(
  p_panels, shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.10)
)

combined_path <- file.path(
  output_dir,
  "ccle-typeI-synthY-30targets-typeI-combined.pdf"
)
ggsave(combined_path, p_combined, width = 12, height = 7.8)
cat("Saved: ", combined_path, "\n", sep = "")

# ---------------------------------------------------------------------------
# 7. Console summary
# ---------------------------------------------------------------------------
cat("\n=== Setups per sparsity facet (one row per method) ===\n")
print(n_setups_per_facet)

cat("\n=== Mean / median Type-I rate by method x sparsity (across setups) ===\n")
rej_long |>
  group_by(sparsity_tag, method_label) |>
  summarise(mean   = mean(rejection_rate),
            median = median(rejection_rate),
            min    = min(rejection_rate),
            max    = max(rejection_rate),
            .groups = "drop") |>
  print(n = Inf)

cat("\n=== Three chosen setups for QQ panel ===\n")
for (k in chosen_keys) {
  r <- results[[k]]
  cat(sprintf("  %-18s  tag=%-15s  freq=%.4f  s=%2d\n",
              r$target, r$sparsity_tag, r$freq, r$s))
  cat("    rejection rates: ",
      paste(sprintf("%s=%.3f", names(r$rejection_rate), r$rejection_rate),
            collapse = "  "), "\n", sep = "")
}
