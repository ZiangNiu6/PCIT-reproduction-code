#!/usr/bin/env Rscript
###############################################################################
# Companion figure to pcm_additive.pdf: calibration under three
# (null, nuisance) regimes.  The 3 x 3 figure is unaffected; this is additional.
#
# Rows:
#   1. Random forest, MI null      (X independent of Z; setting1_xindep, n=100)
#   2. Random forest, CI null      (Heteroscedastic; setting2, n=600)
#   3. Penalised spline, CI null   (Heteroscedastic; setting2, n=200)
#
# Row 1 is where permutation calibration is EXACT in finite samples, and it
# uses the learner a practitioner would reach for, so the guarantee cannot be
# attributed to a convenient nuisance.
#
# Rows 2 and 3 are the same design under two nuisance regressions.  The forest
# is given three times the sample size of the spline, so the comparison is
# conservative: it handicaps the spline.  The forest's asymptotic rejection is
# nevertheless 0.122 at n = 200 and 0.124 at n = 600 -- flat, because with a
# one-dimensional Z ranger must take mtry = 1, the forest degenerates to bagged
# unpruned trees, and its nuisance error does not vanish.  The spline with a
# FIXED smoothing parameter is by contrast a linear smoother, and so lies
# inside the class the permutation-stability theory covers.
#
# One column: Type-I error at alpha = 0.05 with 95% Wilson intervals.  Power is
# not shown because the spline cell is a null-only single-cell specifier, and a
# p-value QQ panel is not shown because the adaptive permutation scheme stops as
# soon as p_AV <= alpha, so permutation p-values pile up at 0.05 and the lower
# tail is censored by construction.  Rejection at alpha is exactly the quantity
# the adaptive scheme estimates correctly.
#
# Usage: Rscript plot-pcm-additive-nuisance.R [--validate-only]
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(cowplot)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_argument))
script_directory <- dirname(script_file)
project_root <- normalizePath(file.path(script_directory, "..", "..", ".."))
source(file.path(project_root, "reproduction-code", "plotting-style.R"))
validate_only <- "--validate-only" %in% commandArgs(trailingOnly = TRUE)

results_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(results_root)) {
  results_root <- file.path(project_root, "results")
} else if (!grepl("^/", results_root)) {
  results_root <- file.path(project_root, results_root)
}
sim_dir <- file.path(results_root, "simulation-results")
diagnostic_directory <- file.path(sim_dir, "pcm-additive-matched-dcrt")

ALPHA <- 0.05

## ---- locate the three cells ------------------------------------------------
read_cell <- function(file, run_method, schema, suffix, dgp, n, label) {
  path <- file.path(sim_dir, dirname(file), basename(file))
  if (!file.exists(path)) stop("Result not found: ", path)
  o <- readRDS(path)
  res <- o$results
  # grid_id -> (dgp, n) mapping comes from the specifier's own grid
  pg <- as.data.frame(o$parameter_grid %||% NULL)
  if (is.null(pg) || !nrow(pg)) {
    grid_file <- file.path(dirname(script_directory), "sim-spec-formal",
                           paste0(sub("_results\\.rds$", "", basename(file)),
                                  "-grid.R"))
    env <- new.env(); sys.source(grid_file, envir = env)
    pg <- as.data.frame(env$parameter_grid)
    if (is.null(pg$grid_id)) pg$grid_id <- seq_len(nrow(pg))
  }
  gid <- pg$grid_id[pg$dgp == dgp & pg$n == n & pg$gamma_signal == 0]
  if (length(gid) != 1L) stop("Cell not uniquely identified for ", label)
  sub <- res[
    res$grid_id == gid & as.character(res$method) == run_method,
    ,
    drop = FALSE
  ]
  if (nrow(sub) == 0L) {
    stop("No `", run_method, "` outputs found for ", label)
  }
  stopifnot(all(vapply(sub$output, function(x)
    identical(x$schema_version, schema), logical(1L))))
  P <- do.call(rbind, lapply(sub$output, function(x) x$p_values))
  methods <- paste(c("perm_gcm", "asym_gcm", "perm_pcm", "asym_pcm"),
                   suffix, sep = "_")
  stopifnot(all(methods %in% colnames(P)))
  as.data.frame(P[, methods, drop = FALSE]) |>
    setNames(c("perm_gcm", "asym_gcm", "perm_pcm", "asym_pcm")) |>
    tidyr::pivot_longer(everything(), names_to = "method",
                        values_to = "p.value") |>
    dplyr::mutate(row_label = label, R = nrow(P))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

ROWS <- c(
  "Forest, MI null (n = 100)",
  "Forest, CI null (n = 600)",
  "Spline, CI null (n = 200)"
)
pv <- dplyr::bind_rows(
  read_cell("pcm-additive-matched-dcrt/pcm-additive_results.rds",
            "matched_gcm_pcm_rf",
            "pcm_additive_matched", "rf", "setting1_xindep", 100L, ROWS[1]),
  read_cell("pcm-additive-matched-dcrt/pcm-additive_results.rds",
            "matched_gcm_pcm_rf",
            "pcm_additive_matched", "rf", "setting2", 600L, ROWS[2]),
  read_cell("pcm-additive-gam/pcm-additive-gam_results.rds",
            "matched_gcm_pcm_gam",
            "pcm_additive_matched_gam", "gam", "setting2", 200L, ROWS[3])
) |>
  dplyr::mutate(row_label = factor(row_label, levels = ROWS)) |>
  add_calibration_test_aesthetics()

if (any(is.na(pv$p.value)) || any(pv$p.value < 0 | pv$p.value > 1)) {
  stop("Non-probability p-values encountered.")
}
if (validate_only) {
  cat(sprintf("Validated %d p-values across %d rows.\n", nrow(pv),
              dplyr::n_distinct(pv$row_label)))
  quit(save = "no", status = 0L)
}

method_display <- c(perm_gcm = "Permutation GCM", asym_gcm = "Asymptotic GCM",
                    perm_pcm = "Permutation PCM", asym_pcm = "Asymptotic PCM")
pv$method_label <- factor(unname(method_display[pv$method]),
                          levels = unname(method_display))

## ---- Type-I error ----------------------------------------------------------
wilson <- function(k, N, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2); ph <- k / N; d <- 1 + z^2 / N
  ctr <- (ph + z^2 / (2 * N)) / d
  hw <- (z / d) * sqrt(ph * (1 - ph) / N + z^2 / (4 * N^2))
  c(lower = max(0, ctr - hw), upper = min(1, ctr + hw))
}
typeI <- pv |>
  dplyr::group_by(.data$row_label, .data$method_label, .data$calibration,
                  .data$test) |>
  dplyr::summarise(k = sum(.data$p.value <= ALPHA), N = dplyr::n(),
                   .groups = "drop") |>
  dplyr::mutate(rate = .data$k / .data$N,
                lower = mapply(function(a, b) wilson(a, b)[["lower"]], k, N),
                upper = mapply(function(a, b) wilson(a, b)[["upper"]], k, N))
typeI_upper <- max(0.15, ceiling(max(typeI$upper) * 1.15 / 0.05) * 0.05)

make_typeI <- function(rl, show_x = FALSE) {
  ggplot(dplyr::filter(typeI, .data$row_label == rl),
         aes(x = .data$method_label, y = .data$rate,
             fill = .data$calibration, linetype = .data$test)) +
    geom_hline(yintercept = ALPHA, linetype = "dashed", colour = "black",
               linewidth = 0.5) +
    geom_col(colour = "black", linewidth = 0.4, width = 0.68) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0,
                  linewidth = 0.45, linetype = "solid") +
    geom_text(aes(y = .data$upper, label = sprintf("%.3f", .data$rate)),
              vjust = -0.6, size = 3) +
    scale_calibration_fill(guide = "none") +
    scale_linetype_manual(values = TEST_LINETYPES, guide = "none") +
    scale_y_continuous(limits = c(0, typeI_upper),
                       breaks = seq(0, typeI_upper, by = 0.05),
                       labels = scales::label_percent(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Type-I error") +
    theme_figure3() +
    theme(legend.position = "none",
          axis.text.x = if (show_x) element_text(angle = 20, hjust = 1)
                        else element_blank())
}

text_panel <- function(label, font_size = 11.5, rotation = 0) {
  patchwork::wrap_elements(full = grid::textGrob(
    label, rot = rotation,
    gp = grid::gpar(fontsize = font_size, fontface = "bold")))
}

legend_src <- ggplot(typeI, aes(x = .data$method_label, y = .data$rate,
                                fill = .data$calibration)) +
  geom_col() +
  scale_calibration_fill(guide = guide_legend(nrow = 1)) +
  theme_figure3() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.direction = "horizontal")
shared_legend <- cowplot::get_plot_component(legend_src, "guide-box-bottom")

joint <- patchwork::wrap_plots(
  A = patchwork::plot_spacer(),
  B = text_panel("Type-I error (alpha = 0.05)"),
  C = text_panel(ROWS[1], 9.5, rotation = 90),
  D = make_typeI(ROWS[1]) + labs(tag = "a"),
  E = text_panel(ROWS[2], 9.5, rotation = 90),
  F = make_typeI(ROWS[2]) + labs(tag = "b"),
  G = text_panel(ROWS[3], 9.5, rotation = 90),
  H = make_typeI(ROWS[3], show_x = TRUE) + labs(tag = "c"),
  I = patchwork::wrap_elements(full = shared_legend),
  design = "AB\nCD\nEF\nGH\nII",
  widths = c(0.16, 1), heights = c(0.10, 1, 1, 1.30, 0.14)
)

dir.create(diagnostic_directory, recursive = TRUE, showWarnings = FALSE)
out <- file.path(diagnostic_directory, "pcm_additive_nuisance.pdf")
ggsave(out, joint, height = 8.4, width = 5.9, device = cairo_pdf)
cat("Plot saved to ", out, "\n", sep = "")
