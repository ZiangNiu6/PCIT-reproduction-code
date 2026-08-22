#!/usr/bin/env Rscript
## Publication figure for the Frangieh low-response-expression self-target
## sensitivity panel. The input is the full B=4,999 exact-one-guide pair CSV.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path")
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)))
project_root <- normalizePath(file.path(script_dir, "..", ".."))

usage <- function() {
  paste(
    "Usage: Rscript --vanilla plot_frangieh_low_expression_pvalue_comparison.R",
    "  [--pairs FILE] [--output-dir DIR] [--prefix NAME]",
    sep = "\n"
  )
}

parse_cli <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    if (identical(token, "--help")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
      key <- gsub("-", "_", bits[[1L]], fixed = TRUE)
      out[[key]] <- paste(bits[-1L], collapse = "=")
      i <- i + 1L
      next
    }
    if (i == length(x)) stop("Missing value after ", token)
    key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(cli$help)) {
  cat(usage(), "\n")
  quit(save = "no", status = 0L)
}

pair_file <- if (!is.null(cli$pairs)) {
  cli$pairs
} else {
  file.path(
    project_root, "results", "frangieh-2021",
    "exact-one-guide-lowest-expression-50-B4999",
    "frangieh-low-expression-self-target-pairs.csv"
  )
}
if (!file.exists(pair_file)) stop("Pair CSV not found: ", pair_file)
pair_file <- normalizePath(pair_file)

output_dir <- if (!is.null(cli$output_dir)) {
  cli$output_dir
} else {
  file.path(project_root, "results", "frangieh-2021", "comparison")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
prefix <- if (!is.null(cli$prefix)) {
  cli$prefix
} else {
  "frangieh-low-expression-b4999-pvalue-comparison"
}
if (!nzchar(prefix) || grepl("[/\\\\]", prefix)) {
  stop("--prefix must be a nonempty filename stem")
}

EXPECTED_N <- 50L
EXPECTED_B <- 4999L
ALPHA <- 0.05
required_columns <- c(
  "gene", "information_stat", "empirical_stat",
  "hits_abs_barry", "hits_abs_empirical",
  "p_two_sided_barry", "p_two_sided_empirical"
)
pairs <- utils::read.csv(pair_file, stringsAsFactors = FALSE,
                         check.names = FALSE)
missing_columns <- setdiff(required_columns, names(pairs))
if (length(missing_columns)) {
  stop("Pair CSV lacks: ", paste(missing_columns, collapse = ", "))
}
if (nrow(pairs) != EXPECTED_N || anyDuplicated(pairs$gene)) {
  stop("Expected exactly 50 unique self-target gene pairs")
}
p_columns <- c("p_two_sided_barry", "p_two_sided_empirical")
if (any(!is.finite(as.matrix(pairs[p_columns]))) ||
    any(as.matrix(pairs[p_columns]) <= 0) ||
    any(as.matrix(pairs[p_columns]) > 1)) {
  stop("Two-sided p-values must be finite and in (0, 1]")
}
expected_barry <- (1 + pairs$hits_abs_barry) / (EXPECTED_B + 1)
expected_empirical <- (1 + pairs$hits_abs_empirical) / (EXPECTED_B + 1)
if (!isTRUE(all.equal(pairs$p_two_sided_barry, expected_barry,
                      tolerance = 1e-14)) ||
    !isTRUE(all.equal(pairs$p_two_sided_empirical, expected_empirical,
                      tolerance = 1e-14))) {
  stop("P-values do not reconstruct from B=4,999 absolute-tail hit counts")
}
if (any(sign(pairs$information_stat) != sign(pairs$empirical_stat))) {
  stop("Barry and empirical statistics must share the adjusted-score sign")
}

pairs$barry_detected <- pairs$p_two_sided_barry <= ALPHA
pairs$empirical_detected <- pairs$p_two_sided_empirical <= ALPHA
pairs$discordant <- xor(pairs$barry_detected, pairs$empirical_detected)
pairs$adjusted_sign <- factor(
  ifelse(pairs$information_stat < 0, "Downregulated", "Upregulated"),
  levels = c("Downregulated", "Upregulated")
)
pairs$x <- -log10(pairs$p_two_sided_barry)
pairs$y <- -log10(pairs$p_two_sided_empirical)

n_barry <- sum(pairs$barry_detected)
n_empirical <- sum(pairs$empirical_detected)
n_both <- sum(pairs$barry_detected & pairs$empirical_detected)
n_barry_only <- sum(pairs$barry_detected & !pairs$empirical_detected)
n_empirical_only <- sum(!pairs$barry_detected & pairs$empirical_detected)
if (!identical(
  c(n_barry, n_empirical, n_both, n_barry_only, n_empirical_only),
  c(7L, 14L, 6L, 1L, 8L)
)) {
  stop("Unexpected alpha=.05 detection counts for the audited B=4,999 panel")
}

threshold <- -log10(ALPHA)
axis_max <- max(4, ceiling(10 * max(c(pairs$x, pairs$y))) / 10 + 0.15)
annotation_data <- data.frame(
  x = c(0.16, axis_max - 0.16),
  y = c(axis_max - 0.16, 0.16),
  label = c("Empirical p smaller", "Barry p smaller"),
  hjust = c(0, 1)
)

subtitle <- paste0(
  "Nominal α = .05: Barry ", n_barry, "/", EXPECTED_N,
  "; empirical ", n_empirical, "/", EXPECTED_N,
  "  (", n_both, " both, ", n_barry_only, " Barry-only, ",
  n_empirical_only, " empirical-only)"
)
caption_text <- paste(
  strwrap(
    paste0(
      "Full bottom-50 response-expression panel; B = 4,999 paired label ",
      "permutations. Circled and labeled genes are discordant at nominal ",
      "alpha = .05. Exact-one-guide empirical sensitivity analysis."
    ),
    width = 104
  ),
  collapse = "\n"
)

plot_object <- ggplot(pairs, aes(x = x, y = y)) +
  annotate(
    "rect", xmin = 0, xmax = threshold, ymin = threshold, ymax = Inf,
    fill = "#6A51A3", alpha = 0.055
  ) +
  annotate(
    "rect", xmin = threshold, xmax = Inf, ymin = 0, ymax = threshold,
    fill = "#E6550D", alpha = 0.045
  ) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.55,
              colour = "grey65") +
  geom_vline(xintercept = threshold, linetype = "22", linewidth = 0.55,
             colour = "grey38") +
  geom_hline(yintercept = threshold, linetype = "22", linewidth = 0.55,
             colour = "grey38") +
  geom_point(
    aes(colour = adjusted_sign), size = 2.8, alpha = 0.88
  ) +
  geom_point(
    data = pairs[pairs$discordant, , drop = FALSE],
    shape = 21, size = 4.25, stroke = 0.72, colour = "black", fill = NA,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = pairs[pairs$discordant, , drop = FALSE],
    aes(label = gene),
    seed = 20260822, size = 3.1, fontface = "plain", colour = "black",
    box.padding = 0.50, point.padding = 0.34, min.segment.length = 0,
    segment.size = 0.35, segment.colour = "grey45", max.overlaps = Inf,
    force = 1.8, force_pull = 0.8, max.time = 3, max.iter = 20000,
    show.legend = FALSE
  ) +
  geom_text(
    data = annotation_data,
    aes(x = x, y = y, label = label, hjust = hjust),
    inherit.aes = FALSE, size = 3.05, colour = "grey38"
  ) +
  annotate(
    "text", x = threshold + 0.06, y = axis_max - 0.05,
    label = "α = .05", hjust = 0, vjust = 1, size = 3.05,
    colour = "grey32"
  ) +
  scale_colour_manual(
    name = "Adjusted score sign",
    values = c("Downregulated" = "#0072B2", "Upregulated" = "#D55E00"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, floor(axis_max), by = 1),
    limits = c(0, axis_max), expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, floor(axis_max), by = 1),
    limits = c(0, axis_max), expand = expansion(mult = 0)
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Low-response-expression self-target controls",
    subtitle = subtitle,
    x = expression(paste("Barry et al.:  ", -log[10], "(two-sided ", italic(p), ")")),
    y = expression(paste("Empirical studentization:  ", -log[10], "(two-sided ", italic(p), ")")),
    caption = caption_text
  ) +
  theme_classic(base_size = 11.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 14.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10.7, colour = "grey25",
                                 margin = margin(b = 11)),
    plot.caption = element_text(size = 8.8, colour = "grey35", hjust = 0,
                                margin = margin(t = 10)),
    axis.title = element_text(size = 11.3),
    axis.text = element_text(colour = "grey20"),
    axis.line = element_line(linewidth = 0.45, colour = "grey25"),
    axis.ticks = element_line(linewidth = 0.4, colour = "grey30"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(1.2, "lines"),
    plot.margin = margin(12, 24, 10, 12)
  )

png_file <- file.path(output_dir, paste0(prefix, ".png"))
pdf_file <- file.path(output_dir, paste0(prefix, ".pdf"))
ggsave(
  png_file, plot_object, width = 7.8, height = 7.15, units = "in",
  dpi = 320, bg = "white"
)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggsave(
  pdf_file, plot_object, width = 7.8, height = 7.15, units = "in",
  device = pdf_device, bg = "white"
)

cat("Frangieh low-expression p-value comparison\n")
cat("  input:       ", pair_file, "\n", sep = "")
cat("  pairs:       ", nrow(pairs), "\n", sep = "")
cat("  discordant:  ", sum(pairs$discordant), "\n", sep = "")
cat("  detections:  Barry ", n_barry, "; empirical ", n_empirical,
    " at alpha=.05\n", sep = "")
cat("  PNG:         ", png_file, "\n", sep = "")
cat("  PDF:         ", pdf_file, "\n", sep = "")
