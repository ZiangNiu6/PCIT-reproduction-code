#!/usr/bin/env Rscript

# Plot the finite-sample normal-reference failure using saved sparse-NB runs.
# Run from any directory with: Rscript --vanilla plot-nb-normal-reference.R

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
result_dir <- file.path(project_root, "results", "simulation-results",
                        "nb-sparse-ci-search")

read_rows <- function(run_name) {
  path <- file.path(result_dir, run_name, "nb-sparse-ci-search-all.rds")
  readRDS(path)$replicates
}

baseline <- read_rows("baseline")
moderate <- read_rows("confirm-theta0p5-beta1-delta0p5-normal-reference")
small <- rbind(
  read_rows("confirm-theta0p1-beta1-delta0p5-poisson-only"),
  read_rows("confirm-theta0p1-beta1-delta0p5-poisson-only-independent")
)

summarize_method <- function(rows, setting, method, rejected) {
  n <- nrow(rows)
  k <- sum(rejected, na.rm = TRUE)
  ci <- stats::binom.test(k, n)$conf.int[1:2]
  data.frame(
    setting = setting, method = method, n = n, rejected = k,
    rejection_rate = k / n, ci_lower = ci[1], ci_upper = ci[2]
  )
}

make_summary <- function(rows, setting) {
  p_normal <- 2 * stats::pnorm(-abs(rows$empirical_statistic))
  rbind(
    summarize_method(rows, setting, "Barry permutation", rows$reject_barry),
    summarize_method(rows, setting, "Our permutation", rows$reject_ours),
    summarize_method(rows, setting, "Our normal reference", p_normal <= 0.05)
  )
}

summary <- rbind(
  make_summary(baseline, "Baseline: size 0.5, effects 0.5/0.5"),
  make_summary(moderate, "Matched: size 0.5, effects 1/0.5"),
  make_summary(small, "Small size: 0.1, effects 1/0.5")
)
csv_path <- file.path(result_dir, "normal-reference-confirmations.csv")
utils::write.csv(summary, csv_path, row.names = FALSE)

png_path <- file.path(result_dir, "normal-reference-diagnostic.png")
grDevices::png(png_path, width = 1800, height = 760, res = 160)
op <- graphics::par(mfrow = c(1, 2), mar = c(4.8, 4.8, 2.5, 1.2),
                    mgp = c(2.8, 0.8, 0), cex = 0.9)

x <- seq(-4.5, 3.5, length.out = 500L)
graphics::plot(x, stats::pnorm(x), type = "l", lwd = 2, lty = 2,
               col = "#4D4D4D", ylim = c(0, 1),
               xlab = "Observed empirical-studentized score",
               ylab = "Cumulative probability",
               main = "Null score distribution (effects 1/0.5)")
graphics::abline(v = c(-1.96, 1.96), col = "#AAAAAA", lty = 3)
graphics::lines(stats::ecdf(moderate$empirical_statistic),
                col = "#277DA1", lwd = 2)
graphics::lines(stats::ecdf(small$empirical_statistic),
                col = "#E76F51", lwd = 2)
graphics::legend("bottomright",
                 legend = c("Standard normal", "Size 0.5", "Size 0.1",
                            "Normal cutoffs"),
                 col = c("#4D4D4D", "#277DA1", "#E76F51", "#AAAAAA"),
                 lty = c(2, 1, 1, 3), lwd = c(2, 2, 2, 1), bty = "n")

method_order <- c("Barry permutation", "Our permutation",
                  "Our normal reference")
setting_order <- unique(summary$setting)[2:3]
rates <- sapply(setting_order, function(s) {
  summary$rejection_rate[summary$setting == s]
})
lower <- sapply(setting_order, function(s) {
  summary$ci_lower[summary$setting == s]
})
upper <- sapply(setting_order, function(s) {
  summary$ci_upper[summary$setting == s]
})
stopifnot(identical(summary$method[seq_along(method_order)], method_order))
colors <- c("#4C78A8", "#54A24B", "#E76F51")
positions <- graphics::barplot(
  rates, beside = TRUE, col = colors, border = NA, ylim = c(0, 0.33),
  names.arg = c("NB size 0.5", "NB size 0.1"),
  ylab = "Null rejection rate",
  main = "Type-I rejection (effects 1/0.5)"
)
graphics::abline(h = 0.05, col = "#333333", lty = 3, lwd = 1.5)
graphics::arrows(positions, lower, positions, upper,
                 angle = 90, code = 3, length = 0.035, col = "#333333")
graphics::legend("topleft", legend = method_order, fill = colors, bty = "n")
graphics::par(op)
grDevices::dev.off()
cat("Saved ", csv_path, " and ", png_path, "\n", sep = "")
