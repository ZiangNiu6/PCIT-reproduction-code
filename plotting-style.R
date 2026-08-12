## Shared visual encoding for GCM/PCM calibration comparisons.
##
## Calibration is encoded by colour, while test type is encoded by line type.
## Simulation and CCLE plotting scripts source this file so their legends do
## not drift apart.

CALIBRATION_LEVELS <- c("Permutation", "Asymptotic", "Pearson")
CALIBRATION_COLORS <- c(
  "Permutation" = "#E41A1C",
  "Asymptotic" = "#377EB8",
  "Pearson" = "#999999"
)

TEST_LEVELS <- c("GCM", "PCM", "Baseline")
TEST_LINETYPES <- c(
  "GCM" = "solid",
  "PCM" = "dashed",
  "Baseline" = "solid"
)

add_calibration_test_aesthetics <- function(data, method_column = "method") {
  method_code <- as.character(data[[method_column]])
  calibration <- dplyr::case_when(
    grepl("^perm_", method_code) ~ "Permutation",
    grepl("^asym_", method_code) ~ "Asymptotic",
    method_code %in% c("naive", "pearson") ~ "Pearson",
    TRUE ~ NA_character_
  )
  test <- dplyr::case_when(
    grepl("_gcm(_|$)", method_code) ~ "GCM",
    grepl("_pcm(_|$)", method_code) ~ "PCM",
    method_code %in% c("naive", "pearson") ~ "Baseline",
    TRUE ~ NA_character_
  )

  if (anyNA(calibration) || anyNA(test)) {
    bad_methods <- unique(method_code[is.na(calibration) | is.na(test)])
    stop(
      "Unknown method code(s) in shared plotting style: ",
      paste(bad_methods, collapse = ", "),
      call. = FALSE
    )
  }

  data$calibration <- factor(calibration, levels = CALIBRATION_LEVELS)
  data$test <- factor(test, levels = TEST_LEVELS)
  data
}

scale_calibration_colour <- function(include_pearson = FALSE,
                                     guide = "legend",
                                     name = "Calibration") {
  breaks <- if (include_pearson) {
    CALIBRATION_LEVELS
  } else {
    CALIBRATION_LEVELS[1:2]
  }
  ggplot2::scale_colour_manual(
    name = name,
    values = CALIBRATION_COLORS,
    breaks = breaks,
    guide = guide
  )
}

scale_calibration_fill <- function(include_pearson = FALSE, guide = "legend") {
  breaks <- if (include_pearson) {
    CALIBRATION_LEVELS
  } else {
    CALIBRATION_LEVELS[1:2]
  }
  ggplot2::scale_fill_manual(
    name = "Calibration",
    values = CALIBRATION_COLORS,
    breaks = breaks,
    guide = guide
  )
}

scale_test_linetype <- function(guide = "legend") {
  ggplot2::scale_linetype_manual(
    name = "Test",
    values = TEST_LINETYPES,
    breaks = TEST_LEVELS[1:2],
    guide = guide
  )
}

theme_figure3 <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 11),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(0.22, "in"),
      legend.spacing.x = grid::unit(0.03, "in"),
      legend.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_text(size = 10.5),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3),
      plot.tag = ggplot2::element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.title = ggplot2::element_text(size = 11, hjust = 0.5)
    )
}
