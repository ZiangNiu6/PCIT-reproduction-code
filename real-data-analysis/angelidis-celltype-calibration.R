#!/usr/bin/env Rscript

# Conditional permutation-orbit calibration for the Angelidis et al. (2019)
# lung-aging cell-composition data (GEO GSE124872).
#
# The 15 mouse-level count vectors are held fixed.  We enumerate all
# choose(15, 7) = 6,435 assignments of seven "old" labels and compare:
#   * the exact, empirically studentized beta-binomial proportional-variance
#     quasi-score test (PCIT);
#   * the same proportional-variance score with a standard-normal reference;
#   * DCATS-style beta-binomial LRTs with a chi-square reference;
#
# This is a conditional random-label calibration diagnostic, not a
# repeated-sampling simulation and not design-based causal inference: age was
# observational in the original study.

parse_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1L]])
}

as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

script_path <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
)
script_dir <- if (!is.na(script_path)) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, "..", ".."))
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- 1L
output_dir <- parse_arg(
  "out-dir",
  file.path(project_root, "data", "angelidis")
)

cfg <- list(
  metadata = parse_arg(
    "metadata",
    "/tmp/GSE124872_Angelidis_2018_metadata.csv.gz"
  ),
  out_dir = output_dir,
  figure_dir = parse_arg(
    "figure-dir",
    file.path(output_dir, "figures")
  ),
  include_33_sensitivity = as_flag(parse_arg("include-33-sensitivity", "true")),
  validate_aod = as_flag(parse_arg("validate-aod", "true")),
  aod_lib = parse_arg("aod-lib", "/private/tmp/aod-lib"),
  cores = as.integer(parse_arg(
    "cores",
    as.character(min(6L, max(1L, detected_cores - 1L)))
  )),
  tolerance = as.numeric(parse_arg("tolerance", "1e-8"))
)

if (!file.exists(cfg$metadata)) {
  stop(
    "Angelidis metadata not found: ", cfg$metadata, "\n",
    "Download GSE124872_Angelidis_2018_metadata.csv.gz from GEO and pass ",
    "--metadata=/absolute/path/to/file.",
    call. = FALSE
  )
}
if (!is.finite(cfg$tolerance) || cfg$tolerance <= 0) {
  stop("--tolerance must be a positive number", call. = FALSE)
}
if (!is.finite(cfg$cores) || cfg$cores < 1L) {
  stop("--cores must be a positive integer", call. = FALSE)
}

make_assignments <- function(n, old) {
  selected <- utils::combn(n, old)
  assignments <- matrix(0L, nrow = n, ncol = ncol(selected))
  assignments[cbind(
    as.vector(selected),
    rep(seq_len(ncol(selected)), each = old)
  )] <- 1L
  assignments
}

read_angelidis <- function(path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path) else path
  metadata <- utils::read.csv(
    connection,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  required <- c("identifier", "grouping", "celltype")
  if (!all(required %in% names(metadata))) {
    stop(
      "Metadata must contain: ", paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  missing_celltype <- is.na(metadata$celltype) | !nzchar(metadata$celltype)
  if (sum(missing_celltype) != 2L) {
    stop("Expected exactly two unannotated cells in the GEO metadata")
  }

  mouse_age <- unique(metadata[c("identifier", "grouping")])
  mouse_age <- mouse_age[
    !is.na(mouse_age$identifier) & nzchar(mouse_age$identifier) &
      !is.na(mouse_age$grouping) & nzchar(mouse_age$grouping),
    , drop = FALSE
  ]
  if (anyDuplicated(mouse_age$identifier)) {
    stop("At least one mouse maps to multiple age groups", call. = FALSE)
  }
  mouse_age <- mouse_age[order(mouse_age$identifier), , drop = FALSE]
  rownames(mouse_age) <- mouse_age$identifier
  if (!identical(sort(table(mouse_age$grouping)), c(`24m` = 7L, `3m` = 8L))) {
    # Compare values rather than names because sort(table(.)) reorders names.
    if (!identical(as.integer(table(mouse_age$grouping)[c("24m", "3m")]), c(7L, 8L))) {
      stop("Expected exactly seven 24m and eight 3m mice", call. = FALSE)
    }
  }
  if (nrow(mouse_age) != 15L) stop("Expected exactly 15 mice", call. = FALSE)

  annotated <- metadata[
    !is.na(metadata$identifier) & nzchar(metadata$identifier) &
      !is.na(metadata$celltype) & nzchar(metadata$celltype),
    , drop = FALSE
  ]
  all_types <- sort(unique(annotated$celltype))
  annotated$identifier <- factor(
    annotated$identifier,
    levels = mouse_age$identifier
  )
  annotated$celltype <- factor(annotated$celltype, levels = all_types)
  counts_33 <- as.matrix(stats::xtabs(~ identifier + celltype, data = annotated))
  storage.mode(counts_33) <- "integer"

  excluded_from_paper_30 <- c(
    "Gamma-Delta_T_cells", "low_quality_cells", "red_blood_cells"
  )
  if (!all(excluded_from_paper_30 %in% colnames(counts_33))) {
    stop("Could not identify the three labels excluded from the paper's 30 types")
  }
  primary_types <- setdiff(colnames(counts_33), excluded_from_paper_30)
  counts_30 <- counts_33[, primary_types, drop = FALSE]

  if (ncol(counts_33) != 33L || ncol(counts_30) != 30L) {
    stop("Expected 33 annotated metadata labels and 30 paper-defined types")
  }
  if (sum(counts_33) != 14811L || sum(counts_30) != 14124L) {
    stop("Unexpected annotated-cell totals; the GEO metadata may have changed")
  }
  if (any(rowSums(counts_30) <= 0L)) stop("A mouse has no retained cells")

  list(
    metadata_rows = nrow(metadata),
    mouse_age = mouse_age,
    observed_old = as.integer(mouse_age$grouping == "24m"),
    counts_30 = counts_30,
    counts_33 = counts_33,
    excluded_from_paper_30 = excluded_from_paper_30
  )
}

# Beta-binomial parameterization used by aod/DCATS:
#   alpha = mu * (1 - phi) / phi,
#   beta  = (1 - mu) * (1 - phi) / phi.
bb_loglik <- function(eta, y, m, phi) {
  mu <- stats::plogis(pmin(pmax(eta, -30), 30))
  phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
  tau <- (1 - phi) / phi
  shape1 <- pmax(mu * tau, 1e-12)
  shape2 <- pmax((1 - mu) * tau, 1e-12)
  sum(
    lchoose(m, y) +
      lbeta(y + shape1, m - y + shape2) -
      lbeta(shape1, shape2)
  )
}

fit_eta <- function(y, m, phi, tol = 1e-10) {
  fit <- stats::optimize(
    function(eta) bb_loglik(eta, y, m, phi),
    interval = c(-25, 25),
    maximum = TRUE,
    tol = tol
  )
  list(eta = fit$maximum, loglik = fit$objective)
}

fit_cell_null <- function(y, m, tol = 1e-8) {
  profile_loglik <- function(logit_phi) {
    phi <- stats::plogis(logit_phi)
    fit_eta(y, m, phi, tol = tol / 10)$loglik
  }
  fit_phi <- stats::optimize(
    profile_loglik,
    interval = c(stats::qlogis(1e-8), stats::qlogis(0.8)),
    maximum = TRUE,
    tol = tol
  )
  phi <- stats::plogis(fit_phi$maximum)
  eta_fit <- fit_eta(y, m, phi, tol = tol / 10)
  list(
    eta = eta_fit$eta,
    phi = phi,
    loglik = eta_fit$loglik,
    phi_at_boundary = phi <= 1.01e-8 || phi >= 0.799
  )
}

fit_pooled_null_phi <- function(counts, totals, tol = 1e-8) {
  profile_loglik <- function(logit_phi) {
    phi <- stats::plogis(logit_phi)
    sum(vapply(seq_len(ncol(counts)), function(j) {
      fit_eta(counts[, j], totals, phi, tol = tol / 10)$loglik
    }, numeric(1L)))
  }
  fit_phi <- stats::optimize(
    profile_loglik,
    interval = c(stats::qlogis(1e-6), stats::qlogis(0.5)),
    maximum = TRUE,
    tol = tol
  )
  phi <- stats::plogis(fit_phi$maximum)
  list(
    phi = phi,
    loglik = fit_phi$objective,
    phi_at_boundary = phi <= 1.01e-6 || phi >= 0.499
  )
}

# For Var(Y_i | M_i) = M_i mu (1-mu) {1 + (M_i-1) phi}, the
# logit-mean quasi-score is proportional to
#   (Y_i - M_i mu) / {1 + (M_i-1) phi}.
# We estimate mu under the null by solving the corresponding intercept
# estimating equation, conditional on the null-only cell-specific phi.
proportional_variance_score <- function(y, m, phi) {
  phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
  variance_multiplier <- 1 + (m - 1) * phi
  mu <- sum(y / variance_multiplier) / sum(m / variance_multiplier)
  if (!is.finite(mu) || mu <= 0 || mu >= 1) {
    stop("Invalid proportional-variance null mean")
  }
  list(
    score = (y - m * mu) / variance_multiplier,
    mu = mu,
    variance_multiplier = variance_multiplier
  )
}

studentized_score_statistics <- function(score, assignments) {
  propensity <- rowMeans(assignments)
  centered_x <- sweep(assignments, 1L, propensity, FUN = "-")
  contributions <- sweep(centered_x, 1L, score, FUN = "*")
  totals <- colSums(contributions)
  centered_ss <- colSums(contributions^2) - totals^2 / nrow(contributions)
  denominator <- sqrt(pmax(centered_ss, 0))
  statistic <- totals / denominator
  statistic[!is.finite(statistic) | denominator <= 0] <- NA_real_
  statistic
}

inclusive_rank_p <- function(statistic, decreasing_extremeness = TRUE) {
  ok <- is.finite(statistic)
  out <- rep(NA_real_, length(statistic))
  if (decreasing_extremeness) {
    out[ok] <- rank(-statistic[ok], ties.method = "max") / sum(ok)
  } else {
    out[ok] <- rank(statistic[ok], ties.method = "max") / sum(ok)
  }
  out
}

proportional_variance_score_p_matrices <- function(
    counts,
    assignments,
    tol = 1e-8) {
  totals <- rowSums(counts)
  null_fits <- lapply(seq_len(ncol(counts)), function(j) {
    fit_cell_null(counts[, j], totals, tol = tol)
  })
  names(null_fits) <- colnames(counts)

  quasi_fits <- lapply(seq_len(ncol(counts)), function(j) {
    fit <- null_fits[[j]]
    proportional_variance_score(counts[, j], totals, fit$phi)
  })
  names(quasi_fits) <- colnames(counts)

  statistics <- vapply(seq_len(ncol(counts)), function(j) {
    score <- quasi_fits[[j]]$score
    if (any(!is.finite(score)) || sum(score^2) <= 0) {
      stop(
        "Degenerate beta-binomial proportional-variance score for ",
        colnames(counts)[j]
      )
    }
    studentized_score_statistics(score, assignments)
  }, numeric(ncol(assignments)))
  if (is.null(dim(statistics))) statistics <- matrix(statistics, ncol = 1L)
  colnames(statistics) <- colnames(counts)

  exact_p <- apply(abs(statistics), 2L, inclusive_rank_p)
  if (is.null(dim(exact_p))) exact_p <- matrix(exact_p, ncol = 1L)
  colnames(exact_p) <- colnames(counts)
  normal_p <- 2 * stats::pnorm(-abs(statistics))

  list(
    statistics = statistics,
    exact_p = exact_p,
    normal_p = normal_p,
    null_fits = null_fits,
    quasi_fits = quasi_fits
  )
}

# Evaluate the score and Hessian of a group-specific intercept for all label
# assignments simultaneously.  `membership` is samples x assignments.
group_score_hessian <- function(eta, y, m, phi, membership) {
  n_assignments <- ncol(membership)
  if (length(phi) == 1L) phi <- rep(phi, n_assignments)
  mu <- stats::plogis(pmin(pmax(eta, -25), 25))
  phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
  tau <- (1 - phi) / phi
  shape1 <- pmax(mu * tau, 1e-12)
  shape2 <- pmax((1 - mu) * tau, 1e-12)
  common <- tau * mu * (1 - mu)
  score <- numeric(n_assignments)
  hessian <- numeric(n_assignments)

  for (i in seq_along(y)) {
    weight <- membership[i, ]
    first <-
      digamma(shape1 + y[i]) - digamma(shape1) -
      digamma(shape2 + m[i] - y[i]) + digamma(shape2)
    second <-
      trigamma(shape1 + y[i]) - trigamma(shape1) +
      trigamma(shape2 + m[i] - y[i]) - trigamma(shape2)
    score <- score + weight * common * first
    hessian <- hessian + weight * (
      common * (1 - 2 * mu) * first + common^2 * second
    )
  }
  list(score = score, hessian = hessian)
}

group_loglik_many <- function(eta, y, m, phi, membership) {
  n_assignments <- ncol(membership)
  if (length(phi) == 1L) phi <- rep(phi, n_assignments)
  mu <- stats::plogis(pmin(pmax(eta, -25), 25))
  phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
  tau <- (1 - phi) / phi
  shape1 <- pmax(mu * tau, 1e-12)
  shape2 <- pmax((1 - mu) * tau, 1e-12)
  loglik <- numeric(n_assignments)
  for (i in seq_along(y)) {
    contribution <-
      lchoose(m[i], y[i]) +
      lbeta(y[i] + shape1, m[i] - y[i] + shape2) -
      lbeta(shape1, shape2)
    loglik <- loglik + membership[i, ] * contribution
  }
  loglik
}

fit_group_eta_many <- function(y, m, phi, membership, tol = 1e-8) {
  group_y <- as.numeric(crossprod(y, membership))
  group_m <- as.numeric(crossprod(m, membership))
  eta <- stats::qlogis((group_y + 0.5) / (group_m + 1))
  boundary_low <- group_y == 0
  boundary_high <- group_y == group_m
  eta[boundary_low] <- -25
  eta[boundary_high] <- 25
  active <- !(boundary_low | boundary_high)

  for (iteration in seq_len(50L)) {
    if (!any(active)) break
    derivative <- group_score_hessian(eta, y, m, phi, membership)
    good <- active & is.finite(derivative$score) &
      is.finite(derivative$hessian) & derivative$hessian < -1e-12
    if (!any(good)) break
    step <- derivative$score[good] / derivative$hessian[good]
    step <- pmin(pmax(step, -5), 5)
    eta[good] <- pmin(pmax(eta[good] - step, -25), 25)
    converged <- good
    converged[good] <- abs(step) <= tol
    active[converged] <- FALSE
  }

  derivative <- group_score_hessian(eta, y, m, phi, membership)
  unresolved <- which(
    !(boundary_low | boundary_high) &
      (!is.finite(eta) | !is.finite(derivative$score) |
         abs(derivative$score) > 1e-5)
  )
  fallback <- 0L
  if (length(unresolved) > 0L) {
    for (b in unresolved) {
      keep <- membership[, b] == 1L
      fit <- fit_eta(y[keep], m[keep], if (length(phi) == 1L) phi else phi[b])
      eta[b] <- fit$eta
      fallback <- fallback + 1L
    }
  }

  list(
    eta = eta,
    loglik = group_loglik_many(eta, y, m, phi, membership),
    fallback = fallback,
    boundary = sum(boundary_low | boundary_high)
  )
}

fixed_phi_lrt <- function(y, m, phi, assignments, tol = 1e-8) {
  null_fit <- fit_eta(y, m, phi, tol = tol / 10)
  old_fit <- fit_group_eta_many(y, m, phi, assignments, tol = tol)
  young_fit <- fit_group_eta_many(y, m, phi, 1L - assignments, tol = tol)
  likelihood_ratio <- 2 * (
    old_fit$loglik + young_fit$loglik - null_fit$loglik
  )
  small_negative <- likelihood_ratio < 0 & likelihood_ratio > -1e-6
  likelihood_ratio[small_negative] <- 0
  likelihood_ratio[likelihood_ratio < -1e-6] <- NA_real_
  list(
    statistic = likelihood_ratio,
    fallback = old_fit$fallback + young_fit$fallback,
    boundary = old_fit$boundary + young_fit$boundary
  )
}

dcats_style_p_matrices <- function(
    counts,
    assignments,
    cell_null_fits,
    pooled_phi,
    cores = 1L,
    tol = 1e-8) {
  totals <- rowSums(counts)
  map_cells <- function(indices, fun) {
    if (.Platform$OS.type != "windows" && cores > 1L) {
      parallel::mclapply(
        indices, fun, mc.cores = cores, mc.preschedule = TRUE
      )
    } else {
      lapply(indices, fun)
    }
  }
  global_results <- map_cells(seq_len(ncol(counts)), function(j) {
    fixed_phi_lrt(counts[, j], totals, pooled_phi, assignments, tol = tol)
  })
  cell_results <- map_cells(seq_len(ncol(counts)), function(j) {
    fixed_phi_lrt(
      counts[, j], totals, cell_null_fits[[j]]$phi, assignments, tol = tol
    )
  })
  names(global_results) <- names(cell_results) <- colnames(counts)

  extract_statistics <- function(results) {
    out <- vapply(results, `[[`, numeric(ncol(assignments)), "statistic")
    if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
    colnames(out) <- colnames(counts)
    out
  }
  global_stat <- extract_statistics(global_results)
  cell_stat <- extract_statistics(cell_results)

  list(
    global_statistic = global_stat,
    cell_statistic = cell_stat,
    global_chisq_p = stats::pchisq(global_stat, df = 1, lower.tail = FALSE),
    cell_chisq_p = stats::pchisq(cell_stat, df = 1, lower.tail = FALSE),
    global_results = global_results,
    cell_results = cell_results
  )
}

calibration_rows <- function(
    analysis_set,
    counts,
    p_matrices,
    cell_null_fits,
    pooled_phi,
    thresholds) {
  totals <- rowSums(counts)
  diagnostics <- data.frame(
    cell_type = colnames(counts),
    total_cells = colSums(counts),
    nonzero_mice = colSums(counts > 0),
    mean_mouse_proportion = colMeans(counts / totals),
    cell_null_phi = vapply(cell_null_fits, `[[`, numeric(1L), "phi"),
    cell_phi_at_boundary = vapply(
      cell_null_fits, `[[`, logical(1L), "phi_at_boundary"
    ),
    pooled_null_phi = pooled_phi,
    stringsAsFactors = FALSE
  )

  rows <- list()
  for (method in names(p_matrices)) {
    p_values <- p_matrices[[method]]
    if (!identical(colnames(p_values), colnames(counts))) {
      stop("P-value columns do not match count columns for ", method)
    }
    for (alpha in thresholds) {
      finite <- is.finite(p_values)
      rejected <- finite & p_values <= alpha
      row <- diagnostics
      row$analysis_set <- analysis_set
      row$method <- method
      row$alpha <- alpha
      row$total_assignments <- nrow(p_values)
      row$evaluated_assignments <- colSums(finite)
      row$failure_rate <- colMeans(!finite)
      row$rejection_rate <- colSums(rejected) / pmax(colSums(finite), 1L)
      row$calibration_error <- row$rejection_rate - alpha
      row$mean_p_value <- colMeans(p_values, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- row
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

observed_discovery_rows <- function(
    analysis_set,
    result,
    observed_assignment_index,
    observed_old) {
  counts <- result$counts
  totals <- rowSums(counts)
  proportions <- sweep(counts, 1L, totals, FUN = "/")
  cell_types <- colnames(counts)
  p_matrices <- result$p_matrices
  get_observed_p <- function(method) {
    as.numeric(p_matrices[[method]][observed_assignment_index, ])
  }
  pcit_p <- get_observed_p("pcit_permuted_proportional_variance_score")
  normal_p <- get_observed_p("proportional_variance_score_normal_reference")
  global_p <- get_observed_p("dcats_lrt_global_phi_chisq")
  cell_p <- get_observed_p("dcats_lrt_cell_phi_chisq")
  score_ess <- vapply(result$score$quasi_fits, function(fit) {
    information <- (7 / 15) * (8 / 15) * fit$score^2
    sum(information)^2 / sum(information^2)
  }, numeric(1L))

  data.frame(
    analysis_set = analysis_set,
    cell_type = cell_types,
    total_cells = colSums(counts),
    nonzero_mice = colSums(counts > 0),
    young_mean_proportion = colMeans(
      proportions[observed_old == 0L, , drop = FALSE]
    ),
    old_mean_proportion = colMeans(
      proportions[observed_old == 1L, , drop = FALSE]
    ),
    old_minus_young_mean_proportion =
      colMeans(proportions[observed_old == 1L, , drop = FALSE]) -
      colMeans(proportions[observed_old == 0L, , drop = FALSE]),
    score_information_ess = score_ess,
    pcit_statistic = as.numeric(
      result$score$statistics[observed_assignment_index, ]
    ),
    pcit_p = pcit_p,
    pcit_q = stats::p.adjust(pcit_p, method = "BH"),
    normal_score_p = normal_p,
    normal_score_q = stats::p.adjust(normal_p, method = "BH"),
    dcats_global_null_phi_p = global_p,
    dcats_global_null_phi_q = stats::p.adjust(global_p, method = "BH"),
    dcats_cell_null_phi_p = cell_p,
    dcats_cell_null_phi_q = stats::p.adjust(cell_p, method = "BH"),
    stringsAsFactors = FALSE
  )
}

bh_calibration_rows <- function(analysis_set, p_matrices, alpha = 0.05) {
  rows <- lapply(names(p_matrices), function(method) {
    p_values <- p_matrices[[method]]
    rejection_count <- apply(p_values, 1L, function(p) {
      sum(stats::p.adjust(p, method = "BH") <= alpha, na.rm = TRUE)
    })
    data.frame(
      analysis_set = analysis_set,
      method = method,
      alpha = alpha,
      total_assignments = nrow(p_values),
      probability_any_bh_rejection = mean(rejection_count > 0L),
      mean_bh_rejections = mean(rejection_count),
      mean_bh_rejections_given_any = if (any(rejection_count > 0L)) {
        mean(rejection_count[rejection_count > 0L])
      } else {
        0
      },
      maximum_bh_rejections = max(rejection_count),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_lrt_against_aod <- function(
    counts,
    assignments,
    pooled_phi,
    custom_statistics) {
  if (!requireNamespace("aod", quietly = TRUE)) {
    return(list(performed = FALSE, reason = "aod is not installed"))
  }
  totals <- rowSums(counts)
  cell_index <- unique(round(seq(1, ncol(counts), length.out = 5)))
  assignment_index <- unique(round(seq(1, ncol(assignments), length.out = 7)))
  differences <- numeric()

  for (j in cell_index) {
    for (b in assignment_index) {
      dat <- data.frame(
        y = counts[, j],
        m = totals,
        condition = assignments[, b]
      )
      null <- suppressWarnings(aod::betabin(
        cbind(y, m - y) ~ 1, ~1, data = dat, warnings = FALSE,
        fixpar = list(2L, pooled_phi)
      ))
      full <- suppressWarnings(aod::betabin(
        cbind(y, m - y) ~ condition, ~1, data = dat, warnings = FALSE,
        fixpar = list(3L, pooled_phi)
      ))
      aod_lr <- null@dev - full@dev
      differences <- c(
        differences,
        custom_statistics[b, j] - aod_lr
      )
    }
  }
  list(
    performed = TRUE,
    comparisons = length(differences),
    max_abs_lrt_difference = max(abs(differences)),
    mean_abs_lrt_difference = mean(abs(differences))
  )
}

pretty_cell_type <- function(x) {
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("Cd4", "CD4", x, fixed = TRUE)
  x <- gsub("Cd11b", "CD11b", x, fixed = TRUE)
  x <- gsub("Cd103", "CD103", x, fixed = TRUE)
  x <- gsub("Type 2", "Type-2", x, fixed = TRUE)
  x <- gsub("Type1", "Type-1", x, fixed = TRUE)
  x
}

plot_primary_calibration <- function(calibration, path_pdf, path_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is unavailable; skipping calibration plot")
    return(invisible(FALSE))
  }
  methods <- c(
    "pcit_permuted_proportional_variance_score",
    "proportional_variance_score_normal_reference",
    "dcats_lrt_global_phi_chisq",
    "dcats_lrt_cell_phi_chisq"
  )
  labels <- c(
    pcit_permuted_proportional_variance_score = "PCIT: permuted proportional-variance score",
    proportional_variance_score_normal_reference = "Same score: normal reference",
    dcats_lrt_global_phi_chisq = "DCATS LRT: pooled phi",
    dcats_lrt_cell_phi_chisq = "DCATS LRT: cell-specific phi"
  )
  colors <- c(
    pcit_permuted_proportional_variance_score = "#0072B2",
    proportional_variance_score_normal_reference = "#CC79A7",
    dcats_lrt_global_phi_chisq = "#D55E00",
    dcats_lrt_cell_phi_chisq = "#009E73"
  )
  dat <- calibration[
    calibration$analysis_set == "paper_30" &
      calibration$alpha == 0.05 &
      calibration$method %in% methods,
    , drop = FALSE
  ]
  pooled <- dat[
    dat$method == "dcats_lrt_global_phi_chisq",
    c("cell_type", "rejection_rate")
  ]
  order_types <- pooled$cell_type[order(pooled$rejection_rate)]
  dat$cell_label <- factor(
    pretty_cell_type(dat$cell_type),
    levels = pretty_cell_type(order_types)
  )
  dat$method <- factor(dat$method, levels = methods, labels = labels[methods])

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = rejection_rate, y = cell_label, color = method)
  ) +
    ggplot2::geom_vline(
      xintercept = 0.05, linetype = "dashed", color = "grey35", linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 2.1, alpha = 0.9, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_color_manual(values = unname(colors[methods])) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(0.105, max(dat$rejection_rate) * 1.04)),
      breaks = scales::breaks_pretty(n = 6),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      x = "Conditional rejection fraction over all 6,435 label assignments",
      y = NULL,
      color = NULL,
      title = "Cell-type-specific calibration in the Angelidis lung-aging data",
      subtitle = "Paper-defined 30 types; nominal two-sided alpha = 0.05"
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey91", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(path_pdf, p, width = 9.2, height = 9.4, units = "in")
  ggplot2::ggsave(path_png, p, width = 9.2, height = 9.4, units = "in", dpi = 180)
  invisible(TRUE)
}

plot_observed_discoveries <- function(observed, path_pdf, path_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    warning("ggplot2/scales unavailable; skipping observed-discovery plot")
    return(invisible(FALSE))
  }
  primary <- observed[observed$analysis_set == "paper_30", , drop = FALSE]
  method_columns <- c(
    pcit_permuted_proportional_variance_score = "pcit_q",
    proportional_variance_score_normal_reference = "normal_score_q",
    dcats_lrt_global_phi_chisq = "dcats_global_null_phi_q",
    dcats_lrt_cell_phi_chisq = "dcats_cell_null_phi_q"
  )
  method_labels <- c(
    pcit_permuted_proportional_variance_score = "PCIT permutation score",
    proportional_variance_score_normal_reference = "Same score, normal reference",
    dcats_lrt_global_phi_chisq = "DCATS LRT, pooled null phi",
    dcats_lrt_cell_phi_chisq = "DCATS LRT, cell-specific null phi"
  )
  dat <- do.call(rbind, lapply(names(method_columns), function(method) {
    data.frame(
      cell_type = primary$cell_type,
      effect = primary$old_minus_young_mean_proportion,
      method = method,
      q_value = primary[[method_columns[[method]]]],
      stringsAsFactors = FALSE
    )
  }))
  order_types <- primary$cell_type[
    order(primary$old_minus_young_mean_proportion)
  ]
  dat$cell_label <- factor(
    pretty_cell_type(dat$cell_type),
    levels = pretty_cell_type(order_types)
  )
  dat$method <- factor(
    dat$method,
    levels = names(method_columns),
    labels = unname(method_labels[names(method_columns)])
  )
  dat$signed_evidence <- sign(dat$effect) * -log10(pmax(dat$q_value, 1e-15))
  colors <- c("#0072B2", "#CC79A7", "#D55E00", "#009E73")
  cutoff <- -log10(0.05)
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(signed_evidence, cell_label, color = method)
  ) +
    ggplot2::geom_vline(
      xintercept = c(-cutoff, cutoff), linetype = "dashed",
      color = "grey35", linewidth = 0.45
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    ggplot2::geom_point(
      size = 2.1, alpha = 0.9,
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(
      x = "Signed -log10(BH q); negative means lower in old mice",
      y = NULL,
      color = NULL,
      title = "Observed age-associated cell-composition discoveries in Angelidis",
      subtitle = "Corrected 7 old / 8 young labels; paper-defined 30 cell types"
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey91", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(path_pdf, p, width = 9.2, height = 9.4, units = "in")
  ggplot2::ggsave(path_png, p, width = 9.2, height = 9.4, units = "in", dpi = 180)
  invisible(TRUE)
}

run_analysis_set <- function(name, counts, assignments, cores, tol) {
  cat("\n", name, ": ", nrow(counts), " mice, ", ncol(counts),
      " cell types, ", sum(counts), " retained cells\n", sep = "")
  cat("  Fitting cell-specific null beta-binomial models...\n")
  score_result <- proportional_variance_score_p_matrices(
    counts, assignments, tol = tol
  )
  cat("  Profiling a pooled null dispersion...\n")
  pooled_fit <- fit_pooled_null_phi(counts, rowSums(counts), tol = tol)
  if (pooled_fit$phi_at_boundary) stop("Pooled phi estimate hit its search boundary")
  cat("  pooled phi = ", format(pooled_fit$phi, digits = 8), "\n", sep = "")
  cat("  Computing complete-orbit DCATS-style LRTs...\n")
  dcats_result <- dcats_style_p_matrices(
    counts,
    assignments,
    score_result$null_fits,
    pooled_fit$phi,
    cores = cores,
    tol = tol
  )

  p_matrices <- list(
    pcit_permuted_proportional_variance_score = score_result$exact_p,
    proportional_variance_score_normal_reference = score_result$normal_p,
    dcats_lrt_global_phi_chisq = dcats_result$global_chisq_p,
    dcats_lrt_cell_phi_chisq = dcats_result$cell_chisq_p
  )
  thresholds <- c(0.01, 0.025, 0.05, 0.10)
  calibration <- calibration_rows(
    name, counts, p_matrices, score_result$null_fits,
    pooled_fit$phi, thresholds
  )
  list(
    name = name,
    counts = counts,
    score = score_result,
    pooled_fit = pooled_fit,
    dcats = dcats_result,
    p_matrices = p_matrices,
    calibration = calibration
  )
}

run_main <- function(cfg) {
  if (dir.exists(cfg$aod_lib)) .libPaths(c(cfg$aod_lib, .libPaths()))
  angelidis <- read_angelidis(cfg$metadata)
  assignments <- make_assignments(15L, 7L)
  if (ncol(assignments) != 6435L || any(colSums(assignments) != 7L)) {
    stop("Assignment-orbit construction failed")
  }
  observed_assignment_index <- which(colSums(
    abs(assignments - angelidis$observed_old)
  ) == 0L)
  if (length(observed_assignment_index) != 1L) {
    stop("Corrected observed age assignment is not unique in the orbit")
  }

  analysis_sets <- list(paper_30 = angelidis$counts_30)
  if (cfg$include_33_sensitivity) {
    analysis_sets$metadata_33 <- angelidis$counts_33
  }

  start_time <- proc.time()[[3L]]
  results <- lapply(names(analysis_sets), function(name) {
    run_analysis_set(
      name, analysis_sets[[name]], assignments, cfg$cores, cfg$tolerance
    )
  })
  names(results) <- names(analysis_sets)
  elapsed_seconds <- proc.time()[[3L]] - start_time

  calibration <- do.call(rbind, lapply(results, `[[`, "calibration"))
  rownames(calibration) <- NULL
  observed <- do.call(rbind, lapply(names(results), function(name) {
    observed_discovery_rows(
      name,
      results[[name]],
      observed_assignment_index,
      angelidis$observed_old
    )
  }))
  rownames(observed) <- NULL
  bh_calibration <- do.call(rbind, lapply(names(results), function(name) {
    bh_calibration_rows(name, results[[name]]$p_matrices, alpha = 0.05)
  }))
  rownames(bh_calibration) <- NULL

  validation <- list(performed = FALSE, reason = "disabled")
  if (cfg$validate_aod) {
    primary <- results$paper_30
    validation <- validate_lrt_against_aod(
      primary$counts,
      assignments,
      primary$pooled_fit$phi,
      primary$dcats$global_statistic
    )
    if (isTRUE(validation$performed) &&
        validation$max_abs_lrt_difference > 1e-4) {
      stop(
        "Custom fixed-phi LRT did not reproduce aod/DCATS closely enough: ",
        validation$max_abs_lrt_difference
      )
    }
  }

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$figure_dir, recursive = TRUE, showWarnings = FALSE)

  calibration_path <- file.path(
    cfg$out_dir, "angelidis-celltype-permutation-calibration.csv"
  )
  summary_path <- file.path(
    cfg$out_dir, "angelidis-celltype-permutation-calibration-alpha05.csv"
  )
  observed_path <- file.path(
    cfg$out_dir, "angelidis-celltype-observed-discoveries.csv"
  )
  bh_calibration_path <- file.path(
    cfg$out_dir, "angelidis-celltype-bh-calibration.csv"
  )
  rds_path <- file.path(
    cfg$out_dir, "angelidis-celltype-permutation-calibration.rds"
  )
  figure_pdf <- file.path(cfg$figure_dir, "angelidis-celltype-calibration.pdf")
  figure_png <- file.path(cfg$figure_dir, "angelidis-celltype-calibration.png")
  discovery_figure_pdf <- file.path(
    cfg$figure_dir, "angelidis-observed-discoveries.pdf"
  )
  discovery_figure_png <- file.path(
    cfg$figure_dir, "angelidis-observed-discoveries.png"
  )

  utils::write.csv(calibration, calibration_path, row.names = FALSE)
  utils::write.csv(
    calibration[calibration$alpha == 0.05, , drop = FALSE],
    summary_path,
    row.names = FALSE
  )
  utils::write.csv(observed, observed_path, row.names = FALSE)
  utils::write.csv(bh_calibration, bh_calibration_path, row.names = FALSE)

  artifact <- list(
    schema_version = "angelidis_celltype_permutation_calibration_v1",
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    metadata_path = normalizePath(cfg$metadata),
    metadata_md5 = unname(tools::md5sum(cfg$metadata)),
    metadata_rows = angelidis$metadata_rows,
    mouse_age = angelidis$mouse_age,
    observed_old = angelidis$observed_old,
    observed_assignment_index = observed_assignment_index,
    excluded_from_paper_30 = angelidis$excluded_from_paper_30,
    assignment_description = "all 7-of-15 old-label allocations",
    x_given_z_probability = 7 / 15,
    nuisance_covariates = character(),
    assignments = assignments,
    results = results,
    calibration = calibration,
    observed = observed,
    bh_calibration = bh_calibration,
    aod_validation = validation,
    elapsed_seconds = elapsed_seconds,
    cores = cfg$cores,
    session_info = utils::sessionInfo()
  )
  saveRDS(artifact, rds_path, compress = "xz")
  plot_primary_calibration(calibration, figure_pdf, figure_png)
  plot_observed_discoveries(
    observed, discovery_figure_pdf, discovery_figure_png
  )

  cat("\nPrimary paper-30 calibration at alpha = 0.05:\n")
  primary_alpha05 <- calibration[
    calibration$analysis_set == "paper_30" & calibration$alpha == 0.05,
    , drop = FALSE
  ]
  method_summary <- do.call(rbind, lapply(
    split(primary_alpha05, primary_alpha05$method),
    function(x) data.frame(
      method = x$method[[1L]],
      mean_rate = mean(x$rejection_rate),
      min_rate = min(x$rejection_rate),
      max_rate = max(x$rejection_rate),
      max_abs_error = max(abs(x$calibration_error)),
      failure_rate = mean(x$failure_rate),
      stringsAsFactors = FALSE
    )
  ))
  print(method_summary, row.names = FALSE, digits = 5)
  cat("\nPrimary paper-30 observed BH discoveries at q <= 0.05:\n")
  primary_observed <- observed[observed$analysis_set == "paper_30", , drop = FALSE]
  discovery_columns <- c(
    pcit_permuted_proportional_variance_score = "pcit_q",
    proportional_variance_score_normal_reference = "normal_score_q",
    dcats_lrt_global_phi_chisq = "dcats_global_null_phi_q",
    dcats_lrt_cell_phi_chisq = "dcats_cell_null_phi_q"
  )
  for (method in names(discovery_columns)) {
    q_value <- primary_observed[[discovery_columns[[method]]]]
    cat("  ", method, ": ", sum(q_value <= 0.05), "\n", sep = "")
  }
  cat("\nPrimary paper-30 family-level BH null calibration:\n")
  print(
    bh_calibration[bh_calibration$analysis_set == "paper_30", ],
    row.names = FALSE,
    digits = 5
  )
  cat("\nAOD validation:\n")
  print(validation)
  cat("\nElapsed: ", format(elapsed_seconds, digits = 5), " seconds\n", sep = "")
  cat("Wrote ", calibration_path, "\n", sep = "")
  cat("Wrote ", summary_path, "\n", sep = "")
  cat("Wrote ", observed_path, "\n", sep = "")
  cat("Wrote ", bh_calibration_path, "\n", sep = "")
  cat("Wrote ", rds_path, "\n", sep = "")
  if (file.exists(figure_pdf)) cat("Wrote ", figure_pdf, "\n", sep = "")
  if (file.exists(figure_png)) cat("Wrote ", figure_png, "\n", sep = "")
  if (file.exists(discovery_figure_pdf)) {
    cat("Wrote ", discovery_figure_pdf, "\n", sep = "")
  }
  if (file.exists(discovery_figure_png)) {
    cat("Wrote ", discovery_figure_png, "\n", sep = "")
  }
  invisible(artifact)
}

if (sys.nframe() == 0L) run_main(cfg)
