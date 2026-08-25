#!/usr/bin/env Rscript

# Covariate-adjusted conditional calibration on the DCATS-packaged Ren et al.
# (2021) PBMC composition data.  The primary contrast is control versus
# severe/critical progression (20 versus 48 libraries; 10 major cell types).
#
# PCIT specification requested for this analysis:
#   * Z = age + sex + fresh/frozen PBMC enters only the null Y | Z fit;
#   * X is randomized marginally, with 48 of 68 labels in every draw;
#   * no X | Z model, propensity fit, restricted permutation, or exposure
#     residualization is used;
#   * the score is the beta-binomial proportional-variance quasi-score;
#   * empirical contribution variance studentizes every randomized statistic.
#
# DCATS comparators retain their native one-df chi-square calibration.  We do
# not rank/permutation-calibrate a DCATS statistic.  Outer randomized labels
# are used only to measure the conditional null rejection rate of its native
# p-values on this fixed count/covariate table.

parse_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1L]])
}

as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

script_path <- sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
)
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- 1L
output_dir <- parse_arg(
  "out-dir",
  file.path(project_root, "data", "ren2021")
)

cfg <- list(
  ren_data = parse_arg("ren-data", "/tmp/DCATS/data/Ren2021.rda"),
  outer_draws = as.integer(parse_arg("outer-draws", "1000")),
  reference_draws = as.integer(parse_arg("reference-draws", "19999")),
  seed = as.integer(parse_arg("seed", "20260822")),
  cores = as.integer(parse_arg(
    "cores", as.character(min(6L, max(1L, detected_cores - 1L)))
  )),
  validate_aod = as_flag(parse_arg("validate-aod", "true")),
  aod_lib = parse_arg("aod-lib", "/private/tmp/aod-lib"),
  out_dir = output_dir,
  figure_dir = parse_arg(
    "figure-dir",
    file.path(output_dir, "figures")
  ),
  tolerance = as.numeric(parse_arg("tolerance", "1e-8"))
)

if (!is.finite(cfg$outer_draws) || cfg$outer_draws < 20L) {
  stop("--outer-draws must be at least 20", call. = FALSE)
}
if (!is.finite(cfg$reference_draws) || cfg$reference_draws < 99L) {
  stop("--reference-draws must be at least 99", call. = FALSE)
}
if (!is.finite(cfg$cores) || cfg$cores < 1L) stop("--cores must be positive")
if (!is.finite(cfg$tolerance) || cfg$tolerance <= 0) {
  stop("--tolerance must be positive")
}

load_ren2021 <- function(path) {
  if (file.exists(path)) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("Ren2021", envir = env, inherits = FALSE)) {
      stop("--ren-data does not contain an object named Ren2021", call. = FALSE)
    }
    return(get("Ren2021", envir = env, inherits = FALSE))
  }
  if (requireNamespace("DCATS", quietly = TRUE)) {
    env <- new.env(parent = emptyenv())
    suppressMessages(utils::data("Ren2021", package = "DCATS", envir = env))
    if (exists("Ren2021", envir = env, inherits = FALSE)) {
      return(get("Ren2021", envir = env, inherits = FALSE))
    }
  }
  stop(
    "Ren2021 data not found. Pass --ren-data=/absolute/path/to/Ren2021.rda ",
    "or install DCATS.",
    call. = FALSE
  )
}

prepare_ren <- function(ren) {
  required <- c("age", "sex", "sample_type", "state")
  if (!all(required %in% names(ren$designM))) {
    stop("Ren2021 design is missing an expected variable", call. = FALSE)
  }
  states <- c("control", "severe/critical_progression")
  keep <- ren$designM$state %in% states
  design <- ren$designM[keep, required, drop = FALSE]
  counts <- as.matrix(ren$countM[rownames(design), , drop = FALSE])
  storage.mode(counts) <- "integer"
  if (!identical(rownames(counts), rownames(design))) {
    stop("Count and design rows are misaligned", call. = FALSE)
  }
  if (anyNA(design) || any(rowSums(counts) <= 0L)) {
    stop("Selected Ren2021 rows contain missing design values or zero totals")
  }
  design$age_scaled <- as.numeric(scale(as.numeric(design$age)))
  design$sex <- droplevels(factor(design$sex))
  design$sample_type <- droplevels(factor(design$sample_type))
  design$x <- as.integer(design$state == "severe/critical_progression")
  z <- stats::model.matrix(~ age_scaled + sex + sample_type, data = design)
  if (qr(z)$rank != ncol(z)) stop("The nuisance design is rank deficient")
  if (!identical(as.integer(table(design$x)), c(20L, 48L))) {
    stop("Expected 20 control and 48 severe-progression libraries")
  }
  list(
    counts = counts,
    totals = rowSums(counts),
    design = design,
    z = z,
    x = design$x
  )
}

sample_fixed_allocations <- function(n, treated, draws) {
  out <- matrix(0L, nrow = n, ncol = draws)
  for (b in seq_len(draws)) out[sample.int(n, treated), b] <- 1L
  out
}

bb_loglik <- function(beta, y, m, phi, design) {
  eta <- as.numeric(design %*% beta)
  mu <- stats::plogis(pmin(pmax(eta, -30), 30))
  phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
  tau <- (1 - phi) / phi
  a <- pmax(mu * tau, 1e-12)
  b <- pmax((1 - mu) * tau, 1e-12)
  sum(lchoose(m, y) + lbeta(y + a, m - y + b) - lbeta(a, b))
}

glm_start <- function(y, m, design) {
  fit <- suppressWarnings(stats::glm.fit(
    x = design,
    y = y / m,
    weights = m,
    family = stats::binomial()
  ))
  beta <- fit$coefficients
  beta[!is.finite(beta)] <- 0
  if (length(beta) != ncol(design)) beta <- numeric(ncol(design))
  beta
}

fit_beta_fixed_phi <- function(
    y, m, phi, design, start = NULL, tol = 1e-8, floor_loglik = -Inf) {
  if (is.null(start)) start <- glm_start(y, m, design)
  objective <- function(beta) -bb_loglik(beta, y, m, phi, design)
  starts <- list(start)
  glm_beta <- glm_start(y, m, design)
  if (max(abs(glm_beta - start)) > 1e-8) starts[[2L]] <- glm_beta
  fits <- lapply(starts, function(initial) {
    try(stats::optim(
      initial,
      objective,
      method = "BFGS",
      control = list(reltol = tol / 10, maxit = 1000L)
    ), silent = TRUE)
  })
  usable <- vapply(fits, function(fit) {
    !inherits(fit, "try-error") && is.finite(fit$value)
  }, logical(1L))
  if (!any(usable)) {
    return(list(beta = start, loglik = floor_loglik, convergence = 99L))
  }
  fits <- fits[usable]
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "value"))]]
  loglik <- -best$value
  if (is.finite(floor_loglik) && loglik < floor_loglik) {
    return(list(beta = start, loglik = floor_loglik, convergence = 98L))
  }
  list(beta = best$par, loglik = loglik, convergence = best$convergence)
}

fit_cell_null <- function(y, m, z, tol = 1e-8) {
  beta_start <- glm_start(y, m, z)
  profile <- function(logit_phi) {
    phi <- stats::plogis(logit_phi)
    fit_beta_fixed_phi(y, m, phi, z, start = beta_start, tol = tol)$loglik
  }
  fit_phi <- stats::optimize(
    profile,
    c(stats::qlogis(1e-6), stats::qlogis(0.6)),
    maximum = TRUE,
    tol = tol
  )
  phi <- stats::plogis(fit_phi$maximum)
  beta_fit <- fit_beta_fixed_phi(y, m, phi, z, start = beta_start, tol = tol)
  list(
    beta = beta_fit$beta,
    phi = phi,
    loglik = beta_fit$loglik,
    convergence = beta_fit$convergence,
    phi_at_boundary = phi <= 1.01e-6 || phi >= 0.599
  )
}

fit_pooled_null <- function(counts, m, z, tol = 1e-8) {
  starts <- lapply(seq_len(ncol(counts)), function(j) glm_start(counts[, j], m, z))
  profile <- function(logit_phi) {
    phi <- stats::plogis(logit_phi)
    sum(vapply(seq_len(ncol(counts)), function(j) {
      fit_beta_fixed_phi(
        counts[, j], m, phi, z, start = starts[[j]], tol = tol
      )$loglik
    }, numeric(1L)))
  }
  fit_phi <- stats::optimize(
    profile,
    c(stats::qlogis(1e-6), stats::qlogis(0.5)),
    maximum = TRUE,
    tol = tol
  )
  phi <- stats::plogis(fit_phi$maximum)
  beta_fits <- lapply(seq_len(ncol(counts)), function(j) {
    fit_beta_fixed_phi(counts[, j], m, phi, z, start = starts[[j]], tol = tol)
  })
  list(phi = phi, loglik = fit_phi$objective, beta_fits = beta_fits)
}

# Solve Z'[(Y-M*mu)/{1+(M-1)phi}] = 0.  This is the requested
# proportional-variance null Y | Z fit; X never enters this routine.
fit_proportional_variance_score <- function(
    y, m, phi, z, beta_start, tol = 1e-10) {
  beta <- beta_start
  d <- 1 + (m - 1) * phi
  converged <- FALSE
  for (iteration in seq_len(100L)) {
    mu <- stats::plogis(pmin(pmax(as.numeric(z %*% beta), -30), 30))
    score <- (y - m * mu) / d
    information <- crossprod(z, z * as.numeric(m * mu * (1 - mu) / d))
    update <- try(solve(information, crossprod(z, score)), silent = TRUE)
    if (inherits(update, "try-error") || any(!is.finite(update))) break
    update <- as.numeric(update)
    if (max(abs(update)) > 2) update <- update * 2 / max(abs(update))
    beta <- beta + update
    if (max(abs(update)) <= tol) {
      converged <- TRUE
      break
    }
  }
  mu <- stats::plogis(pmin(pmax(as.numeric(z %*% beta), -30), 30))
  score <- (y - m * mu) / d
  equation_error <- max(abs(crossprod(z, score)))
  if (!is.finite(equation_error) || equation_error > 1e-5) {
    stop("Proportional-variance null score equation did not converge")
  }
  list(
    beta = beta,
    mu = mu,
    score = score,
    variance_multiplier = d,
    equation_error = equation_error,
    converged = converged,
    iterations = iteration
  )
}

studentized_statistics <- function(score, allocations, treated_fraction) {
  centered <- allocations - treated_fraction
  contributions <- sweep(centered, 1L, score, FUN = "*")
  totals <- colSums(contributions)
  denominator <- sqrt(pmax(
    colSums(contributions^2) - totals^2 / nrow(contributions),
    0
  ))
  statistic <- totals / denominator
  statistic[!is.finite(statistic) | denominator <= 0] <- NA_real_
  statistic
}

reference_p_values <- function(observed, reference) {
  reference <- abs(reference[is.finite(reference)])
  vapply(abs(observed), function(value) {
    if (!is.finite(value)) return(NA_real_)
    (1 + sum(reference >= value - 1e-12)) / (length(reference) + 1)
  }, numeric(1L))
}

fit_lrt_many <- function(
    y, m, z, allocations, phi, null_fit, tol = 1e-8) {
  statistic <- p_value <- rep(NA_real_, ncol(allocations))
  convergence <- integer(ncol(allocations))
  null_start <- c(null_fit$beta, 0)
  for (b in seq_len(ncol(allocations))) {
    full_design <- cbind(z, condition = allocations[, b])
    full <- fit_beta_fixed_phi(
      y,
      m,
      phi,
      full_design,
      start = null_start,
      tol = tol,
      floor_loglik = null_fit$loglik
    )
    lr <- 2 * (full$loglik - null_fit$loglik)
    if (is.finite(lr) && lr >= -1e-6) {
      statistic[b] <- max(lr, 0)
      p_value[b] <- stats::pchisq(statistic[b], 1, lower.tail = FALSE)
    }
    convergence[b] <- full$convergence
  }
  list(statistic = statistic, p_value = p_value, convergence = convergence)
}

map_cells <- function(indices, fun, cores) {
  if (.Platform$OS.type != "windows" && cores > 1L) {
    parallel::mclapply(indices, fun, mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(indices, fun)
  }
}

calibration_table <- function(p_matrices, diagnostics, thresholds) {
  rows <- list()
  for (method in names(p_matrices)) {
    p <- p_matrices[[method]]
    for (alpha in thresholds) {
      finite <- is.finite(p)
      row <- diagnostics
      row$method <- method
      row$alpha <- alpha
      row$evaluated_draws <- colSums(finite)
      row$failure_rate <- colMeans(!finite)
      row$rejection_rate <- colSums(finite & p <= alpha) / pmax(colSums(finite), 1L)
      row$calibration_error <- row$rejection_rate - alpha
      row$mean_p_value <- colMeans(p, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- row
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

validate_against_aod <- function(
    prepared, allocations, pooled_phi, pooled_null_fits, custom_statistics) {
  if (!requireNamespace("aod", quietly = TRUE)) {
    return(list(performed = FALSE, reason = "aod unavailable"))
  }
  cell_index <- c(1L, 5L, 10L)
  draw_index <- c(1L, max(1L, ncol(allocations) %/% 2L), ncol(allocations))
  differences <- numeric()
  for (j in cell_index) {
    for (b in draw_index) {
      dat <- data.frame(
        y = prepared$counts[, j],
        m = prepared$totals,
        age_scaled = prepared$design$age_scaled,
        sex = prepared$design$sex,
        sample_type = prepared$design$sample_type,
        condition = allocations[, b]
      )
      null <- suppressWarnings(aod::betabin(
        cbind(y, m - y) ~ age_scaled + sex + sample_type,
        ~1,
        data = dat,
        warnings = FALSE,
        fixpar = list(5L, pooled_phi)
      ))
      full <- suppressWarnings(aod::betabin(
        cbind(y, m - y) ~ age_scaled + sex + sample_type + condition,
        ~1,
        data = dat,
        warnings = FALSE,
        fixpar = list(6L, pooled_phi)
      ))
      differences <- c(
        differences,
        custom_statistics[b, j] - (null@dev - full@dev)
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

plot_calibration <- function(calibration, path_pdf, path_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) return(invisible(FALSE))
  methods <- c(
    "pcit_permuted_proportional_variance_score",
    "proportional_variance_score_normal_reference",
    "dcats_global_null_phi_lrt_chisq",
    "dcats_cell_null_phi_lrt_chisq"
  )
  labels <- c(
    pcit_permuted_proportional_variance_score = "PCIT permutation score",
    proportional_variance_score_normal_reference = "Same score, normal reference",
    dcats_global_null_phi_lrt_chisq = "DCATS LRT, pooled null phi",
    dcats_cell_null_phi_lrt_chisq = "DCATS LRT, cell-specific null phi"
  )
  colors <- c("#0072B2", "#CC79A7", "#D55E00", "#009E73")
  dat <- calibration[
    calibration$alpha == 0.05 & calibration$method %in% methods,
    , drop = FALSE
  ]
  ordering <- dat$cell_type[
    dat$method == "dcats_global_null_phi_lrt_chisq"
  ]
  ordering <- ordering[order(dat$rejection_rate[
    dat$method == "dcats_global_null_phi_lrt_chisq"
  ])]
  dat$cell_type <- factor(dat$cell_type, levels = ordering)
  dat$method <- factor(dat$method, levels = methods, labels = labels[methods])
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(rejection_rate, cell_type, color = method)
  ) +
    ggplot2::geom_vline(
      xintercept = 0.05, linetype = "dashed", color = "grey35", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      size = 2.4, position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 6),
      limits = c(0, max(0.10, max(dat$rejection_rate) * 1.05)),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      x = "Rejection fraction over marginal fixed-allocation null draws",
      y = NULL,
      color = NULL,
      title = "Covariate-adjusted calibration in the Ren COVID-19 PBMC data",
      subtitle = "Y|Z adjusts age, sex, and fresh/frozen status; nominal alpha = 0.05"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(path_pdf, p, width = 9.2, height = 5.8)
  ggplot2::ggsave(path_png, p, width = 9.2, height = 5.8, dpi = 180)
  invisible(TRUE)
}

run_main <- function(cfg) {
  if (dir.exists(cfg$aod_lib)) .libPaths(c(cfg$aod_lib, .libPaths()))
  set.seed(cfg$seed)
  prepared <- prepare_ren(load_ren2021(cfg$ren_data))
  n <- nrow(prepared$counts)
  treated <- sum(prepared$x)
  treated_fraction <- treated / n
  cat(
    "Ren2021 control vs severe/critical progression: ", n, " libraries, ",
    ncol(prepared$counts), " cell types; allocation ", n - treated, " vs ",
    treated, "\n", sep = ""
  )
  cat("Z: age + sex + fresh/frozen PBMC; no X|Z fit\n")

  outer_x <- sample_fixed_allocations(n, treated, cfg$outer_draws)
  reference_x <- sample_fixed_allocations(n, treated, cfg$reference_draws)

  cat("Fitting cell-specific Y|Z null models and proportional-variance scores...\n")
  cell_null <- lapply(seq_len(ncol(prepared$counts)), function(j) {
    fit_cell_null(
      prepared$counts[, j], prepared$totals, prepared$z, tol = cfg$tolerance
    )
  })
  names(cell_null) <- colnames(prepared$counts)
  if (any(vapply(cell_null, `[[`, logical(1L), "phi_at_boundary"))) {
    warning("At least one cell-specific phi estimate is on the search boundary")
  }
  quasi <- lapply(seq_len(ncol(prepared$counts)), function(j) {
    fit_proportional_variance_score(
      prepared$counts[, j],
      prepared$totals,
      cell_null[[j]]$phi,
      prepared$z,
      cell_null[[j]]$beta
    )
  })
  names(quasi) <- colnames(prepared$counts)

  score_outer <- vapply(quasi, function(fit) {
    studentized_statistics(fit$score, outer_x, treated_fraction)
  }, numeric(cfg$outer_draws))
  score_reference <- vapply(quasi, function(fit) {
    studentized_statistics(fit$score, reference_x, treated_fraction)
  }, numeric(cfg$reference_draws))
  colnames(score_outer) <- colnames(score_reference) <- colnames(prepared$counts)
  pcit_p <- vapply(seq_len(ncol(score_outer)), function(j) {
    reference_p_values(score_outer[, j], score_reference[, j])
  }, numeric(cfg$outer_draws))
  colnames(pcit_p) <- colnames(prepared$counts)
  normal_p <- 2 * stats::pnorm(-abs(score_outer))

  cat("Profiling the pooled label-blind Y|Z dispersion...\n")
  pooled <- fit_pooled_null(
    prepared$counts, prepared$totals, prepared$z, tol = cfg$tolerance
  )
  cat("Pooled null phi: ", format(pooled$phi, digits = 7), "\n", sep = "")

  cat(
    "Computing native chi-square DCATS-style LRTs for ", cfg$outer_draws,
    " marginal null allocations on ", cfg$cores, " core(s)...\n", sep = ""
  )
  start_time <- proc.time()[[3L]]
  lrt_results <- map_cells(seq_len(ncol(prepared$counts)), function(j) {
    global_null <- pooled$beta_fits[[j]]
    cell_fit <- cell_null[[j]]
    list(
      global = fit_lrt_many(
        prepared$counts[, j], prepared$totals, prepared$z, outer_x,
        pooled$phi, global_null, tol = cfg$tolerance
      ),
      cell = fit_lrt_many(
        prepared$counts[, j], prepared$totals, prepared$z, outer_x,
        cell_fit$phi, cell_fit, tol = cfg$tolerance
      )
    )
  }, cfg$cores)
  names(lrt_results) <- colnames(prepared$counts)
  lrt_elapsed <- proc.time()[[3L]] - start_time
  extract_lrt <- function(component, field) {
    out <- vapply(lrt_results, function(result) {
      result[[component]][[field]]
    }, numeric(cfg$outer_draws))
    colnames(out) <- colnames(prepared$counts)
    out
  }
  global_stat <- extract_lrt("global", "statistic")
  global_p <- extract_lrt("global", "p_value")
  cell_stat <- extract_lrt("cell", "statistic")
  cell_p <- extract_lrt("cell", "p_value")

  score_ess <- vapply(quasi, function(fit) {
    information <- treated_fraction * (1 - treated_fraction) * fit$score^2
    sum(information)^2 / sum(information^2)
  }, numeric(1L))
  diagnostics <- data.frame(
    cell_type = colnames(prepared$counts),
    total_cells = colSums(prepared$counts),
    nonzero_libraries = colSums(prepared$counts > 0),
    mean_library_proportion = colMeans(prepared$counts / prepared$totals),
    cell_null_phi = vapply(cell_null, `[[`, numeric(1L), "phi"),
    pooled_null_phi = pooled$phi,
    score_information_ess = score_ess,
    score_equation_max_abs = vapply(quasi, `[[`, numeric(1L), "equation_error"),
    nominal_libraries = n,
    control_libraries = n - treated,
    severe_progression_libraries = treated,
    outer_draws = cfg$outer_draws,
    reference_draws = cfg$reference_draws,
    stringsAsFactors = FALSE
  )
  p_matrices <- list(
    pcit_permuted_proportional_variance_score = pcit_p,
    proportional_variance_score_normal_reference = normal_p,
    dcats_global_null_phi_lrt_chisq = global_p,
    dcats_cell_null_phi_lrt_chisq = cell_p
  )
  calibration <- calibration_table(
    p_matrices, diagnostics, c(0.01, 0.025, 0.05, 0.10)
  )

  # Observed-label association results use the same marginal reference for
  # PCIT and native chi-square references for the two DCATS comparators.
  observed_x <- matrix(prepared$x, ncol = 1L)
  observed_score_stat <- vapply(quasi, function(fit) {
    studentized_statistics(fit$score, observed_x, treated_fraction)
  }, numeric(1L))
  observed_pcit <- vapply(seq_along(quasi), function(j) {
    reference_p_values(observed_score_stat[j], score_reference[, j])
  }, numeric(1L))
  observed_global <- vapply(seq_len(ncol(prepared$counts)), function(j) {
    fit_lrt_many(
      prepared$counts[, j], prepared$totals, prepared$z, observed_x,
      pooled$phi, pooled$beta_fits[[j]], tol = cfg$tolerance
    )$p_value
  }, numeric(1L))
  observed_cell <- vapply(seq_len(ncol(prepared$counts)), function(j) {
    fit_lrt_many(
      prepared$counts[, j], prepared$totals, prepared$z, observed_x,
      cell_null[[j]]$phi, cell_null[[j]], tol = cfg$tolerance
    )$p_value
  }, numeric(1L))
  observed_normal <- 2 * stats::pnorm(-abs(observed_score_stat))
  observed <- data.frame(
    cell_type = colnames(prepared$counts),
    control_mean_proportion = colMeans(
      prepared$counts[prepared$x == 0L, , drop = FALSE] /
        prepared$totals[prepared$x == 0L]
    ),
    severe_progression_mean_proportion = colMeans(
      prepared$counts[prepared$x == 1L, , drop = FALSE] /
        prepared$totals[prepared$x == 1L]
    ),
    pcit_statistic = observed_score_stat,
    pcit_p = observed_pcit,
    pcit_q = stats::p.adjust(observed_pcit, "BH"),
    normal_score_p = observed_normal,
    normal_score_q = stats::p.adjust(observed_normal, "BH"),
    dcats_global_null_phi_p = observed_global,
    dcats_global_null_phi_q = stats::p.adjust(observed_global, "BH"),
    dcats_cell_null_phi_p = observed_cell,
    dcats_cell_null_phi_q = stats::p.adjust(observed_cell, "BH"),
    stringsAsFactors = FALSE
  )

  validation <- list(performed = FALSE, reason = "disabled")
  if (cfg$validate_aod) {
    validation <- validate_against_aod(
      prepared, outer_x, pooled$phi, pooled$beta_fits, global_stat
    )
    if (isTRUE(validation$performed) &&
        validation$max_abs_lrt_difference > 1e-3) {
      stop(
        "Custom DCATS-style LRT differs from aod by more than tolerance: ",
        validation$max_abs_lrt_difference
      )
    }
  }

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$figure_dir, recursive = TRUE, showWarnings = FALSE)
  calibration_path <- file.path(cfg$out_dir, "ren2021-covariate-calibration.csv")
  alpha05_path <- file.path(cfg$out_dir, "ren2021-covariate-calibration-alpha05.csv")
  observed_path <- file.path(cfg$out_dir, "ren2021-covariate-observed.csv")
  rds_path <- file.path(cfg$out_dir, "ren2021-covariate-calibration.rds")
  figure_pdf <- file.path(cfg$figure_dir, "ren2021-covariate-calibration.pdf")
  figure_png <- file.path(cfg$figure_dir, "ren2021-covariate-calibration.png")
  utils::write.csv(calibration, calibration_path, row.names = FALSE)
  utils::write.csv(
    calibration[calibration$alpha == 0.05, , drop = FALSE],
    alpha05_path,
    row.names = FALSE
  )
  utils::write.csv(observed, observed_path, row.names = FALSE)
  artifact <- list(
    schema_version = "ren2021_covariate_calibration_v1",
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    data_path = if (file.exists(cfg$ren_data)) normalizePath(cfg$ren_data) else NA_character_,
    data_md5 = if (file.exists(cfg$ren_data)) unname(tools::md5sum(cfg$ren_data)) else NA_character_,
    contrast = "control_vs_severe_critical_progression",
    nuisance_covariates = c("age", "sex", "sample_type"),
    x_model = "none; unrestricted marginal fixed-allocation permutations",
    treated_fraction = treated_fraction,
    prepared = prepared,
    outer_allocations = outer_x,
    reference_allocations = reference_x,
    cell_null = cell_null,
    quasi = quasi,
    pooled_null = pooled,
    score_outer = score_outer,
    score_reference = score_reference,
    global_lrt_statistic = global_stat,
    cell_lrt_statistic = cell_stat,
    p_matrices = p_matrices,
    calibration = calibration,
    observed = observed,
    aod_validation = validation,
    lrt_elapsed_seconds = lrt_elapsed,
    config = cfg,
    session_info = utils::sessionInfo()
  )
  saveRDS(artifact, rds_path, compress = "xz")
  plot_calibration(calibration, figure_pdf, figure_png)

  alpha05 <- calibration[calibration$alpha == 0.05, , drop = FALSE]
  method_summary <- do.call(rbind, lapply(split(alpha05, alpha05$method), function(x) {
    data.frame(
      method = x$method[1L],
      mean_rate = mean(x$rejection_rate),
      min_rate = min(x$rejection_rate),
      max_rate = max(x$rejection_rate),
      mean_abs_error = mean(abs(x$calibration_error)),
      max_abs_error = max(abs(x$calibration_error)),
      mean_failure_rate = mean(x$failure_rate),
      stringsAsFactors = FALSE
    )
  }))
  cat("\nCalibration at alpha=.05:\n")
  print(method_summary, row.names = FALSE, digits = 5)
  cat(
    "\nScore-information ESS range: ", format(min(score_ess), digits = 5),
    " to ", format(max(score_ess), digits = 5), " (nominal n=", n, ")\n",
    sep = ""
  )
  cat("AOD validation:\n")
  print(validation)
  cat("LRT elapsed: ", format(lrt_elapsed, digits = 5), " seconds\n", sep = "")
  for (path in c(
    calibration_path, alpha05_path, observed_path, rds_path, figure_pdf, figure_png
  )) if (file.exists(path)) cat("Wrote ", path, "\n", sep = "")
  invisible(artifact)
}

if (sys.nframe() == 0L) run_main(cfg)
