#!/usr/bin/env Rscript

# Illustrative working-Poisson score example for manuscript/combine.tex.
#
# DGPs (with the manuscript's omitted marginal law made explicit):
#   Z ~ N(0, 1),
#   mu_Y(Z) = exp(2 Z),
#   (a) marginal independence: X ~ Bernoulli(expit(0)), phi = 1.5;
#   (b) conditional independence: X | Z ~ Bernoulli(expit(Z)),
#       phi = 1 (an exactly Poisson outcome);
#   (c) variance-robustness comparison with binary Z, a heteroskedastic
#       Gaussian exposure, and a multiplicative quasi-Poisson factor
#       c(Z) = 5;
#   (d) the same mean and exposure design with nonmultiplicative variance
#       factors c(0) = 5 and c(1) = 1.
# In the continuous-Z designs, E(Y | Z) = mu_Y(Z) and
# Var(Y | Z) = phi mu_Y(Z). In every conditional-independence design,
# X and Y are drawn independently given Z.
#
# A Poisson log-likelihood is used as the working model. It is correctly
# specified when phi = 1 and has a misspecified variance when phi = 1.5.
# The selected one-sided null tests are drawn from:
#   1. normal calibration of the working-information Poisson score;
#   2. permutation calibration of the unstudentized score;
#   3. permutation calibration of the information-studentized score;
#   4. permutation calibration of the empirically studentized score,
#      recomputing the residual-product variance for every permuted X.
# Figure 1 uses methods 1, 2, and 4 to isolate the need for studentization.
# The Section 4.2 figure uses methods 1, 3, and 4 to isolate variance
# robustness under multiplicative and nonmultiplicative misspecification.
#
# The output figure follows the spacrt-manuscript plotting conventions:
# theme_bw with uncluttered panels and a reverse-log p-value QQ plot using
# katlabutils.

required_packages <- c("ggplot2", "katlabutils", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing plotting packages: ", paste(missing_packages, collapse = ", "),
    ". Restore the repository's renv environment before running this script."
  )
}

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
validate_only <- "--validate-only" %in% args
plot_only <- "--plot-only" %in% args || validate_only
flag_options <- c("--quick", "--plot-only", "--validate-only")
value_options <- c(
  "n", "repetitions", "permutations", "cores", "seed", "alpha",
  "x-intercept", "conditional-x-slope", "y-intercept", "y-slope",
  "dispersion", "robustness-dispersion", "alternative", "output-dir"
)
option_names <- sub("=.*$", "", args)
if (any(!option_names %in% c(flag_options, paste0("--", value_options))) ||
    anyDuplicated(option_names) ||
    any(option_names %in% paste0("--", value_options) & !grepl("=", args))) {
  stop("Unknown, repeated, or malformed option. Use --name=value for numeric options.")
}

read_option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  matches <- startsWith(args, prefix)
  if (sum(matches) > 1L) stop("Option supplied more than once: --", name)
  if (!any(matches)) return(default)
  sub(prefix, "", args[matches], fixed = TRUE)
}

read_integer_option <- function(name, default, lower = 1L) {
  raw_value <- read_option(name, default)
  numeric_value <- suppressWarnings(as.numeric(raw_value))
  if (length(numeric_value) != 1L || !is.finite(numeric_value) ||
      numeric_value != floor(numeric_value) || numeric_value < lower ||
      numeric_value > .Machine$integer.max) {
    stop("--", name, " must be an integer >= ", lower, ".")
  }
  as.integer(numeric_value)
}

read_numeric_option <- function(name, default, lower, upper) {
  value <- suppressWarnings(as.numeric(read_option(name, default)))
  if (length(value) != 1L || !is.finite(value) ||
      value <= lower || value >= upper) {
    stop("--", name, " must lie strictly between ", lower, " and ", upper, ".")
  }
  value
}

read_finite_option <- function(name, default) {
  value <- suppressWarnings(as.numeric(read_option(name, default)))
  if (length(value) != 1L || !is.finite(value)) {
    stop("--", name, " must be a finite number.")
  }
  value
}

read_choice_option <- function(name, default, choices) {
  value <- read_option(name, default)
  if (length(value) != 1L || !value %in% choices) {
    stop(
      "--", name, " must be one of: ", paste(choices, collapse = ", "), "."
    )
  }
  value
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_argument)))
} else {
  normalizePath(getwd())
}
source(file.path(script_directory, "..", "reproduction-helpers.R"))
script_file <- file.path(script_directory, "glm-variance-robustness.R")
package_root <- reproduction_root(script_file)
results_root <- reproduction_results_root(package_root)
default_output_directory <- file.path(
  results_root,
  "simulation-results",
  "glm-variance-robustness"
)

default_repetitions <- if (quick) 100L else 2000L
default_permutations <- if (quick) 199L else 2000L
x_intercept <- read_finite_option("x-intercept", 0)
conditional_x_slope <- read_finite_option("conditional-x-slope", 1)
y_intercept <- read_finite_option("y-intercept", 0)
y_slope <- read_finite_option("y-slope", 2)
misspecified_dispersion <- read_numeric_option(
  "dispersion", 1.5, lower = 1, upper = 100
)

format_parameter <- function(value) {
  format(value, trim = TRUE, scientific = FALSE, digits = 8)
}

format_linear_predictor <- function(intercept, slope) {
  intercept_label <- format_parameter(intercept)
  if (slope == 0) return(intercept_label)
  magnitude <- abs(slope)
  coefficient_label <- if (magnitude == 1) {
    ""
  } else {
    paste0(format_parameter(magnitude), " ")
  }
  paste0(
    intercept_label,
    if (slope > 0) " + " else " - ",
    coefficient_label,
    "Z"
  )
}

x_intercept_label <- format_parameter(x_intercept)
x_linear_predictor_label <- format_linear_predictor(
  x_intercept,
  conditional_x_slope
)

robustness_dispersion <- read_numeric_option(
  "robustness-dispersion", 5, lower = 1, upper = 100
)

robustness_probability_z1 <- 3 / 4
robustness_x_sd_z0 <- 1
robustness_x_sd_z1 <- 1 / 3
robustness_y_mean_z0 <- 20
robustness_y_mean_z1 <- 30

scenario_specs <- data.frame(
  scenario = c(
    "marginal_independence",
    "conditional_correct",
    "conditional_multiplicative",
    "conditional_arbitrary"
  ),
  scenario_label = c(
    paste0("X ~ Bernoulli(expit(", x_intercept_label, "))"),
    paste0("X | Z ~ Bernoulli(expit(", x_linear_predictor_label, "))"),
    "Heteroskedastic Gaussian X; multiplicative variance",
    "Heteroskedastic Gaussian X; arbitrary variance"
  ),
  design_type = c("continuous", "continuous", "robustness", "robustness"),
  x_slope = c(0, conditional_x_slope, NA_real_, NA_real_),
  probability_z1 = c(
    NA_real_, NA_real_, robustness_probability_z1, robustness_probability_z1
  ),
  x_sd_z0 = c(
    NA_real_, NA_real_, robustness_x_sd_z0, robustness_x_sd_z0
  ),
  x_sd_z1 = c(
    NA_real_, NA_real_, robustness_x_sd_z1, robustness_x_sd_z1
  ),
  y_mean_z0 = c(
    NA_real_, NA_real_, robustness_y_mean_z0, robustness_y_mean_z0
  ),
  y_mean_z1 = c(
    NA_real_, NA_real_, robustness_y_mean_z1, robustness_y_mean_z1
  ),
  dispersion_z0 = c(
    misspecified_dispersion,
    1,
    robustness_dispersion,
    robustness_dispersion
  ),
  dispersion_z1 = c(
    misspecified_dispersion,
    1,
    robustness_dispersion,
    1
  ),
  variance_specification = c(
    "multiplicative quasi-Poisson variance",
    "correct Poisson variance",
    "multiplicative variance misspecification",
    "nonmultiplicative variance misspecification"
  ),
  stringsAsFactors = FALSE
)

config <- list(
  n = read_integer_option("n", 500L, lower = 20L),
  repetitions = read_integer_option(
    "repetitions", default_repetitions, lower = 10L
  ),
  permutations = read_integer_option(
    "permutations", default_permutations, lower = 19L
  ),
  cores = read_integer_option("cores", 1L, lower = 1L),
  seed = read_integer_option("seed", 20260805L, lower = 1L),
  alpha = read_numeric_option("alpha", 0.01, lower = 0, upper = 1),
  z_mean = 0,
  z_sd = 1,
  x_intercept = x_intercept,
  conditional_x_slope = conditional_x_slope,
  y_intercept = y_intercept,
  y_slope = y_slope,
  misspecified_dispersion = misspecified_dispersion,
  robustness_dispersion = robustness_dispersion,
  robustness_probability_z1 = robustness_probability_z1,
  robustness_x_sd_z0 = robustness_x_sd_z0,
  robustness_x_sd_z1 = robustness_x_sd_z1,
  robustness_y_mean_z0 = robustness_y_mean_z0,
  robustness_y_mean_z1 = robustness_y_mean_z1,
  simulation_schema_version = 2L,
  scenarios = scenario_specs,
  alternative = read_choice_option(
    "alternative", "greater", c("greater", "less")
  )
)

output_directory <- read_option("output-dir", default_output_directory)
if (!plot_only) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}
output_directory <- normalizePath(output_directory, mustWork = FALSE)

if (.Platform$OS.type == "windows" && config$cores > 1L) {
  warning("Parallel fork execution is unavailable on Windows; using one core.")
  config$cores <- 1L
}

mc_pvalue <- function(reference, observed, alternative) {
  exceedances <- switch(
    alternative,
    greater = sum(reference >= observed),
    less = sum(reference <= observed)
  )
  (1 + exceedances) / (length(reference) + 1)
}

wilson_interval <- function(rejections, repetitions, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  estimate <- rejections / repetitions
  denominator <- 1 + z^2 / repetitions
  center <- (estimate + z^2 / (2 * repetitions)) / denominator
  half_width <- z * sqrt(
    estimate * (1 - estimate) / repetitions +
      z^2 / (4 * repetitions^2)
  ) / denominator
  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

one_replicate_once <- function(config, scenario_key) {
  n <- config$n
  B <- config$permutations

  scenario_index <- match(scenario_key, config$scenarios$scenario)
  if (is.na(scenario_index)) stop("Unknown scenario: ", scenario_key)
  scenario <- config$scenarios[scenario_index, , drop = FALSE]

  if (scenario$design_type == "continuous") {
    z <- stats::rnorm(n, mean = config$z_mean, sd = config$z_sd)
    probability_x <- stats::plogis(
      config$x_intercept + scenario$x_slope * z
    )
    conditional_mean_y <- exp(config$y_intercept + config$y_slope * z)
  } else if (scenario$design_type == "robustness") {
    z <- stats::rbinom(n, size = 1L, prob = scenario$probability_z1)
    probability_x <- rep(0.5, n)
    conditional_mean_y <- ifelse(
      z == 0,
      scenario$y_mean_z0,
      scenario$y_mean_z1
    )
  } else {
    stop("Unsupported design type: ", scenario$design_type)
  }
  dispersion <- ifelse(
    z == 0,
    scenario$dispersion_z0,
    scenario$dispersion_z1
  )
  if (scenario$design_type == "continuous") {
    dispersion[] <- scenario$dispersion_z0
  }
  if (any(!is.finite(dispersion)) || any(dispersion < 1)) {
    stop("Scenario dispersions must be finite and at least one.")
  }

  # At dispersion one, draw from the correctly specified Poisson likelihood.
  # Above one, the observation-specific negative-binomial size gives
  # Var(Y | Z) = mu + mu^2 / size = dispersion * mu exactly. Separate draws
  # enforce X independent of Y conditional on Z.
  if (scenario$design_type == "continuous") {
    x <- stats::rbinom(n, size = 1L, prob = probability_x)
  } else {
    conditional_sd_x <- ifelse(
      z == 0,
      scenario$x_sd_z0,
      scenario$x_sd_z1
    )
    x <- stats::rnorm(n, mean = 0, sd = conditional_sd_x)
  }
  poisson_indices <- dispersion == 1
  y <- numeric(n)
  if (any(poisson_indices)) {
    y[poisson_indices] <- stats::rpois(
      sum(poisson_indices), lambda = conditional_mean_y[poisson_indices]
    )
  }
  if (any(!poisson_indices)) {
    negative_binomial_size <- conditional_mean_y[!poisson_indices] /
      (dispersion[!poisson_indices] - 1)
    y[!poisson_indices] <- stats::rnbinom(
      sum(!poisson_indices),
      size = negative_binomial_size,
      mu = conditional_mean_y[!poisson_indices]
    )
  }
  if (any(!is.finite(y)) || sum(y) == 0L) return(NULL)

  design <- cbind(intercept = 1, z = z)
  fit <- suppressWarnings(stats::glm.fit(
    x = design,
    y = y,
    family = stats::poisson(link = "log"),
    control = stats::glm.control(epsilon = 1e-10, maxit = 100L)
  ))
  fitted_mean <- as.numeric(fit$fitted.values)
  if (!isTRUE(fit$converged) || any(!is.finite(fitted_mean))) return(NULL)

  response_residual <- y - fitted_mean
  fisher_weight <- fitted_mean
  if (any(!is.finite(fisher_weight)) || sum(fisher_weight) <= 0) return(NULL)

  weighted_gram <- crossprod(design, fisher_weight * design)
  gram_cholesky <- tryCatch(chol(weighted_gram), error = function(e) NULL)
  if (is.null(gram_cholesky)) return(NULL)

  permutation_indices <- replicate(B, sample.int(n), simplify = "matrix")
  x_bank <- cbind(
    observed = x,
    matrix(x[permutation_indices], nrow = n, ncol = B)
  )

  # Profile information, vectorized over observed X and all B permutations:
  #   X' W X - X' W D (D' W D)^(-1) D' W X.
  weighted_cross <- crossprod(design, fisher_weight * x_bank)
  projection_coefficients <- backsolve(
    gram_cholesky,
    forwardsolve(t(gram_cholesky), weighted_cross)
  )
  raw_profile_information <-
    colSums(fisher_weight * x_bank^2) -
    colSums(weighted_cross * projection_coefficients)
  information_per_observation <- raw_profile_information / n
  score_u <- as.vector(crossprod(response_residual, x_bank)) / sqrt(n)
  residualized_x_bank <- x_bank - design %*% projection_coefficients
  score_product_bank <- residualized_x_bank * response_residual
  score_product_means <- colMeans(score_product_bank)
  empirical_variance_per_observation <-
    colMeans(score_product_bank^2) - score_product_means^2
  if (any(empirical_variance_per_observation < -1e-10)) return(NULL)
  empirical_variance_per_observation <- pmax(
    empirical_variance_per_observation,
    0
  )
  positive_empirical_variance <- empirical_variance_per_observation > 1e-14
  constant_x <- all(x == x[1L])
  if (constant_x) {
    # With no exposure variation the score direction is unidentified. The
    # natural nonrejection convention is T = 0; all permutation statistics
    # are tied and therefore have plus-one p-value 1. Keeping the dataset is
    # important for finite-sample Type-I error rather than redrawing it.
    information_per_observation[] <- 0
    empirical_variance_per_observation[] <- 0
    score_u[] <- 0
    score_t_information <- rep(0, B + 1L)
    score_t_empirical <- rep(0, B + 1L)
  } else {
    if (any(!is.finite(information_per_observation)) ||
        any(information_per_observation <= 1e-14)) {
      return(NULL)
    }
    if (any(!is.finite(empirical_variance_per_observation))) return(NULL)
    score_t_information <- score_u / sqrt(information_per_observation)
    score_t_empirical <- numeric(B + 1L)
    score_t_empirical[positive_empirical_variance] <-
      score_u[positive_empirical_variance] /
      sqrt(empirical_variance_per_observation[positive_empirical_variance])
  }
  if (any(!is.finite(score_u)) ||
      any(!is.finite(score_t_information)) ||
      any(!is.finite(score_t_empirical))) return(NULL)

  observed <- 1L
  reference <- 2L:(B + 1L)
  p_normal_greater <- unname(stats::pnorm(
    score_t_information[observed], lower.tail = FALSE
  ))
  p_normal_less <- unname(stats::pnorm(
    score_t_information[observed], lower.tail = TRUE
  ))
  p_permutation_unstudentized_greater <- mc_pvalue(
    score_u[reference], score_u[observed], alternative = "greater"
  )
  p_permutation_unstudentized_less <- mc_pvalue(
    score_u[reference], score_u[observed], alternative = "less"
  )
  p_permutation_information_greater <- mc_pvalue(
    score_t_information[reference],
    score_t_information[observed],
    alternative = "greater"
  )
  p_permutation_information_less <- mc_pvalue(
    score_t_information[reference],
    score_t_information[observed],
    alternative = "less"
  )
  p_permutation_empirical_greater <- mc_pvalue(
    score_t_empirical[reference],
    score_t_empirical[observed],
    alternative = "greater"
  )
  p_permutation_empirical_less <- mc_pvalue(
    score_t_empirical[reference],
    score_t_empirical[observed],
    alternative = "less"
  )
  selected_suffix <- config$alternative
  p_normal <- get(paste0("p_normal_", selected_suffix))
  p_permutation_unstudentized <- get(paste0(
    "p_permutation_unstudentized_", selected_suffix
  ))
  p_permutation_information <- get(paste0(
    "p_permutation_information_", selected_suffix
  ))
  p_permutation_empirical <- get(paste0(
    "p_permutation_empirical_", selected_suffix
  ))

  observed_residualized_x <- residualized_x_bank[, observed]
  residualized_score_u <- if (constant_x) {
    0
  } else {
    sum(observed_residualized_x * response_residual) / sqrt(n)
  }

  c(
    p_normal = p_normal,
    p_permutation_unstudentized = p_permutation_unstudentized,
    p_permutation_information = p_permutation_information,
    p_permutation_empirical = p_permutation_empirical,
    p_normal_greater = p_normal_greater,
    p_normal_less = p_normal_less,
    p_permutation_unstudentized_greater =
      p_permutation_unstudentized_greater,
    p_permutation_unstudentized_less = p_permutation_unstudentized_less,
    p_permutation_information_greater =
      p_permutation_information_greater,
    p_permutation_information_less = p_permutation_information_less,
    p_permutation_empirical_greater = p_permutation_empirical_greater,
    p_permutation_empirical_less = p_permutation_empirical_less,
    u_observed = unname(score_u[observed]),
    t_information_observed = unname(score_t_information[observed]),
    t_empirical_observed = unname(score_t_empirical[observed]),
    u_permutation_mean = mean(score_u[reference]),
    u_permutation_sd = stats::sd(score_u[reference]),
    t_information_permutation_mean = mean(score_t_information[reference]),
    t_information_permutation_sd = stats::sd(score_t_information[reference]),
    t_empirical_permutation_mean = mean(score_t_empirical[reference]),
    t_empirical_permutation_sd = stats::sd(score_t_empirical[reference]),
    u_permutation_q99 = unname(stats::quantile(
      score_u[reference], probs = 1 - config$alpha, type = 1
    )),
    t_information_permutation_q99 = unname(stats::quantile(
      score_t_information[reference], probs = 1 - config$alpha, type = 1
    )),
    u_permutation_q01 = unname(stats::quantile(
      score_u[reference], probs = config$alpha, type = 1
    )),
    t_information_permutation_q01 = unname(stats::quantile(
      score_t_information[reference], probs = config$alpha, type = 1
    )),
    t_empirical_permutation_q99 = unname(stats::quantile(
      score_t_empirical[reference], probs = 1 - config$alpha, type = 1
    )),
    t_empirical_permutation_q01 = unname(stats::quantile(
      score_t_empirical[reference], probs = config$alpha, type = 1
    )),
    profile_information_observed = unname(
      information_per_observation[observed]
    ),
    empirical_variance_observed = unname(
      empirical_variance_per_observation[observed]
    ),
    zero_empirical_permutation_studentizers = sum(
      !positive_empirical_variance[reference]
    ),
    number_x_positive = sum(x > 0),
    total_y_count = sum(y),
    number_y_positive = sum(y > 0),
    exposed_y_count = sum(y[x > 0]),
    number_xy_positive = sum(x > 0 & y > 0),
    mean_probability_x = mean(probability_x),
    mean_conditional_mean_y = mean(conditional_mean_y),
    pearson_dispersion_estimate = sum(
      response_residual^2 / fitted_mean
    ) / (n - ncol(design)),
    score_equation_error = max(abs(crossprod(design, response_residual))),
    numerator_equivalence_error = abs(score_u[observed] - residualized_score_u),
    permutation_carrier_count_error = max(abs(colSums(x_bank) - sum(x))),
    glm_iterations = fit$iter,
    constant_x = as.integer(constant_x)
  )
}

one_replicate <- function(seed, config, scenario_key) {
  set.seed(seed)
  failure_message <- NULL
  answer <- tryCatch(
    one_replicate_once(
      config,
      scenario_key = scenario_key
    ),
    error = function(e) {
      failure_message <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(answer) || !all(is.finite(answer))) {
    detail <- if (is.null(failure_message)) {
      "degenerate data, a failed null fit, or an invalid denominator"
    } else {
      failure_message
    }
    stop(
      "Replicate seed ", seed, " failed: ", detail,
      ". The dataset was not redrawn."
    )
  }
  c(attempts = 1, answer)
}

summarize_pvalues <- function(replicate_results, alpha) {
  method_keys <- c(
    "normal",
    "permutation_unstudentized",
    "permutation_information",
    "permutation_empirical"
  )
  method_labels <- c(
    "Poisson normal cutoff",
    "Unstudentized permutation",
    "Barry et al.",
    "Empirical-studentized permutation"
  )
  pvalue_columns <- paste0("p_", method_keys)
  scenario_keys <- unique(replicate_results$scenario)

  rows <- unlist(lapply(scenario_keys, function(scenario_key) {
    scenario_results <- replicate_results[
      replicate_results$scenario == scenario_key,
      ,
      drop = FALSE
    ]
    repetitions <- nrow(scenario_results)
    lapply(seq_along(method_keys), function(index) {
      pvalues <- scenario_results[[pvalue_columns[index]]]
      rejections <- sum(pvalues <= alpha)
      interval <- wilson_interval(rejections, repetitions)
      data.frame(
        scenario = scenario_key,
        scenario_label = scenario_results$scenario_label[1L],
        x_slope = scenario_results$x_slope[1L],
        dispersion_z0 = scenario_results$dispersion_z0[1L],
        dispersion_z1 = scenario_results$dispersion_z1[1L],
        variance_specification =
          scenario_results$variance_specification[1L],
        method = method_keys[index],
        method_label = method_labels[index],
        alpha = alpha,
        rejections = rejections,
        repetitions = repetitions,
        rejection_rate = rejections / repetitions,
        mc_standard_error = sqrt(
          (rejections / repetitions) * (1 - rejections / repetitions) /
            repetitions
        ),
        ci_lower = unname(interval["lower"]),
        ci_upper = unname(interval["upper"]),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  do.call(rbind, rows)
}

summarize_scale_diagnostics <- function(replicate_results) {
  scenario_keys <- unique(replicate_results$scenario)
  rows <- lapply(scenario_keys, function(scenario_key) {
    scenario_results <- replicate_results[
      replicate_results$scenario == scenario_key,
      ,
      drop = FALSE
    ]
    data.frame(
      scenario = scenario_key,
      scenario_label = scenario_results$scenario_label[1L],
      dispersion_z0 = rep(scenario_results$dispersion_z0[1L], 3L),
      dispersion_z1 = rep(scenario_results$dispersion_z1[1L], 3L),
      variance_specification = rep(
        scenario_results$variance_specification[1L],
        3L
      ),
      statistic = c(
        "U (unstudentized)",
        "T (information-studentized)",
        "T (empirically studentized)"
      ),
      observed_sampling_mean = c(
        mean(scenario_results$u_observed),
        mean(scenario_results$t_information_observed),
        mean(scenario_results$t_empirical_observed)
      ),
      observed_sampling_sd = c(
        stats::sd(scenario_results$u_observed),
        stats::sd(scenario_results$t_information_observed),
        stats::sd(scenario_results$t_empirical_observed)
      ),
      mean_conditional_permutation_mean = c(
        mean(scenario_results$u_permutation_mean),
        mean(scenario_results$t_information_permutation_mean),
        mean(scenario_results$t_empirical_permutation_mean)
      ),
      root_mean_conditional_permutation_variance = c(
        sqrt(mean(scenario_results$u_permutation_sd^2)),
        sqrt(mean(scenario_results$t_information_permutation_sd^2)),
        sqrt(mean(scenario_results$t_empirical_permutation_sd^2))
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_dgp_diagnostics <- function(replicate_results) {
  scenario_keys <- unique(replicate_results$scenario)
  rows <- lapply(scenario_keys, function(scenario_key) {
    scenario_results <- replicate_results[
      replicate_results$scenario == scenario_key,
      ,
      drop = FALSE
    ]
    positive_information <-
      scenario_results$profile_information_observed > 0
    mean_empirical_to_information_ratio <- if (any(positive_information)) {
      mean(
        scenario_results$empirical_variance_observed[positive_information] /
          scenario_results$profile_information_observed[positive_information]
      )
    } else {
      NA_real_
    }
    data.frame(
      scenario = scenario_key,
      scenario_label = scenario_results$scenario_label[1L],
      dispersion_z0 = scenario_results$dispersion_z0[1L],
      dispersion_z1 = scenario_results$dispersion_z1[1L],
      variance_specification =
        scenario_results$variance_specification[1L],
      repetitions = nrow(scenario_results),
      mean_number_x_positive = mean(scenario_results$number_x_positive),
      mean_total_y_count = mean(scenario_results$total_y_count),
      mean_number_y_positive = mean(scenario_results$number_y_positive),
      mean_number_xy_positive = mean(scenario_results$number_xy_positive),
      mean_exposed_y_count = mean(scenario_results$exposed_y_count),
      mean_pearson_dispersion = mean(
        scenario_results$pearson_dispersion_estimate
      ),
      mean_empirical_to_information_ratio =
        mean_empirical_to_information_ratio,
      zero_empirical_permutation_studentizers = sum(
        scenario_results$zero_empirical_permutation_studentizers
      ),
      constant_x_datasets = sum(scenario_results$constant_x),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_figure <- function(replicate_results, summary, config, output_directory) {
  required_scenarios <- c(
    "marginal_independence",
    "conditional_correct",
    "conditional_multiplicative",
    "conditional_arbitrary"
  )
  missing_scenarios <- setdiff(
    required_scenarios,
    unique(replicate_results$scenario)
  )
  if (length(missing_scenarios) > 0L) {
    stop(
      "The saved results predate the four-scenario design. Missing: ",
      paste(missing_scenarios, collapse = ", "),
      ". Rerun the simulation without --plot-only."
    )
  }

  all_method_keys <- c(
    "normal",
    "permutation_unstudentized",
    "permutation_information",
    "permutation_empirical"
  )
  all_method_labels <- c(
    "Poisson normal cutoff",
    "Unstudentized permutation",
    "Barry et al.",
    "Empirical-studentized permutation"
  )
  method_colors <- c(
    normal = "#F8766D",
    permutation_unstudentized = "#00BA38",
    permutation_information = "#C77CFF",
    permutation_empirical = "#619CFF"
  )
  method_labels <- stats::setNames(all_method_labels, all_method_keys)
  method_color_labels <- stats::setNames(
    method_colors,
    all_method_labels
  )
  summary$method <- factor(summary$method, levels = all_method_keys)

  theme_spacrt <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
      axis.title = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8),
      plot.tag = ggplot2::element_text(face = "bold", size = 11),
      plot.tag.position = c(0, 1)
    )

  bar_upper_for <- function(scenario_keys, selected_method_keys) {
    scenario_summary <- summary[
      summary$scenario %in% scenario_keys &
        as.character(summary$method) %in% selected_method_keys,
      ,
      drop = FALSE
    ]
    max(scenario_summary$ci_upper, config$alpha) * 1.25
  }
  make_bar_plot <- function(
    scenario_key,
    selected_method_keys,
    bar_upper,
    show_legend = FALSE
  ) {
    scenario_summary <- summary[
      summary$scenario == scenario_key &
        as.character(summary$method) %in% selected_method_keys,
      ,
      drop = FALSE
    ]
    scenario_summary$method <- factor(
      as.character(scenario_summary$method),
      levels = selected_method_keys
    )
    ggplot2::ggplot(
      scenario_summary,
      ggplot2::aes(x = method, y = rejection_rate, fill = method)
    ) +
      ggplot2::geom_col(width = 0.64) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        width = 0.13,
        linewidth = 0.55
      ) +
      ggplot2::geom_hline(
        yintercept = config$alpha,
        linetype = "dashed",
        colour = "black",
        linewidth = 0.55
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          y = ci_upper,
          label = sprintf("%.3f", rejection_rate)
        ),
        vjust = -0.5,
        size = 2.5
      ) +
      ggplot2::scale_fill_manual(
        values = method_colors[selected_method_keys],
        breaks = selected_method_keys,
        labels = unname(method_labels[selected_method_keys]),
        guide = if (show_legend) "legend" else "none"
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, bar_upper),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = "Methods",
        y = "Rejection rate"
      ) +
      theme_spacrt +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.width = grid::unit(0.18, "in"),
        plot.margin = ggplot2::margin(3, 2, 8, 2)
      )
  }

  qq_floor <- 1 / (config$permutations + 1)
  qq_minimum_exponent <- floor(log10(qq_floor))
  qq_exponents <- seq(0, qq_minimum_exponent, by = -1)
  qq_breaks <- 10^qq_exponents
  qq_labels <- as.expression(lapply(qq_exponents, function(exponent) {
    if (exponent == 0) 1 else bquote(10^.(exponent))
  }))
  # Extend the observed-p-value range slightly past the permutation floor so
  # points clamped at qq_floor remain fully visible instead of being clipped.
  qq_panel_floor <- qq_floor / 10^0.15

  make_qq_plot <- function(scenario_key, selected_method_keys) {
    scenario_results <- replicate_results[
      replicate_results$scenario == scenario_key,
      ,
      drop = FALSE
    ]
    qq_data <- data.frame(
      method = rep(selected_method_keys, each = nrow(scenario_results)),
      p_value = unlist(lapply(
        paste0("p_", selected_method_keys),
        function(column_name) scenario_results[[column_name]]
      ), use.names = FALSE),
      stringsAsFactors = FALSE
    )
    qq_data$method <- factor(
      qq_data$method,
      levels = selected_method_keys,
      labels = unname(method_labels[selected_method_keys])
    )
    qq_band_data <- qq_data[
      qq_data$method == unname(method_labels[selected_method_keys[1L]]),
      ,
      drop = FALSE
    ]

    ggplot2::ggplot(
      qq_data,
      ggplot2::aes(y = p_value, colour = method)
    ) +
      katlabutils::stat_qq_band(
        data = qq_band_data,
        mapping = ggplot2::aes(y = p_value),
        inherit.aes = FALSE,
        fill = "grey75",
        colour = NA
      ) +
      katlabutils::stat_qq_points(
        ymin = qq_floor,
        size = 1.25,
        max_pts_to_plot = 1000
      ) +
      ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        colour = "black",
        linewidth = 0.55
      ) +
      ggplot2::scale_colour_manual(
        values = method_color_labels[
          unname(method_labels[selected_method_keys])
        ]
      ) +
      ggplot2::scale_x_continuous(
        trans = katlabutils::revlog_trans(10),
        breaks = qq_breaks,
        labels = qq_labels
      ) +
      ggplot2::scale_y_continuous(
        trans = katlabutils::revlog_trans(10),
        breaks = qq_breaks,
        labels = qq_labels
      ) +
      ggplot2::coord_cartesian(
        xlim = c(qq_floor, 1),
        ylim = c(qq_panel_floor, 1),
        expand = FALSE
      ) +
      ggplot2::labs(
        x = "Expected",
        y = "Observed"
      ) +
      theme_spacrt +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.width = grid::unit(0.18, "in"),
        plot.margin = ggplot2::margin(3, 2, 8, 2)
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          nrow = 1,
          override.aes = list(size = 2)
        )
      )
  }

  figure_one_methods <- c(
    "normal",
    "permutation_unstudentized",
    "permutation_empirical"
  )
  variance_robustness_methods <- c(
    "normal",
    "permutation_information",
    "permutation_empirical"
  )
  figure_one_bar_upper <- bar_upper_for(
    c("marginal_independence", "conditional_correct"),
    figure_one_methods
  )
  variance_bar_upper <- bar_upper_for(
    c("conditional_multiplicative", "conditional_arbitrary"),
    variance_robustness_methods
  )

  marginal_pair <- list(
    bar = make_bar_plot(
      "marginal_independence",
      figure_one_methods,
      figure_one_bar_upper
    ) +
      ggplot2::labs(tag = "a"),
    qq = make_qq_plot("marginal_independence", figure_one_methods) +
      ggplot2::labs(tag = "b")
  )
  conditional_pair <- list(
    bar = make_bar_plot(
      "conditional_correct",
      figure_one_methods,
      figure_one_bar_upper
    ) +
      ggplot2::labs(tag = "c"),
    qq = make_qq_plot("conditional_correct", figure_one_methods) +
      ggplot2::labs(tag = "d")
  )
  text_panel <- function(label, font_size = 12, font_face = "plain") {
    patchwork::wrap_elements(full = grid::textGrob(
      label,
      gp = grid::gpar(fontsize = font_size, fontface = font_face)
    ))
  }
  marginal_header <- text_panel(
    "Marginal independence",
    font_size = 10.5,
    font_face = "bold"
  )
  conditional_header <- text_panel(
    "Conditional independence",
    font_size = 10.5,
    font_face = "bold"
  )
  combined <- patchwork::wrap_plots(
    A = marginal_header,
    B = conditional_header,
    S = patchwork::plot_spacer(),
    G = marginal_pair$bar,
    H = marginal_pair$qq,
    I = conditional_pair$bar,
    J = conditional_pair$qq,
    K = patchwork::guide_area(),
    design = "
      AABB
      SSSS
      GHIJ
      KKKK
    ",
    widths = rep(1, 4),
    heights = c(0.12, 0.06, 1, 0.20),
    guides = "collect"
  )

  pdf_path <- file.path(output_directory, "illustrative-example.pdf")
  png_path <- file.path(output_directory, "illustrative-example.png")
  plot_seed <- as.integer(
    (as.double(config$seed) + 104729) %% .Machine$integer.max
  )
  set.seed(plot_seed)
  ggplot2::ggsave(
    filename = pdf_path,
    plot = combined,
    device = grDevices::cairo_pdf,
    width = 7.2,
    height = 2.85,
    units = "in"
  )
  set.seed(plot_seed)
  ggplot2::ggsave(
    filename = png_path,
    plot = combined,
    width = 7.2,
    height = 2.85,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  multiplicative_pair <- list(
    bar = make_bar_plot(
      "conditional_multiplicative",
      variance_robustness_methods,
      variance_bar_upper
    ) + ggplot2::labs(tag = "a"),
    qq = make_qq_plot(
      "conditional_multiplicative",
      variance_robustness_methods
    ) + ggplot2::labs(tag = "b")
  )
  arbitrary_pair <- list(
    bar = make_bar_plot(
      "conditional_arbitrary",
      variance_robustness_methods,
      variance_bar_upper
    ) + ggplot2::labs(tag = "c"),
    qq = make_qq_plot(
      "conditional_arbitrary",
      variance_robustness_methods
    ) + ggplot2::labs(tag = "d")
  )
  multiplicative_header <- text_panel(
    "Multiplicative misspecification",
    font_size = 10,
    font_face = "bold"
  )
  arbitrary_header <- text_panel(
    "Arbitrary misspecification",
    font_size = 10,
    font_face = "bold"
  )
  misspecified_combined <- patchwork::wrap_plots(
    A = multiplicative_header,
    B = arbitrary_header,
    S = patchwork::plot_spacer(),
    G = multiplicative_pair$bar,
    H = multiplicative_pair$qq,
    I = arbitrary_pair$bar,
    J = arbitrary_pair$qq,
    K = patchwork::guide_area(),
    design = "
      AA#BB
      SSSSS
      GH#IJ
      KKKKK
    ",
    widths = c(1, 1, 0.15, 1, 1),
    heights = c(0.18, 0.06, 1, 0.20),
    guides = "collect"
  )
  misspecified_pdf_path <- file.path(
    output_directory,
    "illustrative-example-ci-misspecified.pdf"
  )
  misspecified_png_path <- file.path(
    output_directory,
    "illustrative-example-ci-misspecified.png"
  )
  misspecified_plot_seed <- as.integer(
    (as.double(plot_seed) + 1) %% .Machine$integer.max
  )
  set.seed(misspecified_plot_seed)
  ggplot2::ggsave(
    filename = misspecified_pdf_path,
    plot = misspecified_combined,
    device = grDevices::cairo_pdf,
    width = 7.2,
    height = 2.85,
    units = "in"
  )
  set.seed(misspecified_plot_seed)
  ggplot2::ggsave(
    filename = misspecified_png_path,
    plot = misspecified_combined,
    width = 7.2,
    height = 2.85,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  invisible(c(
    main_pdf = pdf_path,
    main_png = png_path,
    misspecified_pdf = misspecified_pdf_path,
    misspecified_png = misspecified_png_path
  ))
}

results_path <- Sys.getenv("PCIT_GLM_RESULTS_FILE", unset = file.path(
  output_directory, "illustrative-example-results.rds"
))

if (plot_only) {
  if (!file.exists(results_path)) {
    stop("--plot-only requested, but no results file exists at ", results_path)
  }
  saved <- readRDS(results_path)
  saved_schema_version <- saved$config$simulation_schema_version
  if (is.null(saved_schema_version) ||
      !identical(saved_schema_version, config$simulation_schema_version)) {
    stop(
      "The saved results use an outdated simulation design. Rerun without ",
      "--plot-only."
    )
  }
  canonical_config <- config
  canonical_fields <- c(
    "n", "repetitions", "permutations", "seed", "alpha",
    "z_mean", "z_sd", "x_intercept", "conditional_x_slope",
    "y_intercept", "y_slope", "misspecified_dispersion",
    "robustness_dispersion", "robustness_probability_z1",
    "robustness_x_sd_z0", "robustness_x_sd_z1",
    "robustness_y_mean_z0", "robustness_y_mean_z1",
    "simulation_schema_version", "scenarios", "alternative"
  )
  noncanonical_fields <- canonical_fields[!vapply(
    canonical_fields,
    function(field) {
      isTRUE(all.equal(
        saved$config[[field]], canonical_config[[field]],
        tolerance = 1e-12, check.attributes = TRUE
      ))
    },
    logical(1L)
  )]
  if (length(noncanonical_fields) > 0L) {
    stop(
      "The cached GLM result is a pilot or uses a noncanonical design ",
      "and cannot generate the manuscript figure. Mismatched fields: ",
      paste(noncanonical_fields, collapse = ", "), "."
    )
  }
  config <- saved$config
  replicate_results <- saved$replicate_results
  required_result_columns <- c(
    "p_permutation_information",
    "p_permutation_empirical",
    "t_information_observed",
    "t_empirical_observed"
  )
  missing_result_columns <- setdiff(
    required_result_columns,
    names(replicate_results)
  )
  if (length(missing_result_columns) > 0L) {
    stop(
      "The saved results predate empirical studentization. Rerun without ",
      "--plot-only. Missing columns: ",
      paste(missing_result_columns, collapse = ", "), "."
    )
  }
  pvalue_columns <- c(
    "p_normal", "p_permutation_unstudentized",
    "p_permutation_information", "p_permutation_empirical"
  )
  expected_keys <- unlist(lapply(config$scenarios$scenario, function(scenario) {
    paste(scenario, seq_len(config$repetitions), sep = ":")
  }), use.names = FALSE)
  if (!is.data.frame(replicate_results) ||
      !all(c("scenario", "replicate", pvalue_columns) %in%
           names(replicate_results))) {
    stop("Saved GLM results are missing replicate identifiers or p-values.")
  }
  actual_keys <- paste(replicate_results$scenario, replicate_results$replicate,
                       sep = ":")
  pvalues <- as.matrix(replicate_results[pvalue_columns])
  if (length(actual_keys) != length(expected_keys) ||
      anyDuplicated(actual_keys) || !setequal(actual_keys, expected_keys) ||
      any(!is.finite(pvalues)) || any(pvalues < 0 | pvalues > 1)) {
    stop("Saved GLM results have incomplete/duplicate replicates or invalid p-values.")
  }
  recomputed_summary <- summarize_pvalues(replicate_results, config$alpha)
  # Labels may change editorially without changing the simulation. Validate
  # identifiers and numerical summaries, then render with the current labels.
  summary_fields <- setdiff(names(recomputed_summary),
                            c("scenario_label", "method_label"))
  if (!is.data.frame(saved$summary) ||
      !all(summary_fields %in% names(saved$summary)) ||
      !isTRUE(all.equal(saved$summary[summary_fields],
                       recomputed_summary[summary_fields],
                       tolerance = 1e-12, check.attributes = FALSE))) {
    stop("Saved GLM summary disagrees with the replicate results.")
  }
  summary <- recomputed_summary
  scale_diagnostics <- saved$scale_diagnostics
  dgp_diagnostics <- saved$dgp_diagnostics
  message("Loaded existing results from ", results_path)
} else {
  if (config$alpha * (config$permutations + 1L) < 2) {
    warning(
      "Permutation resolution is coarse at alpha=", config$alpha,
      "; use more permutations for inferential runs."
    )
  }

  job_grid <- do.call(rbind, lapply(
    seq_len(nrow(config$scenarios)),
    function(scenario_index) {
      data.frame(
        scenario = config$scenarios$scenario[scenario_index],
        scenario_label = config$scenarios$scenario_label[scenario_index],
        design_type = config$scenarios$design_type[scenario_index],
        x_slope = config$scenarios$x_slope[scenario_index],
        probability_z1 = config$scenarios$probability_z1[scenario_index],
        x_sd_z0 = config$scenarios$x_sd_z0[scenario_index],
        x_sd_z1 = config$scenarios$x_sd_z1[scenario_index],
        y_mean_z0 = config$scenarios$y_mean_z0[scenario_index],
        y_mean_z1 = config$scenarios$y_mean_z1[scenario_index],
        dispersion_z0 = config$scenarios$dispersion_z0[scenario_index],
        dispersion_z1 = config$scenarios$dispersion_z1[scenario_index],
        variance_specification =
          config$scenarios$variance_specification[scenario_index],
        replicate = seq_len(config$repetitions),
        stringsAsFactors = FALSE
      )
    }
  ))
  scenario_seed_offsets <- stats::setNames(
    c(0, 104729, 209759, 314573),
    c(
      "marginal_independence",
      "conditional_correct",
      "conditional_multiplicative",
      "conditional_arbitrary"
    )
  )
  scenario_seed_values <- vapply(
    config$scenarios$scenario,
    function(scenario_key) {
      seed_value <- as.integer(
        (as.double(config$seed) + scenario_seed_offsets[[scenario_key]]) %%
          .Machine$integer.max
      )
      if (seed_value == 0L) 1L else seed_value
    },
    integer(1L)
  )
  replicate_seeds <- unlist(lapply(
    scenario_seed_values,
    function(scenario_seed) {
      set.seed(scenario_seed)
      sample.int(
        .Machine$integer.max,
        size = config$repetitions,
        replace = FALSE
      )
    }
  ), use.names = FALSE)
  job_grid$replicate_seed <- replicate_seeds

  run_job <- function(job_index) {
    one_replicate(
      seed = job_grid$replicate_seed[job_index],
      config = config,
      scenario_key = job_grid$scenario[job_index]
    )
  }

  message(
    "Running n=", config$n,
    ", scenarios=", nrow(config$scenarios),
    ", repetitions per scenario=", config$repetitions,
    ", total datasets=", nrow(job_grid),
    ", permutations=", config$permutations,
    ", alpha=", config$alpha,
    ", alternative=", config$alternative,
    ", cores=", config$cores
  )
  start_time <- proc.time()["elapsed"]
  if (config$cores > 1L) {
    replicate_list <- parallel::mclapply(
      seq_len(nrow(job_grid)),
      run_job,
      mc.cores = config$cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    replicate_list <- lapply(seq_len(nrow(job_grid)), run_job)
  }
  elapsed_seconds <- unname(proc.time()["elapsed"] - start_time)
  failed_jobs <- which(vapply(
    replicate_list,
    inherits,
    logical(1L),
    what = "try-error"
  ))
  if (length(failed_jobs) > 0L) {
    first_failure <- failed_jobs[1L]
    stop(
      "Parallel simulation failed for job ", first_failure,
      " (scenario ", job_grid$scenario[first_failure],
      ", replicate ", job_grid$replicate[first_failure], "): ",
      as.character(replicate_list[[first_failure]])
    )
  }
  replicate_results <- cbind(
    job_grid,
    as.data.frame(do.call(rbind, replicate_list)),
    stringsAsFactors = FALSE
  )

  pvalue_columns <- c(
    "p_normal",
    "p_permutation_unstudentized",
    "p_permutation_information",
    "p_permutation_empirical"
  )
  if (any(!is.finite(as.matrix(replicate_results[pvalue_columns]))) ||
      any(as.matrix(replicate_results[pvalue_columns]) < 0) ||
      any(as.matrix(replicate_results[pvalue_columns]) > 1)) {
    stop("Internal check failed: invalid p-value detected.")
  }
  if (max(replicate_results$score_equation_error) > 1e-5) {
    stop("Internal check failed: null Poisson score equation is inaccurate.")
  }
  if (max(replicate_results$numerator_equivalence_error) > 1e-8) {
    stop("Internal check failed: raw and residualized score numerators differ.")
  }
  if (max(replicate_results$permutation_carrier_count_error) > 1e-10) {
    stop("Internal check failed: a permutation changed the exposure total.")
  }

  summary <- summarize_pvalues(replicate_results, config$alpha)
  scale_diagnostics <- summarize_scale_diagnostics(replicate_results)
  dgp_diagnostics <- summarize_dgp_diagnostics(replicate_results)
  saved <- list(
    config = config,
    elapsed_seconds = elapsed_seconds,
    job_grid = job_grid,
    scenario_seed_values = scenario_seed_values,
    replicate_seeds = replicate_seeds,
    replicate_results = replicate_results,
    summary = summary,
    scale_diagnostics = scale_diagnostics,
    dgp_diagnostics = dgp_diagnostics,
    style_reference = "https://github.com/Katsevich-Lab/spacrt-manuscript",
    session_info = utils::sessionInfo(),
    provenance = list(
      code_md5 = unname(tools::md5sum(script_file)),
      lockfile_md5 = unname(tools::md5sum(file.path(package_root, "renv.lock")))
    )
  )
  temporary_result <- tempfile(".glm-result-", tmpdir = dirname(results_path))
  saveRDS(saved, temporary_result)
  if (!file.rename(temporary_result, results_path)) {
    unlink(temporary_result)
    stop("Could not atomically publish GLM result: ", results_path)
  }
  utils::write.csv(
    replicate_results,
    file.path(output_directory, "replicate-results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(output_directory, "type-I-error-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    scale_diagnostics,
    file.path(output_directory, "scale-diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    dgp_diagnostics,
    file.path(output_directory, "dgp-diagnostics.csv"),
    row.names = FALSE
  )
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(output_directory, "session-info.txt")
  )
  message(sprintf("Simulation completed in %.1f seconds.", elapsed_seconds))
}

if (validate_only) {
  cat("Validated canonical GLM result: ", nrow(replicate_results),
      " complete scenario/replicate outputs.\n", sep = "")
  quit(save = "no", status = 0L)
}

figure_paths <- make_figure(
  replicate_results = replicate_results,
  summary = summary,
  config = config,
  output_directory = output_directory
)

print(summary, row.names = FALSE)
cat("\nScale diagnostics:\n")
print(scale_diagnostics, row.names = FALSE)
if (!is.null(dgp_diagnostics)) {
  cat("\nDGP diagnostics:\n")
  print(dgp_diagnostics, row.names = FALSE)
}
cat("\nSaved results: ", results_path, "\n", sep = "")
cat("Saved main figure:  ", unname(figure_paths["main_pdf"]), "\n", sep = "")
cat("Saved main preview: ", unname(figure_paths["main_png"]), "\n", sep = "")
cat(
  "Saved misspecified figure:  ",
  unname(figure_paths["misspecified_pdf"]),
  "\n",
  sep = ""
)
cat(
  "Saved misspecified preview: ",
  unname(figure_paths["misspecified_png"]),
  "\n",
  sep = ""
)
