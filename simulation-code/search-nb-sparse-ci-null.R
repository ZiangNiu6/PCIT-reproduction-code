#!/usr/bin/env Rscript

# Grid search for Type-I error in the sparse conditional-independence NB setup.
# The DGP and NB score implementations come from the formal Section 5.1
# specifier.  A Poisson-only test uses SCEPTRE's initial glm.fit step without
# NB size.  A normal-reference test uses our observed empirical-studentized
# NB score but substitutes 2 * Phi(-abs(T)) for permutation calibration.
# Within a run, a grid varies the true constant NB size theta, the outcome's
# coefficient beta_z, and the exposure's coefficient x_z_slope. Sample size
# and marginal exposure prevalence can also be set on the command line. NB
# variance is mu + mu^2 / theta, so smaller theta means more overdispersion.
#
# Example (from the repository root):
# Rscript --vanilla reproduction-code/simulation-code/search-nb-sparse-ci-null.R \
#   --nb-sizes=0.25,0.5,1,2 --beta-z=0,0.5,1 \
#   --x-z-slopes=0.5,1,1.5 --reps=200 --permutations=999 --cores=4 \
#   --output-dir=results/simulation-results/nb-sparse-ci-search
# For 10,000 observations and 0.5% marginal exposure, add
#   --n=10000 --x-prevalence=0.005 --permutation-engine=sparse
# The sparse engine samples exposed subsets directly, preserving the exact
# binary-label permutation law without allocating n-by-B matrices.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_arg))
script_dir <- dirname(script_file)
project_root <- normalizePath(file.path(script_dir, "..", ".."))

parse_numeric_grid <- function(value, name, positive = FALSE) {
  entries <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  numbers <- suppressWarnings(as.numeric(entries))
  if (length(numbers) == 0L || any(!nzchar(entries)) ||
      any(!is.finite(numbers)) || (positive && any(numbers <= 0)) ||
      anyDuplicated(numbers)) {
    stop("Invalid --", name, " grid: ", value, call. = FALSE)
  }
  numbers
}

parse_cli <- function(arguments) {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(detected_cores) || detected_cores < 1L) {
    detected_cores <- 1L
  }
  config <- list(
    nb_sizes = c(0.25, 0.5, 1, 2),
    beta_z = c(0, 0.5, 1),
    x_z_slopes = c(0.5, 1, 1.5),
    n = 400L,
    x_prevalence = 0.025,
    reps = 200L,
    permutations = 999L,
    permutation_engine = "dense",
    cores = min(4L, detected_cores),
    seed = 20261010L,
    alpha = 0.05,
    output_dir = file.path(
      project_root, "results", "simulation-results", "nb-sparse-ci-search"
    )
  )
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop("Arguments must have the form --name=value: ", argument,
           call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    value <- paste(pieces[-1L], collapse = "=")
    if (!key %in% names(config)) {
      stop("Unknown argument: --", pieces[[1L]], call. = FALSE)
    }
    config[[key]] <- if (key %in% c("nb_sizes", "beta_z", "x_z_slopes")) {
      parse_numeric_grid(value, pieces[[1L]], positive = key == "nb_sizes")
    } else if (key %in% c("output_dir", "permutation_engine")) {
      value
    } else if (key %in% c("alpha", "x_prevalence")) {
      suppressWarnings(as.numeric(value))
    } else {
      suppressWarnings(as.integer(value))
    }
  }
  if (length(config$n) != 1L || is.na(config$n) || config$n < 4L ||
      length(config$x_prevalence) != 1L ||
      !is.finite(config$x_prevalence) ||
      config$x_prevalence <= 0 || config$x_prevalence >= 1 ||
      length(config$reps) != 1L || is.na(config$reps) || config$reps < 1L ||
      length(config$permutations) != 1L || is.na(config$permutations) ||
      config$permutations < 1L ||
      length(config$cores) != 1L || is.na(config$cores) || config$cores < 1L ||
      length(config$seed) != 1L || is.na(config$seed) || config$seed < 1L ||
      length(config$alpha) != 1L || !is.finite(config$alpha) ||
      config$alpha <= 0 || config$alpha >= 1 ||
      !config$permutation_engine %in% c("dense", "sparse") ||
      !nzchar(config$output_dir)) {
    stop("Invalid simulation configuration.", call. = FALSE)
  }
  config$cores <- min(config$cores, config$reps)
  config$output_dir <- normalizePath(config$output_dir, mustWork = FALSE)
  config
}

config <- parse_cli(commandArgs(trailingOnly = TRUE))

output_paths <- setNames(
  file.path(config$output_dir, paste0(
    "nb-sparse-ci-search-",
    c("grid.csv", "replicates.csv", "methods.csv", "paired.csv", "all.rds")
  )),
  c("grid", "replicates", "methods", "paired", "rds")
)
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite existing search output in ",
       config$output_dir, call. = FALSE)
}

specifier_file <- file.path(
  script_dir, "sim-spec-formal", "nb-permutation-signflip.R"
)
specifier <- new.env(parent = globalenv())
source(specifier_file, local = specifier)
template <- subset(
  specifier$parameter_grid,
  null_design == "CI" & scenario_id == "rare_nb" & !ground_truth
)
if (nrow(template) != 1L || template$n != 400L ||
    abs(template$x_prevalence - 0.025) > 1e-12 ||
    template$mu0 != 10 || template$gamma_signal != 0 ||
    template$working_family != "negative_binomial" ||
    template$working_dispersion_mode != "estimated") {
  stop("The formal sparse CI null no longer matches this search.",
       call. = FALSE)
}

grid <- expand.grid(
  nb_size = config$nb_sizes,
  beta_z = config$beta_z,
  x_z_slope = config$x_z_slopes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$grid_id <- seq_len(nrow(grid))
grid$n <- config$n
grid$x_prevalence <- config$x_prevalence
grid$nb_dispersion <- 1 / grid$nb_size
grid$x_intercept <- vapply(
  grid$x_z_slope,
  function(slope) specifier$calibrate_logistic_intercept(
    config$x_prevalence, slope
  ),
  numeric(1L)
)
grid <- grid[, c(
  "grid_id", "n", "x_prevalence", "nb_size", "nb_dispersion",
  "beta_z", "x_z_slope", "x_intercept"
)]

generate_args <- formalArgs(specifier$generate_data)
information_statistic <- get(
  ".information_statistic", envir = specifier$.nb_signflip_env,
  inherits = FALSE
)
empirical_statistic <- get(
  ".studentized_statistic", envir = specifier$.nb_signflip_env,
  inherits = FALSE
)
if (identical(config$permutation_engine, "sparse")) {
  source(file.path(script_dir, "sparse-binary-permutation.R"))
}

prepare_poisson_only_cache <- function(data) {
  D <- cbind(`(Intercept)` = 1, Z = data$Z)
  # This is the first step of sceptre:::perform_response_precomputation.
  fit <- suppressWarnings(stats::glm.fit(
    y = data$Y, x = D, family = stats::poisson()
  ))
  mu_hat <- as.numeric(fit$fitted.values)
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients)) ||
      length(mu_hat) != data$n || any(!is.finite(mu_hat)) ||
      any(mu_hat <= 0)) {
    stop("Poisson null fit returned invalid estimates")
  }

  W <- mu_hat
  gram_inverse <- solve(crossprod(D, D * W))
  projection_coef <- gram_inverse %*% crossprod(D, W * data$X)
  X_tilde <- as.numeric(data$X - D %*% projection_coef)
  outcome_factor <- data$Y - mu_hat
  score_contrib <- X_tilde * outcome_factor
  information_norm_sq <- sum(X_tilde^2 * W)
  score_norm_sq <- sum(score_contrib^2)
  degenerate <- !is.finite(information_norm_sq) ||
    !is.finite(score_norm_sq) ||
    information_norm_sq <= .Machine$double.eps ||
    score_norm_sq <= .Machine$double.eps
  list(
    D = D, W = W, gram_inverse = gram_inverse,
    outcome_factor = outcome_factor, X_tilde = X_tilde,
    score_contrib = score_contrib, degenerate = degenerate
  )
}

set.seed(config$seed)
seed_matrix <- matrix(
  sample.int(
    .Machine$integer.max,
    size = nrow(grid) * config$reps,
    replace = FALSE
  ),
  nrow = nrow(grid), ncol = config$reps
)

# Compute both reference distributions using exactly the same permutations.
# This reproduces the two-sided, plus-one tests in the formal specifier while
# avoiding a second projection and permutation pass for our method.
paired_p_values <- function(data) {
  cache <- data$score_cache
  observed_information <- information_statistic(
    cache$score_contrib, cache$X_tilde, cache$W
  )[1L]
  observed_empirical <- empirical_statistic(cache$score_contrib)[1L]
  hits_information <- 0L
  hits_empirical <- 0L
  completed <- 0L
  set.seed(data$perm_seed)
  while (completed < config$permutations) {
    batch <- min(128L, config$permutations - completed)
    permutation_indices <- vapply(
      seq_len(batch), function(i) sample.int(data$n), integer(data$n)
    )
    X_permuted <- matrix(
      data$X[permutation_indices], nrow = data$n, ncol = batch
    )
    projection_coef <- cache$gram_inverse %*%
      crossprod(cache$D, X_permuted * cache$W)
    X_tilde_permuted <- X_permuted - cache$D %*% projection_coef
    score_permuted <- X_tilde_permuted * cache$outcome_factor
    reference_information <- information_statistic(
      score_permuted, X_tilde_permuted, cache$W
    )
    reference_empirical <- empirical_statistic(score_permuted)
    hits_information <- hits_information +
      sum(abs(reference_information) >= abs(observed_information))
    hits_empirical <- hits_empirical +
      sum(abs(reference_empirical) >= abs(observed_empirical))
    completed <- completed + batch
  }
  list(
    information_statistic = observed_information,
    empirical_statistic = observed_empirical,
    information_hits = hits_information,
    empirical_hits = hits_empirical,
    p_barry = (1 + hits_information) / (config$permutations + 1),
    p_ours = (1 + hits_empirical) / (config$permutations + 1)
  )
}

# Same raw-label permutations and information statistic as Barry, now with
# H = 1 and W = mu_hat.  Resetting perm_seed gives the exact same draws as
# paired_p_values even when the NB fit fails.
poisson_only_p_value <- function(data, cache) {
  if (isTRUE(cache$degenerate)) {
    return(list(
      poisson_only_statistic = 0,
      poisson_only_hits = config$permutations,
      p_poisson_only = 1
    ))
  }
  observed <- information_statistic(
    cache$score_contrib, cache$X_tilde, cache$W
  )[1L]
  hits <- 0L
  completed <- 0L
  set.seed(data$perm_seed)
  while (completed < config$permutations) {
    batch <- min(128L, config$permutations - completed)
    permutation_indices <- vapply(
      seq_len(batch), function(i) sample.int(data$n), integer(data$n)
    )
    X_permuted <- matrix(
      data$X[permutation_indices], nrow = data$n, ncol = batch
    )
    projection_coef <- cache$gram_inverse %*%
      crossprod(cache$D, X_permuted * cache$W)
    X_tilde_permuted <- X_permuted - cache$D %*% projection_coef
    score_permuted <- X_tilde_permuted * cache$outcome_factor
    reference <- information_statistic(
      score_permuted, X_tilde_permuted, cache$W
    )
    hits <- hits + sum(abs(reference) >= abs(observed))
    completed <- completed + batch
  }
  list(
    poisson_only_statistic = observed,
    poisson_only_hits = hits,
    p_poisson_only = (1 + hits) / (config$permutations + 1)
  )
}

run_one <- function(grid_index, replicate_index) {
  parameters <- grid[grid_index, , drop = FALSE]
  data_seed <- seed_matrix[grid_index, replicate_index]
  result <- data.frame(
    grid_id = grid_index,
    replicate = replicate_index,
    data_seed = data_seed,
    permutation_seed = NA_integer_,
    success = FALSE,
    failure_reason = NA_character_,
    success_ours_normal = FALSE,
    failure_reason_ours_normal = NA_character_,
    success_poisson_only = FALSE,
    failure_reason_poisson_only = NA_character_,
    degenerate_poisson_only = FALSE,
    degenerate = FALSE,
    theta_hat = NA_real_,
    n_exposed = NA_integer_,
    n_exposed_nonzero = NA_integer_,
    score_eff_n = NA_real_,
    max_score_share = NA_real_,
    information_statistic = NA_real_,
    empirical_statistic = NA_real_,
    information_hits = NA_integer_,
    empirical_hits = NA_integer_,
    p_barry = NA_real_,
    p_ours = NA_real_,
    p_ours_normal = NA_real_,
    poisson_only_statistic = NA_real_,
    poisson_only_hits = NA_integer_,
    p_poisson_only = NA_real_,
    reject_barry = FALSE,
    reject_ours = FALSE,
    reject_ours_normal = FALSE,
    reject_poisson_only = FALSE,
    stringsAsFactors = FALSE
  )

  arguments <- as.list(template[1L, generate_args, drop = FALSE])
  arguments$n <- config$n
  arguments$x_prevalence <- config$x_prevalence
  arguments$nb_size <- parameters$nb_size
  arguments$nb_size_intercept <- NA_real_
  arguments$nb_size_z_slope <- NA_real_
  arguments$beta_z <- parameters$beta_z
  arguments$x_z_slope <- parameters$x_z_slope
  arguments$x_intercept <- parameters$x_intercept
  arguments$B_resample <- config$permutations
  set.seed(data_seed)
  data <- tryCatch(
    do.call(specifier$generate_data, arguments),
    error = function(e) e
  )
  if (inherits(data, "error")) {
    result$failure_reason <- paste0(
      "data generation failed: ", conditionMessage(data)
    )
    result$failure_reason_ours_normal <- result$failure_reason
    result$failure_reason_poisson_only <- result$failure_reason
    return(result)
  }
  result$permutation_seed <- data$perm_seed
  result$theta_hat <- data$theta_hat
  result$n_exposed <- data$n_exposed
  result$n_exposed_nonzero <- data$n_exposed_nonzero

  poisson_cache <- tryCatch(
    prepare_poisson_only_cache(data), error = function(e) e
  )
  if (identical(config$permutation_engine, "sparse")) {
    if (inherits(poisson_cache, "error")) {
      result$failure_reason_poisson_only <- paste0(
        "Poisson score preparation failed: ", conditionMessage(poisson_cache)
      )
      poisson_cache <- NULL
    }
    if (!isTRUE(data$success)) {
      result$failure_reason <- data$failure_reason
      result$failure_reason_ours_normal <- data$failure_reason
      data$score_cache <- NULL
    } else {
      result$score_eff_n <- data$score_cache$score_eff_n
      result$max_score_share <- data$score_cache$max_score_share
      result$degenerate <- isTRUE(data$score_cache$degenerate)
      observed_empirical <- if (result$degenerate) {
        0
      } else {
        tryCatch(
          empirical_statistic(data$score_cache$score_contrib)[1L],
          error = function(e) e
        )
      }
      if (inherits(observed_empirical, "error")) {
        result$failure_reason_ours_normal <- paste0(
          "empirical statistic failed: ", conditionMessage(observed_empirical)
        )
      } else {
        result$empirical_statistic <- observed_empirical
        normal_p <- 2 * stats::pnorm(-abs(observed_empirical))
        if (length(normal_p) == 1L && is.finite(normal_p) &&
            normal_p >= 0 && normal_p <= 1) {
          result$p_ours_normal <- normal_p
          result$success_ours_normal <- TRUE
          result$reject_ours_normal <- normal_p <= config$alpha
        } else {
          result$failure_reason_ours_normal <-
            "normal reference returned an invalid p-value"
        }
      }
    }
    if (is.null(data$score_cache) && is.null(poisson_cache)) {
      return(result)
    }
    sparse <- tryCatch(
      sparse_binary_permutation_p_values(
        data, config$permutations, poisson_cache
      ),
      error = function(e) e
    )
    if (inherits(sparse, "error")) {
      failure <- paste0(
        "sparse permutation calibration failed: ", conditionMessage(sparse)
      )
      if (isTRUE(data$success)) result$failure_reason <- failure
      if (!is.null(poisson_cache)) result$failure_reason_poisson_only <- failure
      return(result)
    }
    for (name in names(sparse)) result[[name]] <- sparse[[name]]
    if (isTRUE(data$success)) {
      valid_nb <- all(is.finite(c(result$p_barry, result$p_ours))) &&
        all(c(result$p_barry, result$p_ours) >= 0) &&
        all(c(result$p_barry, result$p_ours) <= 1)
      if (valid_nb) {
        result$success <- TRUE
        result$reject_barry <- result$p_barry <= config$alpha
        result$reject_ours <- result$p_ours <= config$alpha
      } else {
        result$failure_reason <- "sparse calibration returned invalid NB p-values"
        result$p_barry <- NA_real_
        result$p_ours <- NA_real_
      }
    }
    if (!is.null(poisson_cache)) {
      valid_poisson <- is.finite(result$p_poisson_only) &&
        result$p_poisson_only >= 0 && result$p_poisson_only <= 1
      if (valid_poisson) {
        result$success_poisson_only <- TRUE
        result$degenerate_poisson_only <- isTRUE(poisson_cache$degenerate)
        result$reject_poisson_only <- result$p_poisson_only <= config$alpha
      } else {
        result$failure_reason_poisson_only <-
          "sparse calibration returned an invalid Poisson p-value"
        result$p_poisson_only <- NA_real_
      }
    }
    return(result)
  }
  if (inherits(poisson_cache, "error")) {
    result$failure_reason_poisson_only <- paste0(
      "Poisson score preparation failed: ", conditionMessage(poisson_cache)
    )
  } else {
    poisson_result <- tryCatch(
      poisson_only_p_value(data, poisson_cache), error = function(e) e
    )
    if (inherits(poisson_result, "error")) {
      result$failure_reason_poisson_only <- paste0(
        "Poisson permutation calibration failed: ",
        conditionMessage(poisson_result)
      )
    } else {
      for (name in names(poisson_result)) {
        result[[name]] <- poisson_result[[name]]
      }
      if (is.finite(result$p_poisson_only) &&
          result$p_poisson_only >= 0 && result$p_poisson_only <= 1) {
        result$success_poisson_only <- TRUE
        result$degenerate_poisson_only <- poisson_cache$degenerate
        result$reject_poisson_only <- result$p_poisson_only <= config$alpha
      } else {
        result$failure_reason_poisson_only <-
          "Poisson calibration returned an invalid p-value"
        result$p_poisson_only <- NA_real_
      }
    }
  }

  if (!isTRUE(data$success)) {
    result$failure_reason <- data$failure_reason
    result$failure_reason_ours_normal <- data$failure_reason
    return(result)
  }

  result$score_eff_n <- data$score_cache$score_eff_n
  result$max_score_share <- data$score_cache$max_score_share
  if (isTRUE(data$score_cache$degenerate)) {
    result$success <- TRUE
    result$degenerate <- TRUE
    result$information_statistic <- 0
    result$empirical_statistic <- 0
    result$information_hits <- config$permutations
    result$empirical_hits <- config$permutations
    result$p_barry <- 1
    result$p_ours <- 1
    result$p_ours_normal <- 1
    result$success_ours_normal <- TRUE
    return(result)
  }

  # This uses exactly our observed statistic, with only its reference law
  # changed from paired permutations to the standard normal distribution.
  observed_empirical <- tryCatch(
    empirical_statistic(data$score_cache$score_contrib)[1L],
    error = function(e) e
  )
  if (inherits(observed_empirical, "error")) {
    result$failure_reason_ours_normal <- paste0(
      "empirical statistic failed: ", conditionMessage(observed_empirical)
    )
  } else {
    result$empirical_statistic <- observed_empirical
    normal_p <- 2 * stats::pnorm(-abs(observed_empirical))
    if (length(normal_p) == 1L && is.finite(normal_p) &&
        normal_p >= 0 && normal_p <= 1) {
      result$p_ours_normal <- normal_p
      result$success_ours_normal <- TRUE
      result$reject_ours_normal <- normal_p <= config$alpha
    } else {
      result$failure_reason_ours_normal <-
        "normal reference returned an invalid p-value"
    }
  }

  paired <- tryCatch(paired_p_values(data), error = function(e) e)
  if (inherits(paired, "error")) {
    result$failure_reason <- paste0(
      "permutation calibration failed: ", conditionMessage(paired)
    )
    return(result)
  }
  for (name in names(paired)) {
    result[[name]] <- paired[[name]]
  }
  if (any(!is.finite(c(result$p_barry, result$p_ours))) ||
      any(c(result$p_barry, result$p_ours) < 0) ||
      any(c(result$p_barry, result$p_ours) > 1)) {
    result$failure_reason <- "calibration returned invalid p-values"
    result$p_barry <- NA_real_
    result$p_ours <- NA_real_
    return(result)
  }
  result$success <- TRUE
  result$reject_barry <- result$p_barry <= config$alpha
  result$reject_ours <- result$p_ours <= config$alpha
  result
}

binomial_interval <- function(rejected, total) {
  if (total == 0L) {
    return(c(NA_real_, NA_real_))
  }
  stats::binom.test(rejected, total)$conf.int[1:2]
}

method_summary <- function(rows, method) {
  rejected <- rows[[paste0("reject_", method)]]
  valid <- if (identical(method, "ours_normal")) {
    rows$success_ours_normal
  } else if (identical(method, "poisson_only")) {
    rows$success_poisson_only
  } else {
    rows$success
  }
  all_interval <- binomial_interval(sum(rejected), nrow(rows))
  valid_interval <- binomial_interval(sum(rejected[valid]), sum(valid))
  data.frame(
    method = method,
    n_all = nrow(rows),
    n_valid = sum(valid),
    n_failed = sum(!valid),
    n_rejected = sum(rejected),
    rate_all = mean(rejected),
    ci_all_lower = all_interval[1L],
    ci_all_upper = all_interval[2L],
    rate_valid = if (any(valid)) mean(rejected[valid]) else NA_real_,
    ci_valid_lower = valid_interval[1L],
    ci_valid_upper = valid_interval[2L],
    stringsAsFactors = FALSE
  )
}

paired_summary <- function(rows, valid_only) {
  if (valid_only) {
    rows <- rows[rows$success, , drop = FALSE]
  }
  n <- nrow(rows)
  if (n == 0L) {
    return(data.frame(
      subset = "valid", n = 0L, both = 0L, barry_only = 0L,
      ours_only = 0L, neither = 0L, rate_barry = NA_real_,
      rate_ours = NA_real_, ours_minus_barry = NA_real_,
      difference_ci_lower = NA_real_, difference_ci_upper = NA_real_,
      exact_discordance_p = NA_real_
    ))
  }
  barry <- rows$reject_barry
  ours <- rows$reject_ours
  difference <- as.numeric(ours) - as.numeric(barry)
  barry_only <- sum(barry & !ours)
  ours_only <- sum(ours & !barry)
  difference_se <- if (n > 1L) stats::sd(difference) / sqrt(n) else NA_real_
  difference_ci <- mean(difference) + c(-1, 1) *
    stats::qnorm(0.975) * difference_se
  difference_ci <- pmin(1, pmax(-1, difference_ci))
  discordance_p <- if (barry_only + ours_only == 0L) {
    1
  } else {
    stats::binom.test(
      ours_only, barry_only + ours_only, p = 0.5
    )$p.value
  }
  data.frame(
    subset = if (valid_only) "valid" else "all",
    n = n,
    both = sum(barry & ours),
    barry_only = barry_only,
    ours_only = ours_only,
    neither = sum(!barry & !ours),
    rate_barry = mean(barry),
    rate_ours = mean(ours),
    ours_minus_barry = mean(difference),
    difference_ci_lower = difference_ci[1L],
    difference_ci_upper = difference_ci[2L],
    exact_discordance_p = discordance_p,
    stringsAsFactors = FALSE
  )
}

cat(sprintf(
  "Sparse CI null grid: %d cells x %d replicates, n=%d, Pr(X=1)=%.4g, B=%d, engine=%s, cores=%d, seed=%d\n",
  nrow(grid), config$reps, config$n, config$x_prevalence,
  config$permutations, config$permutation_engine, config$cores, config$seed
))
started <- Sys.time()
all_results <- vector("list", nrow(grid))
for (grid_index in seq_len(nrow(grid))) {
  cell_started <- Sys.time()
  rows <- if (.Platform$OS.type == "unix" && config$cores > 1L) {
    parallel::mclapply(
      seq_len(config$reps),
      function(replicate_index) run_one(grid_index, replicate_index),
      mc.cores = config$cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(
      seq_len(config$reps),
      function(replicate_index) run_one(grid_index, replicate_index)
    )
  }
  if (any(vapply(rows, inherits, logical(1L), what = "try-error"))) {
    stop("A worker failed in grid cell ", grid_index, call. = FALSE)
  }
  all_results[[grid_index]] <- do.call(rbind, rows)
  cell <- grid[grid_index, ]
  cell_rows <- all_results[[grid_index]]
  cat(sprintf(
    "  %d/%d theta=%.4g beta_z=%.4g x_z_slope=%.4g: Barry %.3f, ours %.3f, ours-normal %.3f, Poisson-only %.3f, NB failures %d, Poisson failures %d (%.1fs)\n",
    grid_index, nrow(grid), cell$nb_size, cell$beta_z, cell$x_z_slope,
    mean(cell_rows$reject_barry), mean(cell_rows$reject_ours),
    mean(cell_rows$reject_ours_normal),
    mean(cell_rows$reject_poisson_only), sum(!cell_rows$success),
    sum(!cell_rows$success_poisson_only),
    as.numeric(difftime(Sys.time(), cell_started, units = "secs"))
  ))
  flush.console()
}

replicates <- do.call(rbind, all_results)
rownames(replicates) <- NULL
methods <- do.call(rbind, lapply(seq_len(nrow(grid)), function(index) {
  rows <- replicates[replicates$grid_id == index, , drop = FALSE]
  cbind(
    grid[rep(index, 4L), , drop = FALSE],
    rbind(
      method_summary(rows, "barry"),
      method_summary(rows, "ours"),
      method_summary(rows, "ours_normal"),
      method_summary(rows, "poisson_only")
    )
  )
}))
paired <- do.call(rbind, lapply(seq_len(nrow(grid)), function(index) {
  rows <- replicates[replicates$grid_id == index, , drop = FALSE]
  cbind(
    grid[rep(index, 2L), , drop = FALSE],
    rbind(paired_summary(rows, FALSE), paired_summary(rows, TRUE))
  )
}))
rownames(methods) <- NULL
rownames(paired) <- NULL

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
if (any(file.exists(output_paths))) {
  stop("Output files appeared while the search was running; refusing overwrite.",
       call. = FALSE)
}
utils::write.csv(grid, output_paths[["grid"]], row.names = FALSE)
utils::write.csv(replicates, output_paths[["replicates"]], row.names = FALSE,
                 na = "")
utils::write.csv(methods, output_paths[["methods"]], row.names = FALSE,
                 na = "")
utils::write.csv(paired, output_paths[["paired"]], row.names = FALSE,
                 na = "")
saveRDS(list(
  configuration = c(config, list(
    mu0 = 10,
    gamma_signal = 0,
    null_design = "CI",
    test = paste(
      "two-sided plus-one raw-label permutation;",
      "ours_normal uses a two-sided standard-normal reference"
    ),
    specifier = normalizePath(specifier_file),
    started = format(started, tz = "UTC", usetz = TRUE),
    completed = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )),
  grid = grid,
  replicates = replicates,
  methods = methods,
  paired = paired
), output_paths[["rds"]])
cat("Saved search output to ", config$output_dir, "\n", sep = "")
