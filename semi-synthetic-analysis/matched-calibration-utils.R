# Utilities shared by the CCLE matched-calibration runner and its smoke test.
#
# The normal and permutation p-values for each method must be calibrations of
# one observed statistic.  Accordingly, each pcmt permutation test is called
# exactly once below; its returned statistic supplies the normal calibration,
# and its returned p.value supplies the permutation calibration.

CCLE_SEED_COMPONENT_OFFSETS <- c(
  data_seed = 1L,
  pcm_split_seed = 2L,
  gcm_permutation_seed = 3L,
  pcm_permutation_seed = 4L,
  gcm_YonZ_seed = 5L,
  gcm_XonZ_seed = 6L,
  pcm_YonXZ_seed = 7L,
  pcm_YonZ_seed = 8L,
  pcm_YhatonZ_seed = 9L,
  pcm_RonZ_seed = 10L
)

ccle_seed_schedule <- function(setup_idx, rep_idx, seed_base = 20260000L) {
  setup_idx <- as.integer(setup_idx)
  rep_idx <- as.integer(rep_idx)
  seed_base <- as.integer(seed_base)

  if (length(setup_idx) != 1L || is.na(setup_idx) || setup_idx < 1L) {
    stop("setup_idx must be one positive integer", call. = FALSE)
  }
  if (length(rep_idx) != 1L || is.na(rep_idx) || rep_idx < 1L) {
    stop("rep_idx must be one positive integer", call. = FALSE)
  }
  if (length(seed_base) != 1L || is.na(seed_base) || seed_base < 1L) {
    stop("seed_base must be one positive integer", call. = FALSE)
  }

  # A 100-seed block is reserved for every replicate, and a 100,000-seed
  # block for every setup.  The current design (60 setups x 500 replicates)
  # therefore has disjoint streams with ample room for additional components.
  seeds <- as.double(seed_base) +
    as.double(setup_idx) * 100000 +
    as.double(rep_idx) * 100 +
    as.double(CCLE_SEED_COMPONENT_OFFSETS)

  if (any(seeds > .Machine$integer.max)) {
    stop("Derived seed exceeds R's integer seed range", call. = FALSE)
  }

  stats::setNames(as.integer(seeds), names(CCLE_SEED_COMPONENT_OFFSETS))
}

ccle_with_seed <- function(seed, code) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed) || seed < 1L) {
    stop("seed must be one positive integer", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

ccle_make_pcm_split <- function(n, split_seed, frac = 0.5) {
  n <- as.integer(n)
  frac <- as.numeric(frac)
  if (length(n) != 1L || is.na(n) || n < 4L) {
    stop("n must be one integer of at least 4", call. = FALSE)
  }
  if (length(frac) != 1L || !is.finite(frac) || frac <= 0 || frac >= 1) {
    stop("frac must lie strictly between 0 and 1", call. = FALSE)
  }

  n_D2 <- ceiling(frac * n)
  if (n_D2 < 1L || n_D2 >= n) {
    stop("frac produces an empty PCM split", call. = FALSE)
  }

  ccle_with_seed(split_seed, sort(sample.int(n, size = n_D2)))
}

ccle_normal_pvalue <- function(statistic,
                               alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  statistic <- as.numeric(statistic)
  if (length(statistic) != 1L) {
    stop("statistic must have length one", call. = FALSE)
  }

  pvalue <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
    greater = stats::pnorm(statistic, lower.tail = FALSE),
    less = stats::pnorm(statistic)
  )
  pmin(1, pvalue)
}

ccle_adaptive_calibration <- function(
    statistic,
    perm_stats,
    h,
    alpha,
    alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  statistic <- as.numeric(statistic)
  perm_stats <- as.numeric(perm_stats)
  h <- as.integer(h)
  alpha <- as.numeric(alpha)

  if (length(statistic) != 1L || !is.finite(statistic)) {
    stop("statistic must be one finite number", call. = FALSE)
  }
  if (!length(perm_stats) || any(!is.finite(perm_stats))) {
    stop("perm_stats must be a nonempty finite numeric vector", call. = FALSE)
  }
  if (length(h) != 1L || is.na(h) || h < 1L) {
    stop("h must be one positive integer", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between 0 and 1", call. = FALSE)
  }

  score <- function(value) {
    switch(
      alternative,
      two.sided = abs(value),
      greater = value,
      less = -value
    )
  }

  S0 <- score(statistic)
  S <- score(perm_stats)
  K <- 0L
  t_star <- Inf
  actual_B <- length(perm_stats)
  stop_reason <- "cap"
  trajectory_K <- integer(length(perm_stats))
  trajectory_p <- numeric(length(perm_stats))

  for (b in seq_along(perm_stats)) {
    if (S[[b]] >= S0) {
      K <- K + 1L
      if (K == h && is.infinite(t_star)) {
        t_star <- b
      }
    }

    pvalue <- if (K < h) {
      h / (b + h - K)
    } else {
      h / t_star
    }
    trajectory_K[[b]] <- K
    trajectory_p[[b]] <- pvalue

    if (pvalue <= alpha) {
      actual_B <- b
      stop_reason <- "rejection"
      break
    }
    if (K == h && pvalue > alpha) {
      actual_B <- b
      stop_reason <- "futility"
      break
    }
  }

  list(
    p.value = trajectory_p[[actual_B]],
    perm_stats = perm_stats[seq_len(actual_B)],
    adaptive = list(
      h = h,
      alpha = alpha,
      K = K,
      t_star = t_star,
      stopped_early = actual_B < length(perm_stats),
      stop_reason = stop_reason,
      trajectory = data.frame(
        t = seq_len(actual_B),
        S = S[seq_len(actual_B)],
        K = trajectory_K[seq_len(actual_B)],
        p_AV = trajectory_p[seq_len(actual_B)]
      )
    )
  )
}

ccle_rf_args <- function(seed, num_trees = 500L) {
  num_trees <- as.integer(num_trees)
  if (length(num_trees) != 1L || is.na(num_trees) || num_trees < 1L) {
    stop("num_trees must be one positive integer", call. = FALSE)
  }
  list(
    seed = as.integer(seed),
    num.threads = 1L,
    num.trees = num_trees
  )
}

ccle_empty_diagnostic <- function() {
  list(
    actual_B = NA_integer_,
    stop_reason = "not_run",
    error = NA_character_,
    matched_statistic = FALSE,
    normal_statistic = NA_real_,
    permutation_statistic = NA_real_,
    loss_count = NA_integer_,
    t_star = NA_real_,
    stopped_early = NA,
    adaptive_backend = NA_character_
  )
}

ccle_run_matched_calibration <- function(
    Y,
    X,
    Z,
    seeds,
    pcm_D2_indices,
    B_max = 499L,
    adaptive_h = 15L,
    alpha = 0.05,
    n_cores = 1L,
    alternative = c("two.sided", "greater", "less"),
    pcm_frac = 0.5,
    rf_num_trees = 500L,
    test_order = c("gcm", "pcm")) {
  alternative <- match.arg(alternative)
  B_max <- as.integer(B_max)
  adaptive_h <- as.integer(adaptive_h)
  alpha <- as.numeric(alpha)
  n_cores <- as.integer(n_cores)
  pcm_frac <- as.numeric(pcm_frac)

  if (length(B_max) != 1L || is.na(B_max) || B_max < 1L) {
    stop("B_max must be one positive integer", call. = FALSE)
  }
  if (length(adaptive_h) != 1L || is.na(adaptive_h) || adaptive_h < 1L) {
    stop("adaptive_h must be one positive integer", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between 0 and 1", call. = FALSE)
  }
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores != 1L) {
    stop(
      "Matched CCLE runs require n_cores = 1 for deterministic RF refits",
      call. = FALSE
    )
  }
  if (!setequal(test_order, c("gcm", "pcm")) || length(test_order) != 2L) {
    stop("test_order must contain gcm and pcm exactly once", call. = FALSE)
  }

  required_seeds <- names(CCLE_SEED_COMPONENT_OFFSETS)
  if (!all(required_seeds %in% names(seeds))) {
    stop(
      "Missing component seeds: ",
      paste(setdiff(required_seeds, names(seeds)), collapse = ", "),
      call. = FALSE
    )
  }
  seeds <- as.integer(seeds[required_seeds])
  names(seeds) <- required_seeds

  pcm_D2_indices <- as.integer(pcm_D2_indices)
  if (!length(pcm_D2_indices) || anyNA(pcm_D2_indices)) {
    stop("pcm_D2_indices must be an explicit nonempty integer vector", call. = FALSE)
  }

  pvalues <- c(
    asym_gcm = NA_real_,
    asym_pcm = NA_real_,
    perm_gcm = NA_real_,
    perm_pcm = NA_real_
  )
  statistics <- c(gcm = NA_real_, pcm = NA_real_)
  diagnostics <- list(
    gcm = ccle_empty_diagnostic(),
    pcm = ccle_empty_diagnostic()
  )

  rf_args <- list(
    gcm_YonZ = ccle_rf_args(seeds[["gcm_YonZ_seed"]], rf_num_trees),
    # The same role seed is reused by the observed and every permuted X|Z fit.
    gcm_XonZ = ccle_rf_args(seeds[["gcm_XonZ_seed"]], rf_num_trees),
    pcm_YonXZ = ccle_rf_args(seeds[["pcm_YonXZ_seed"]], rf_num_trees),
    pcm_YonZ = ccle_rf_args(seeds[["pcm_YonZ_seed"]], rf_num_trees),
    pcm_YhatonZ = ccle_rf_args(seeds[["pcm_YhatonZ_seed"]], rf_num_trees),
    # The same role seed is reused by the observed and every permuted R|Z fit.
    pcm_RonZ = ccle_rf_args(seeds[["pcm_RonZ_seed"]], rf_num_trees)
  )

  apply_adaptive_fallback <- function(fit) {
    calibrated <- ccle_adaptive_calibration(
      statistic = fit$statistic,
      perm_stats = fit$perm_stats,
      h = adaptive_h,
      alpha = alpha,
      alternative = alternative
    )
    fit$p.value <- calibrated$p.value
    fit$perm_stats <- calibrated$perm_stats
    fit$parameter[["B"]] <- length(calibrated$perm_stats)
    fit$adaptive <- calibrated$adaptive
    fit
  }

  adaptive_call_args <- function(fn) {
    fn_formals <- names(formals(fn))
    if (all(c("h", "alpha") %in% fn_formals)) {
      list(h = adaptive_h, alpha = alpha, backend = "pcmt_native")
    } else {
      list(backend = "ordered_perm_stats_fallback")
    }
  }

  for (test_name in test_order) {
    if (identical(test_name, "gcm")) {
      adaptive_args <- adaptive_call_args(pcmt::perm_gcm)
      attempt <- tryCatch(
        {
          fit <- ccle_with_seed(
            seeds[["gcm_permutation_seed"]],
            do.call(
              pcmt::perm_gcm,
              c(
                list(
                  Y = Y,
                  X = X,
                  Z = Z,
                  reg_YonZ = "rf",
                  reg_XonZ = "rf",
                  args_YonZ = rf_args$gcm_YonZ,
                  args_XonZ = rf_args$gcm_XonZ,
                  B = B_max,
                  n_cores = n_cores,
                  alternative = alternative
                ),
                adaptive_args[intersect(names(adaptive_args), c("h", "alpha"))]
              )
            )
          )
          if (identical(adaptive_args$backend, "ordered_perm_stats_fallback")) {
            fit <- apply_adaptive_fallback(fit)
          }
          if (!identical(fit$alternative, alternative)) {
            stop("perm_gcm returned a different alternative")
          }
          list(fit = fit, error = NULL)
        },
        error = function(e) list(fit = NULL, error = conditionMessage(e))
      )

      if (is.null(attempt$error)) {
        fit <- attempt$fit
        T_obs <- as.numeric(fit$statistic)
        actual_B <- length(fit$perm_stats)
        statistics[["gcm"]] <- T_obs
        pvalues[["asym_gcm"]] <- ccle_normal_pvalue(T_obs, alternative)
        pvalues[["perm_gcm"]] <- as.numeric(fit$p.value)
        diagnostics$gcm <- list(
          actual_B = as.integer(actual_B),
          stop_reason = fit$adaptive$stop_reason,
          error = NA_character_,
          matched_statistic = TRUE,
          normal_statistic = T_obs,
          permutation_statistic = T_obs,
          loss_count = as.integer(fit$adaptive$K),
          t_star = as.numeric(fit$adaptive$t_star),
          stopped_early = isTRUE(fit$adaptive$stopped_early),
          adaptive_backend = adaptive_args$backend
        )
      } else {
        diagnostics$gcm$stop_reason <- "error"
        diagnostics$gcm$error <- attempt$error
      }
    } else {
      adaptive_args <- adaptive_call_args(pcmt::perm_pcm)
      attempt <- tryCatch(
        {
          fit <- ccle_with_seed(
            seeds[["pcm_permutation_seed"]],
            do.call(
              pcmt::perm_pcm,
              c(
                list(
                  Y = Y,
                  X = X,
                  Z = Z,
                  reg_YonXZ = "rf",
                  reg_YonZ = "rf",
                  reg_YhatonZ = "rf",
                  reg_RonZ = "rf",
                  args_YonXZ = rf_args$pcm_YonXZ,
                  args_YonZ = rf_args$pcm_YonZ,
                  args_YhatonZ = rf_args$pcm_YhatonZ,
                  args_RonZ = rf_args$pcm_RonZ,
                  frac = pcm_frac,
                  indices = pcm_D2_indices,
                  B = B_max,
                  n_cores = n_cores,
                  alternative = alternative
                ),
                adaptive_args[intersect(names(adaptive_args), c("h", "alpha"))]
              )
            )
          )
          if (identical(adaptive_args$backend, "ordered_perm_stats_fallback")) {
            fit <- apply_adaptive_fallback(fit)
          }
          if (!identical(fit$alternative, alternative)) {
            stop("perm_pcm returned a different alternative")
          }
          if (!identical(as.integer(fit$D2_indices), pcm_D2_indices)) {
            stop("perm_pcm did not use the supplied D2 indices")
          }
          list(fit = fit, error = NULL)
        },
        error = function(e) list(fit = NULL, error = conditionMessage(e))
      )

      if (is.null(attempt$error)) {
        fit <- attempt$fit
        T_obs <- as.numeric(fit$statistic)
        actual_B <- length(fit$perm_stats)
        statistics[["pcm"]] <- T_obs
        pvalues[["asym_pcm"]] <- ccle_normal_pvalue(T_obs, alternative)
        pvalues[["perm_pcm"]] <- as.numeric(fit$p.value)
        diagnostics$pcm <- list(
          actual_B = as.integer(actual_B),
          stop_reason = fit$adaptive$stop_reason,
          error = NA_character_,
          matched_statistic = TRUE,
          normal_statistic = T_obs,
          permutation_statistic = T_obs,
          loss_count = as.integer(fit$adaptive$K),
          t_star = as.numeric(fit$adaptive$t_star),
          stopped_early = isTRUE(fit$adaptive$stopped_early),
          adaptive_backend = adaptive_args$backend
        )
      } else {
        diagnostics$pcm$stop_reason <- "error"
        diagnostics$pcm$error <- attempt$error
      }
    }
  }

  list(
    pvalues = pvalues,
    statistics = statistics,
    diagnostics = diagnostics,
    matched_calibration = isTRUE(diagnostics$gcm$matched_statistic) &&
      isTRUE(diagnostics$pcm$matched_statistic),
    alternative = alternative
  )
}
