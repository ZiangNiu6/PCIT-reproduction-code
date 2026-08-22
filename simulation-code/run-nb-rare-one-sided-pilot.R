#!/usr/bin/env Rscript

###############################################################################
# Paired one-sided pilot for the rare-exposure conditional-independence design
# in Section 5.1 of combine.pdf. The default uses the manuscript's negative
# signals and lower tail; --direction=up mirrors the signals and uses the
# upper tail.
#
# This deliberately reports nominal-cutoff rejection rates. It does not use an
# empirical-null cutoff or any other size adjustment. Barry et al.'s
# information studentization and empirical studentization share each dataset,
# fitted null, and raw-label permutation stream.
###############################################################################

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this script with Rscript.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", file_arg))
script_dir <- dirname(script_file)
project_root <- normalizePath(file.path(script_dir, "..", ".."))

parse_cli <- function(arguments) {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(detected_cores) || detected_cores < 1L) {
    detected_cores <- 1L
  }
  config <- list(
    reps = 500L,
    permutations = 1999L,
    cores = min(4L, detected_cores),
    seed = 20260822L,
    alpha = 0.05,
    direction = "down",
    output_dir = file.path(
      project_root, "results", "simulation-results",
      "nb-rare-one-sided-pilot"
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
    config[[key]] <- if (key %in% c("output_dir", "direction")) {
      value
    } else if (key == "alpha") {
      as.numeric(value)
    } else {
      as.integer(value)
    }
  }
  if (!is.finite(config$reps) || config$reps < 1L ||
      !is.finite(config$permutations) || config$permutations < 1L ||
      !is.finite(config$cores) || config$cores < 1L ||
      !is.finite(config$seed) || config$seed < 1L ||
      !is.finite(config$alpha) || config$alpha <= 0 || config$alpha >= 1 ||
      !config$direction %in% c("down", "up")) {
    stop("Invalid simulation configuration.", call. = FALSE)
  }
  config$cores <- min(config$cores, config$reps)
  config$output_dir <- normalizePath(
    config$output_dir, mustWork = FALSE
  )
  config
}

config <- parse_cli(commandArgs(trailingOnly = TRUE))

specifier_file <- file.path(
  script_dir, "sim-spec-formal", "nb-permutation-signflip.R"
)
if (!file.exists(specifier_file)) {
  stop("Missing formal simulation specifier: ", specifier_file,
       call. = FALSE)
}
specifier <- new.env(parent = globalenv())
source(specifier_file, local = specifier)

pilot_grid <- subset(
  specifier$parameter_grid,
  null_design == "CI" & scenario_id == "rare_nb"
)
pilot_grid <- pilot_grid[order(pilot_grid$gamma_signal, decreasing = TRUE), ]
expected_negative_gamma <- c(0, -0.2125, -0.425, -0.6375, -0.85)
if (nrow(pilot_grid) != 5L ||
    !isTRUE(all.equal(pilot_grid$gamma_signal, expected_negative_gamma,
                      tolerance = 1e-12)) ||
    any(pilot_grid$n != 400L)) {
  stop("The formal CI rare-NB grid does not match Section 5.1.",
       call. = FALSE)
}
pilot_grid$template_gamma_signal <- pilot_grid$gamma_signal
if (config$direction == "up") {
  pilot_grid$gamma_signal <- -expected_negative_gamma
  pilot_grid$fold_change <- exp(pilot_grid$gamma_signal)
}
expected_gamma <- pilot_grid$gamma_signal
primary_tail <- config$direction
primary_barry_column <- if (config$direction == "down") {
  "p_down_barry"
} else {
  "p_up_barry"
}
primary_empirical_column <- if (config$direction == "down") {
  "p_down_empirical"
} else {
  "p_up_empirical"
}

generate_args <- formalArgs(specifier$generate_data)
information_statistic <- get(
  ".information_statistic", envir = specifier$.nb_signflip_env,
  inherits = FALSE
)
empirical_statistic <- get(
  ".studentized_statistic", envir = specifier$.nb_signflip_env,
  inherits = FALSE
)

set.seed(config$seed)
seed_matrix <- matrix(
  sample.int(
    .Machine$integer.max,
    size = nrow(pilot_grid) * config$reps,
    replace = FALSE
  ),
  nrow = nrow(pilot_grid), ncol = config$reps
)

score_references <- function(data, observed_information,
                             observed_empirical) {
  cache <- data$score_cache
  B <- config$permutations
  n <- data$n
  batch_size <- 128L
  completed <- 0L
  hits <- c(
    information_down = 0L,
    empirical_down = 0L,
    information_up = 0L,
    empirical_up = 0L,
    information_absolute = 0L,
    empirical_absolute = 0L
  )

  set.seed(data$perm_seed)
  while (completed < B) {
    current_batch <- min(batch_size, B - completed)
    permutation_indices <- vapply(
      seq_len(current_batch),
      function(i) sample.int(n),
      FUN.VALUE = integer(n)
    )
    X_permuted <- matrix(
      data$X[permutation_indices], nrow = n, ncol = current_batch
    )
    projection_coef <- cache$gram_inverse %*%
      crossprod(cache$D, X_permuted * cache$W)
    X_tilde_permuted <- X_permuted - cache$D %*% projection_coef
    score_permuted <- X_tilde_permuted * cache$outcome_factor
    reference_information <- information_statistic(
      score_permuted, X_tilde_permuted, cache$W
    )
    reference_empirical <- empirical_statistic(score_permuted)

    hits[["information_down"]] <- hits[["information_down"]] +
      sum(reference_information <= observed_information)
    hits[["empirical_down"]] <- hits[["empirical_down"]] +
      sum(reference_empirical <= observed_empirical)
    hits[["information_up"]] <- hits[["information_up"]] +
      sum(reference_information >= observed_information)
    hits[["empirical_up"]] <- hits[["empirical_up"]] +
      sum(reference_empirical >= observed_empirical)
    hits[["information_absolute"]] <-
      hits[["information_absolute"]] +
      sum(abs(reference_information) >= abs(observed_information))
    hits[["empirical_absolute"]] <- hits[["empirical_absolute"]] +
      sum(abs(reference_empirical) >= abs(observed_empirical))
    completed <- completed + current_batch
  }
  hits
}

run_one <- function(grid_index, replicate_index) {
  grid_row <- pilot_grid[grid_index, , drop = FALSE]
  data_seed <- seed_matrix[grid_index, replicate_index]
  set.seed(data_seed)
  arguments <- as.list(grid_row[1L, generate_args, drop = FALSE])
  arguments$B_resample <- config$permutations
  data <- do.call(specifier$generate_data, arguments)

  base_row <- data.frame(
    grid_id = grid_row$grid_id,
    replicate = replicate_index,
    direction = config$direction,
    gamma_signal = grid_row$gamma_signal,
    fold_change = grid_row$fold_change,
    data_seed = data_seed,
    permutation_seed = data$perm_seed,
    success = isTRUE(data$success),
    failure_reason = if (isTRUE(data$success)) NA_character_ else
      as.character(data$failure_reason),
    degenerate = FALSE,
    theta_hat = data$theta_hat,
    n_exposed = data$n_exposed,
    n_exposed_nonzero = data$n_exposed_nonzero,
    information_statistic = NA_real_,
    empirical_statistic = NA_real_,
    information_down_hits = NA_integer_,
    empirical_down_hits = NA_integer_,
    information_up_hits = NA_integer_,
    empirical_up_hits = NA_integer_,
    information_absolute_hits = NA_integer_,
    empirical_absolute_hits = NA_integer_,
    p_down_barry = NA_real_,
    p_down_empirical = NA_real_,
    p_up_barry = NA_real_,
    p_up_empirical = NA_real_,
    p_two_sided_barry = NA_real_,
    p_two_sided_empirical = NA_real_,
    score_eff_n = if (is.null(data$score_cache)) NA_real_ else
      data$score_cache$score_eff_n,
    max_score_share = if (is.null(data$score_cache)) NA_real_ else
      data$score_cache$max_score_share,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(data$success)) {
    return(base_row)
  }
  if (isTRUE(data$score_cache$degenerate)) {
    base_row$degenerate <- TRUE
    base_row$information_statistic <- 0
    base_row$empirical_statistic <- 0
    hit_columns <- grep("_hits$", names(base_row), value = TRUE)
    base_row[hit_columns] <- config$permutations
    p_columns <- grep("^p_", names(base_row), value = TRUE)
    base_row[p_columns] <- 1
    return(base_row)
  }

  cache <- data$score_cache
  observed_information <- information_statistic(
    cache$score_contrib, cache$X_tilde, cache$W
  )[1L]
  observed_empirical <- empirical_statistic(cache$score_contrib)[1L]
  hits <- score_references(
    data, observed_information, observed_empirical
  )
  denominator <- config$permutations + 1

  base_row$information_statistic <- observed_information
  base_row$empirical_statistic <- observed_empirical
  base_row$information_down_hits <- hits[["information_down"]]
  base_row$empirical_down_hits <- hits[["empirical_down"]]
  base_row$information_up_hits <- hits[["information_up"]]
  base_row$empirical_up_hits <- hits[["empirical_up"]]
  base_row$information_absolute_hits <- hits[["information_absolute"]]
  base_row$empirical_absolute_hits <- hits[["empirical_absolute"]]
  base_row$p_down_barry <- (1 + hits[["information_down"]]) / denominator
  base_row$p_down_empirical <- (1 + hits[["empirical_down"]]) / denominator
  base_row$p_up_barry <- (1 + hits[["information_up"]]) / denominator
  base_row$p_up_empirical <- (1 + hits[["empirical_up"]]) / denominator
  base_row$p_two_sided_barry <-
    (1 + hits[["information_absolute"]]) / denominator
  base_row$p_two_sided_empirical <-
    (1 + hits[["empirical_absolute"]]) / denominator
  base_row
}

cat(sprintf(
  paste0(
    "CI rare-NB one-sided %s-tail pilot: 5 signals x %d reps, n=400, ",
    "B=%d, cores=%d, seed=%d\n"
  ),
  config$direction, config$reps, config$permutations, config$cores,
  config$seed
))
start_time <- Sys.time()
all_results <- vector("list", nrow(pilot_grid))
for (grid_index in seq_len(nrow(pilot_grid))) {
  grid_start <- Sys.time()
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
  all_results[[grid_index]] <- do.call(rbind, rows)
  cat(sprintf(
    "  gamma=%7.4f complete in %.1f seconds\n",
    pilot_grid$gamma_signal[grid_index],
    as.numeric(difftime(Sys.time(), grid_start, units = "secs"))
  ))
  flush.console()
}
replicates <- do.call(rbind, all_results)
rownames(replicates) <- NULL

valid <- replicates$success &
  !is.na(replicates[[primary_barry_column]]) &
  !is.na(replicates[[primary_empirical_column]])
if (!all(valid)) {
  warning(sum(!valid), " pilot rows failed and are excluded from summaries.")
}

make_method_summary <- function(gamma, method, p_column, tail) {
  keep <- valid & replicates$gamma_signal == gamma
  p_values <- replicates[[p_column]][keep]
  rate <- mean(p_values <= config$alpha)
  data.frame(
    gamma_signal = gamma,
    fold_change = exp(gamma),
    tail = tail,
    method = method,
    n_valid = length(p_values),
    n_rejected = sum(p_values <= config$alpha),
    rejection_rate = rate,
    monte_carlo_se = sqrt(rate * (1 - rate) / length(p_values)),
    stringsAsFactors = FALSE
  )
}

method_summary <- do.call(rbind, lapply(expected_gamma, function(gamma) {
  rbind(
    make_method_summary(
      gamma, "barry", primary_barry_column, primary_tail
    ),
    make_method_summary(
      gamma, "empirical", primary_empirical_column, primary_tail
    ),
    make_method_summary(
      gamma, "barry", "p_two_sided_barry", "two_sided"
    ),
    make_method_summary(
      gamma, "empirical", "p_two_sided_empirical", "two_sided"
    )
  )
}))

paired_summary <- do.call(rbind, lapply(expected_gamma, function(gamma) {
  keep <- valid & replicates$gamma_signal == gamma
  barry <- replicates[[primary_barry_column]][keep] <= config$alpha
  empirical <- replicates[[primary_empirical_column]][keep] <= config$alpha
  empirical_only <- sum(empirical & !barry)
  barry_only <- sum(barry & !empirical)
  difference <- as.numeric(empirical) - as.numeric(barry)
  discordance_p <- if (empirical_only + barry_only == 0L) {
    1
  } else {
    stats::binom.test(
      empirical_only, empirical_only + barry_only,
      p = 0.5, alternative = "two.sided"
    )$p.value
  }
  data.frame(
    gamma_signal = gamma,
    n_valid = length(barry),
    both = sum(barry & empirical),
    barry_only = barry_only,
    empirical_only = empirical_only,
    neither = sum(!barry & !empirical),
    barry_rejection_rate = mean(barry),
    empirical_rejection_rate = mean(empirical),
    empirical_minus_barry = mean(difference),
    paired_monte_carlo_se = stats::sd(difference) / sqrt(length(difference)),
    exact_discordance_p = discordance_p,
    stringsAsFactors = FALSE
  )
}))

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
stem_prefix <- if (config$direction == "down") {
  "nb-rare-one-sided-pilot"
} else {
  "nb-rare-one-sided-right-signal-pilot"
}
stem <- sprintf(
  "%s-R%d-B%d-seed%d", stem_prefix,
  config$reps, config$permutations, config$seed
)
paths <- list(
  rds = file.path(config$output_dir, paste0(stem, ".rds")),
  replicates = file.path(
    config$output_dir, paste0(stem, "-replicates.csv")
  ),
  summary = file.path(config$output_dir, paste0(stem, "-summary.csv")),
  paired = file.path(
    config$output_dir, paste0(stem, "-paired-summary.csv")
  )
)
if (any(file.exists(unlist(paths)))) {
  stop("Refusing to overwrite an existing pilot output for this configuration.",
       call. = FALSE)
}

saved <- list(
  configuration = c(
    config,
    list(
      n = 400L,
      null_design = "CI",
      scenario_id = "rare_nb",
      test = paste0("one-sided ", config$direction, "-tail"),
      calibration = "nominal unadjusted",
      formula = if (config$direction == "down") {
        "(1 + sum(T_perm <= T_obs)) / (B + 1)"
      } else {
        "(1 + sum(T_perm >= T_obs)) / (B + 1)"
      },
      specifier = normalizePath(specifier_file),
      started = format(start_time, tz = "UTC", usetz = TRUE),
      completed = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  ),
  grid = pilot_grid,
  replicates = replicates,
  method_summary = method_summary,
  paired_summary = paired_summary
)
saveRDS(saved, paths$rds)
write.csv(replicates, paths$replicates, row.names = FALSE, na = "")
write.csv(method_summary, paths$summary, row.names = FALSE, na = "")
write.csv(paired_summary, paths$paired, row.names = FALSE, na = "")

cat(sprintf(
  paste0(
    "Completed in %.1f seconds.\n\nOne-sided %s-tail nominal ",
    "alpha=%.3f:\n"
  ),
  as.numeric(difftime(Sys.time(), start_time, units = "secs")),
  config$direction, config$alpha
))
print(paired_summary, row.names = FALSE, digits = 4)
cat("\nSaved:\n", paste0("  ", normalizePath(unlist(paths)),
                           collapse = "\n"), "\n", sep = "")
