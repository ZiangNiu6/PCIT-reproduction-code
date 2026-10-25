#!/usr/bin/env Rscript

# Compare four nuisance fits on the same NB datasets (constant true size):
# formula Poisson GLM, SCEPTRE's Poisson glm.fit stage alone, SCEPTRE's full
# response precomputation (Poisson coefficients plus estimated NB size), and
# MASS::glm.nb (alternating mean and size estimation).
# One sample size per R process. Usage: Rscript --vanilla script theta n reps
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected theta, n, repetitions")
theta <- as.numeric(args[[1L]])
n <- as.integer(args[[2L]])
reps <- as.integer(args[[3L]])
if (!is.finite(theta) || theta <= 0 ||
    is.na(n) || n < 1L || is.na(reps) || reps < 2L) {
  stop("Invalid benchmark parameters")
}

set.seed(20261019L + n + as.integer(1000 * theta))
datasets <- replicate(reps, {
  Z <- stats::rnorm(n)
  Y <- stats::rnbinom(n, mu = 10 * exp(0.5 * Z), size = theta)
  Xmat <- cbind(Intercept = 1, Z = Z)
  list(Y = Y, Z = Z, Xmat = Xmat)
}, simplify = FALSE)
sceptre_precomp <- get(
  "perform_response_precomputation", envir = asNamespace("sceptre")
)
control <- stats::glm.control(maxit = 100L)

fitters <- list(
  poisson_glm = function(i) {
    Y <- datasets[[i]]$Y
    Z <- datasets[[i]]$Z
    stats::glm(Y ~ Z, family = stats::poisson(link = "log"),
               control = control)
  },
  sceptre_poisson_only = function(i) {
    stats::glm.fit(
      x = datasets[[i]]$Xmat, y = datasets[[i]]$Y,
      family = stats::poisson(link = "log")
    )
  },
  sceptre = function(i) {
    sceptre_precomp(
      expressions = datasets[[i]]$Y,
      covariate_matrix = datasets[[i]]$Xmat
    )
  },
  mass_glm_nb = function(i) {
    Y <- datasets[[i]]$Y
    Z <- datasets[[i]]$Z
    MASS::glm.nb(Y ~ Z, control = control)
  }
)
invisible(lapply(fitters, function(f) suppressWarnings(f(1L))))

is_bad <- function(fit, method) {
  if (inherits(fit, "error")) return(TRUE)
  if (method == "sceptre") {
    return(!is.finite(fit$theta) || fit$theta <= 0 ||
           any(!is.finite(fit$fitted_coefs)))
  }
  !isTRUE(fit$converged) ||
    (method == "mass_glm_nb" && !is.null(fit$th.warn))
}

elapsed <- setNames(vector("list", length(fitters)), names(fitters))
failed <- setNames(integer(length(fitters)), names(fitters))
for (pass in seq_len(2L)) {
  order <- if (pass == 1L) names(fitters) else rev(names(fitters))
  for (method in order) {
    gc(FALSE)
    fitter <- fitters[[method]]
    bad <- 0L
    if (n >= 40000L) {
      trial_time <- numeric(reps)
      for (i in seq_len(reps)) {
        gc(FALSE)
        trial_time[[i]] <- system.time({
          fit <- tryCatch(suppressWarnings(fitter(i)), error = identity)
          if (is_bad(fit, method)) bad <- bad + 1L
        })[["elapsed"]]
        rm(fit)
        gc(FALSE)
      }
      timing <- sum(trial_time)
    } else {
      timing <- system.time({
        for (i in seq_len(reps)) {
          fit <- tryCatch(suppressWarnings(fitter(i)), error = identity)
          if (is_bad(fit, method)) bad <- bad + 1L
        }
      })[["elapsed"]]
    }
    elapsed[[method]] <- c(elapsed[[method]], timing)
    if (pass == 1L) failed[[method]] <- bad
  }
}

ms <- vapply(elapsed, mean, numeric(1L)) * 1000 / reps
row <- data.frame(
  theta = theta, n = n, reps_per_pass = reps,
  poisson_glm_ms = ms[["poisson_glm"]],
  sceptre_poisson_only_ms = ms[["sceptre_poisson_only"]],
  sceptre_ms = ms[["sceptre"]],
  mass_glm_nb_ms = ms[["mass_glm_nb"]],
  sceptre_to_poisson_only = ms[["sceptre"]] / ms[["sceptre_poisson_only"]],
  mass_to_sceptre = ms[["mass_glm_nb"]] / ms[["sceptre"]],
  failed_poisson_glm = failed[["poisson_glm"]],
  failed_sceptre_poisson_only = failed[["sceptre_poisson_only"]],
  failed_sceptre = failed[["sceptre"]],
  failed_mass_glm_nb = failed[["mass_glm_nb"]]
)
utils::write.table(row, file = stdout(), sep = ",", row.names = FALSE,
                   col.names = FALSE, quote = FALSE)
