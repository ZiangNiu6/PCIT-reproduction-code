###############################################################################
# Focused parameter grid: five score calibrations under sparse exposure and
# variance misspecification
#
# Each MI/CI figure has three rows and five signal points per row.  The first
# point is the null (|gamma| = 0), followed by four equally spaced
# downregulation alternatives.  The first two rows use a correctly specified,
# high-dispersion NB model with common or rare binary X.  The third row keeps
# the conditional mean exactly specified but makes the outcome variance vary
# nonmultiplicatively across a binary Z stratum while fitting a Poisson working
# score model.
###############################################################################

n <- 400L
B_resample <- 1999L
alpha <- 0.05
ci_x_z_slope <- 0.5

scenario_levels <- c("common_nb", "rare_nb", "severe_variance")

# Five total points per row: the null and four alternatives.  The endpoints
# target about 80% size-adjusted rejection for common X, about 40% for rare X,
# and about 80% in each variance-stress design.  The variance-stress endpoints
# differ between MI and CI because their binary-stratum exposure information is
# intentionally different.
signal_endpoints <- data.frame(
  null_design = rep(c("MI", "CI"), each = 3L),
  scenario_id = rep(scenario_levels, times = 2L),
  endpoint = c(0.45, 0.85, 0.105, 0.45, 0.85, 0.18),
  stringsAsFactors = FALSE
)

calibrate_logistic_intercept <- function(target_prevalence, z_slope) {
  marginal_difference <- function(intercept) {
    stats::integrate(
      function(z) stats::plogis(intercept + z_slope * z) * stats::dnorm(z),
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-12
    )$value - target_prevalence
  }

  stats::uniroot(
    marginal_difference,
    interval = c(-20, 20),
    tol = 1e-12
  )$root
}

marginal_logistic_probability <- function(intercept, z_slope) {
  stats::integrate(
    function(z) stats::plogis(intercept + z_slope * z) * stats::dnorm(z),
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-12
  )$value
}

make_correct_nb_design <- function(scenario_id, x_prevalence, null_design) {
  x_slope <- if (null_design == "MI") 0 else ci_x_z_slope
  x_intercept <- if (null_design == "MI") {
    stats::qlogis(x_prevalence)
  } else {
    calibrate_logistic_intercept(x_prevalence, x_slope)
  }

  data.frame(
    null_design = null_design,
    scenario_id = scenario_id,
    z_distribution = "normal",
    z_probability = NA_real_,
    x_prevalence = x_prevalence,
    x_intercept = x_intercept,
    x_z_slope = x_slope,
    x_probability_z0 = NA_real_,
    x_probability_z1 = NA_real_,
    mu0 = 10,
    beta_z = 0.5,
    variance_model = "constant_nb",
    working_family = "negative_binomial",
    nb_dispersion = 2,
    nb_size = 0.5,
    variance_factor_z0 = NA_real_,
    variance_factor_z1 = NA_real_,
    stringsAsFactors = FALSE
  )
}

# The severe-variance row adapts the manuscript's constructed
# nonmultiplicative quasi-Poisson example to a binary exposure.  With
# P(Z=1)=0.75, both MI and CI have marginal P(X=1)=0.14.  Under CI, the
# stratum probabilities (0.40, 0.053333...) make Var(X|Z) align with the
# high-variance outcome stratum, separating the sampling and Barry-permutation
# information scales without making the robust procedures conservative.
make_severe_variance_design <- function(null_design) {
  z_probability <- 0.75
  x_probability_z0 <- if (null_design == "MI") 0.14 else 0.40
  x_probability_z1 <- if (null_design == "MI") {
    0.14
  } else {
    (0.14 - (1 - z_probability) * x_probability_z0) / z_probability
  }
  x_prevalence <-
    (1 - z_probability) * x_probability_z0 +
    z_probability * x_probability_z1

  data.frame(
    null_design = null_design,
    scenario_id = "severe_variance",
    z_distribution = "binary",
    z_probability = z_probability,
    x_prevalence = x_prevalence,
    x_intercept = NA_real_,
    x_z_slope = NA_real_,
    x_probability_z0 = x_probability_z0,
    x_probability_z1 = x_probability_z1,
    mu0 = 20,
    beta_z = log(1.5),
    variance_model = "stratum_quasi_poisson",
    working_family = "poisson",
    nb_dispersion = NA_real_,
    nb_size = NA_real_,
    variance_factor_z0 = 5,
    variance_factor_z1 = 1,
    stringsAsFactors = FALSE
  )
}

design_grid <- do.call(
  rbind,
  c(
    lapply(c("MI", "CI"), function(null_design) {
      make_correct_nb_design("common_nb", 0.50, null_design)
    }),
    lapply(c("MI", "CI"), function(null_design) {
      make_correct_nb_design("rare_nb", 0.025, null_design)
    }),
    lapply(c("MI", "CI"), make_severe_variance_design)
  )
)

parameter_grid <- do.call(
  rbind,
  lapply(seq_len(nrow(design_grid)), function(index) {
    design <- design_grid[index, , drop = FALSE]
    endpoint <- signal_endpoints$endpoint[
      signal_endpoints$null_design == design$null_design &
        signal_endpoints$scenario_id == design$scenario_id
    ]
    if (length(endpoint) != 1L) {
      stop("Each design/scenario combination must have one signal endpoint")
    }
    magnitudes <- seq(0, endpoint, length.out = 5L)
    design[rep(1L, length(magnitudes)), , drop = FALSE] |>
      transform(gamma_signal = -magnitudes)
  })
)

parameter_grid$nuisance_mode <- "fitted"
parameter_grid$fold_change <- exp(parameter_grid$gamma_signal)
parameter_grid$n <- n
parameter_grid$B_resample <- B_resample
parameter_grid$alpha <- alpha
parameter_grid$ground_truth <- parameter_grid$gamma_signal != 0

# Keep a stable order: design, scenario row, and increasing signal magnitude.
parameter_grid <- parameter_grid[
  order(
    match(parameter_grid$null_design, c("MI", "CI")),
    match(parameter_grid$scenario_id, scenario_levels),
    abs(parameter_grid$gamma_signal)
  ),
]
rownames(parameter_grid) <- NULL
parameter_grid$grid_id <- seq_len(nrow(parameter_grid))

parameter_grid <- parameter_grid[, c(
  "grid_id", "n", "null_design", "scenario_id", "nuisance_mode",
  "z_distribution", "z_probability",
  "x_prevalence", "x_intercept", "x_z_slope",
  "x_probability_z0", "x_probability_z1",
  "mu0", "beta_z", "variance_model", "working_family",
  "nb_dispersion", "nb_size",
  "variance_factor_z0", "variance_factor_z1",
  "gamma_signal", "fold_change", "B_resample", "alpha", "ground_truth"
)]

expected_rows <- 2L * length(scenario_levels) * 5L
if (nrow(parameter_grid) != expected_rows ||
    any(parameter_grid$nuisance_mode != "fitted")) {
  stop("Unexpected parameter-grid size or nuisance mode")
}
factor_columns <- c("null_design", "scenario_id", "gamma_signal")
if (anyDuplicated(parameter_grid[, factor_columns])) {
  stop("Duplicated design/scenario/signal cell")
}
if (!setequal(unique(parameter_grid$scenario_id), scenario_levels) ||
    any(parameter_grid$gamma_signal > 0) ||
    any(parameter_grid$ground_truth != (parameter_grid$gamma_signal != 0))) {
  stop("Focused scenario or signal constraints were violated")
}

correct_nb_rows <- parameter_grid$scenario_id %in% c("common_nb", "rare_nb")
if (any(!is.finite(parameter_grid$nb_dispersion[correct_nb_rows])) ||
    any(!is.finite(parameter_grid$nb_size[correct_nb_rows])) ||
    max(abs(
      parameter_grid$nb_size[correct_nb_rows] *
        parameter_grid$nb_dispersion[correct_nb_rows] - 1
    )) > 1e-12 ||
    any(parameter_grid$working_family[correct_nb_rows] != "negative_binomial") ||
    any(parameter_grid$variance_model[correct_nb_rows] != "constant_nb")) {
  stop("Correct-NB scenario parameters are inconsistent")
}

normal_designs <- unique(parameter_grid[
  correct_nb_rows,
  c("null_design", "scenario_id", "x_prevalence", "x_intercept", "x_z_slope")
])
normal_prevalence <- mapply(
  marginal_logistic_probability,
  intercept = normal_designs$x_intercept,
  z_slope = normal_designs$x_z_slope
)
if (max(abs(normal_prevalence - normal_designs$x_prevalence)) > 1e-10) {
  stop("Continuous-Z exposure intercept calibration failed")
}

severe_rows <- parameter_grid$scenario_id == "severe_variance"
severe_designs <- unique(parameter_grid[
  severe_rows,
  c(
    "null_design", "z_probability", "x_prevalence",
    "x_probability_z0", "x_probability_z1", "working_family",
    "variance_model", "variance_factor_z0", "variance_factor_z1"
  )
])
severe_prevalence <- with(
  severe_designs,
  (1 - z_probability) * x_probability_z0 +
    z_probability * x_probability_z1
)
if (max(abs(severe_prevalence - severe_designs$x_prevalence)) > 1e-12 ||
    any(abs(severe_designs$x_prevalence - 0.14) > 1e-12) ||
    any(severe_designs$working_family != "poisson") ||
    any(severe_designs$variance_model != "stratum_quasi_poisson") ||
    any(severe_designs$variance_factor_z0 != 5) ||
    any(severe_designs$variance_factor_z1 != 1)) {
  stop("Severe-variance scenario parameters are inconsistent")
}

design_scenario_groups <- split(
  parameter_grid,
  interaction(
    parameter_grid$null_design,
    parameter_grid$scenario_id,
    drop = TRUE
  )
)
for (group in design_scenario_groups) {
  endpoint <- signal_endpoints$endpoint[
    signal_endpoints$null_design == group$null_design[[1L]] &
      signal_endpoints$scenario_id == group$scenario_id[[1L]]
  ]
  expected_magnitudes <- seq(0, endpoint, length.out = 5L)
  if (nrow(group) != 5L ||
      sum(!group$ground_truth) != 1L ||
      sum(group$ground_truth) != 4L ||
      !isTRUE(all.equal(
        sort(abs(group$gamma_signal)),
        expected_magnitudes,
        tolerance = 1e-12
      ))) {
    stop("Each design/scenario cell must have one null and four alternatives")
  }
}

# This file is sourced by nb-permutation-signflip.R. Keeping the grid as code
# removes the old generated code/parameter-grid layer and makes the simulation
# specifier self-contained.
if (sys.nframe() == 0L) {
  cat(
    nrow(parameter_grid), "rows:",
    "2 null designs x 3 scenario rows x 5 signal values\n"
  )
  print(unique(parameter_grid[, c(
    "null_design", "scenario_id", "z_distribution", "x_prevalence",
    "working_family", "variance_model", "nb_dispersion", "nb_size",
    "variance_factor_z0", "variance_factor_z1"
  )]))
}
