###############################################################################
# Parameter grid for the three-row continuous-regression comparison.
#
# Let s(Z) = sin(2*pi*Z), with Z, eps independent standard Gaussian.
#
#   1. X independent of Z:
#        X = xi,  Y = s(Z) + gamma * X^2 + eps.
#   2. Additive Gaussian:
#        X = s(Z) + xi,  Y = s(Z) + gamma * X^2 + eps.
#   3. Heteroscedastic:
#        xi + 1 ~ Exp(1),  X = s(Z) * (1 - xi),
#        Y = s(Z) + gamma * X + eps.
#
# The first two quadratic alternatives are deliberately PCM-favourable.  The
# heteroscedastic alternative is aligned with the GCM residual-covariance
# moment and supplies the complementary GCM-favourable example.
#
# Every displayed row has its own gamma=0 sample-size grid.  Its largest-n null
# cell is shared with the first point of its five-point power curve, so the
# final grid contains 27 unique cells rather than duplicated null cells.
###############################################################################

n_xindep_null <- c(20L, 40L, 60L, 80L, 100L)
n_ci_null <- c(200L, 300L, 400L, 500L, 600L)

signal_grids <- list(
  setting1_xindep = seq(0, 0.8, length.out = 5L),
  setting1 = seq(0, 0.3, length.out = 5L),
  setting2 = seq(0, 0.2, length.out = 5L)
)

make_design_grid <- function(dgp, null_n, power_n, signal_grid) {
  null_grid <- data.frame(
    n = as.integer(null_n),
    dgp = dgp,
    gamma_signal = 0,
    stringsAsFactors = FALSE
  )
  alternative_grid <- data.frame(
    n = as.integer(power_n),
    dgp = dgp,
    gamma_signal = signal_grid[signal_grid > 0],
    stringsAsFactors = FALSE
  )
  rbind(null_grid, alternative_grid)
}

parameter_grid <- rbind(
  make_design_grid(
    "setting1_xindep", n_xindep_null, 100L,
    signal_grids$setting1_xindep
  ),
  make_design_grid(
    "setting1", n_ci_null, 600L,
    signal_grids$setting1
  ),
  make_design_grid(
    "setting2", n_ci_null, 600L,
    signal_grids$setting2
  )
)

parameter_grid$B_perm <- 499L
parameter_grid$alpha <- 0.05
parameter_grid$ground_truth <- parameter_grid$gamma_signal != 0
parameter_grid$grid_id <- seq_len(nrow(parameter_grid))
parameter_grid <- parameter_grid[, c(
  "grid_id", "n", "dgp", "gamma_signal",
  "B_perm", "alpha", "ground_truth"
)]

expected_dgps <- names(signal_grids)
if (nrow(parameter_grid) != 27L ||
    anyDuplicated(parameter_grid[c("n", "dgp", "gamma_signal")]) ||
    !setequal(parameter_grid$dgp, expected_dgps) ||
    any(parameter_grid$ground_truth != (parameter_grid$gamma_signal != 0))) {
  stop("The pcm-additive parameter-grid contract was violated")
}

for (dgp_value in expected_dgps) {
  power_n <- if (identical(dgp_value, "setting1_xindep")) 100L else 600L
  found_signal <- parameter_grid$gamma_signal[
    parameter_grid$dgp == dgp_value & parameter_grid$n == power_n
  ]
  if (!isTRUE(all.equal(
    sort(found_signal), signal_grids[[dgp_value]], tolerance = 1e-12
  ))) {
    stop("Unexpected five-point signal grid for ", dgp_value)
  }
}

# This file is sourced by pcm-additive.R. Keeping the grid as code removes the
# old generated code/parameter-grid layer and makes the specifier portable.
if (sys.nframe() == 0L) {
  cat("pcm-additive parameter grid (", nrow(parameter_grid), " rows):\n",
      sep = "")
  print(parameter_grid)
}
