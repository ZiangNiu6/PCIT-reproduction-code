###############################################################################
# Parameter grid for the spline-nuisance companion cell.
#
# A single cell: the Heteroscedastic design at n = 200 under the null.  This
# exists to sit beside the random-forest row of the same design and sample
# size, so that the two differ only in the nuisance regression.  Nothing else
# is needed -- the surrounding rows come from pcm-additive.R.
###############################################################################

parameter_grid <- data.frame(
  grid_id = 1L,
  n = 200L,
  dgp = "setting2",
  gamma_signal = 0,
  B_perm = 499L,
  alpha = 0.05,
  ground_truth = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(parameter_grid) != 1L ||
    parameter_grid$dgp != "setting2" ||
    parameter_grid$n != 200L ||
    parameter_grid$gamma_signal != 0 ||
    parameter_grid$B_perm != 499L) {
  stop("The pcm-additive-gam parameter-grid contract was violated")
}

if (sys.nframe() == 0L) print(parameter_grid)
