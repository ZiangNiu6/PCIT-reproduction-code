###############################################################################
# Simulation specifier: three-row GCM/PCM continuous-regression comparison
#
# DGPs (Z1 and ε are independent standard Gaussian; ξ is row-specific):
#   setting1_xindep:
#       X = ξ,                         Y = sin(2πZ1) + γ X² + ε
#   setting1 (Additive Gaussian):
#       X = sin(2πZ1) + ξ,            Y = sin(2πZ1) + γ X² + ε
#   setting2 (Heteroscedastic; ξ+1~Exp(1)):
#       X = sin(2πZ1)(1-ξ),          Y = sin(2πZ1) + γ X + ε
#
# Each row has its own gamma=0 sample-size grid.  Power uses n=600 for
# settings 1-2 and n=100 for setting1_xindep.  The two quadratic rows
# demonstrate PCM's advantage, whereas the linear heteroscedastic alternative
# makes GCM powerful.
#
# Methods:
#   1. Perm GCM-RF   (pcmt::perm_gcm, adaptive h=15, two-sided)
#   2. Perm PCM-RF   (pcmt::perm_pcm, adaptive h=15, upper-tailed)
#   3. Asym GCM-RF   (comets::gcm, two-sided)
#   4. Asym PCM-RF   (comets::pcm, its native upper-tailed test)
###############################################################################

#####################
# 1. Parameter grid
#####################
specifier_file <- normalizePath(sys.frame(1L)$ofile)
sys.source(
  file.path(dirname(specifier_file), "pcm-additive-grid.R"),
  envir = environment()
)
parameter_grid <- as.data.frame(parameter_grid)

expected_dgps <- c("setting1_xindep", "setting1", "setting2")
if (nrow(parameter_grid) != 27L ||
    anyDuplicated(parameter_grid[c("n", "dgp", "gamma_signal")]) ||
    !setequal(parameter_grid$dgp, expected_dgps) ||
    any(parameter_grid$ground_truth != (parameter_grid$gamma_signal != 0))) {
  stop(
    "The pcm-additive parameter grid is incompatible with this specifier."
  )
}

#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  B = 1000,          # number of dataset replications per grid point
  seed = 2024
)

###################
# 3. Generate data
###################
generate_data <- function(n, dgp, gamma_signal, B_perm) {
  Z1 <- rnorm(n)
  Z <- matrix(Z1, ncol = 1)
  colnames(Z) <- "Z1"
  eps <- rnorm(n)

  if (dgp == "setting1") {
    xi <- rnorm(n)
    X <- sin(2 * pi * Z1) + xi
    Y <- sin(2 * pi * Z1) + gamma_signal * X^2 + eps
  } else if (dgp == "setting2") {
    xi <- rexp(n, 1) - 1  # ξ+1 ~ Exp(1)
    X <- sin(2 * pi * Z1) - sin(2 * pi * Z1) * xi
    Y <- sin(2 * pi * Z1) + gamma_signal * X + eps
  } else if (dgp == "setting1_xindep") {
    xi <- rnorm(n)
    X <- xi  # X independent of Z
    Y <- sin(2 * pi * Z1) + gamma_signal * X^2 + eps
  } else {
    stop("Unknown pcm-additive DGP: ", dgp)
  }

  list(X = X, Y = Y, Z = Z, n = n, B_perm = B_perm)
}

generate_data_function <- simulatr::simulatr_function(
  f = generate_data,
  arg_names = formalArgs(generate_data),
  loop = TRUE
)

######################
# 4. Method functions
######################

# ---- Method 1: Perm GCM-RF (adaptive) ----
perm_gcm_rf_f <- function(data) {
  suppressWarnings({
    res <- tryCatch(
      pcmt::perm_gcm(
        Y = data$Y, X = data$X, Z = data$Z,
        reg_YonZ = "rf", reg_XonZ = "rf",
        B = data$B_perm, h = 15L, alpha = 0.05,
        alternative = "two.sided"
      ),
      error = function(e) list(p.value = NA)
    )
    list(p.value = res$p.value)
  })
}

# ---- Method 2: Perm PCM-RF (adaptive) ----
perm_pcm_rf_f <- function(data) {
  suppressWarnings({
    res <- tryCatch(
      pcmt::perm_pcm(
        Y = data$Y, X = data$X, Z = data$Z,
        reg_YonXZ = "rf", reg_YonZ = "rf",
        reg_YhatonZ = "rf", reg_RonZ = "rf",
        B = data$B_perm, h = 15L, alpha = 0.05,
        alternative = "greater"
      ),
      error = function(e) list(p.value = NA)
    )
    list(p.value = res$p.value)
  })
}

# ---- Method 3: Asym GCM-RF ----
asym_gcm_rf_f <- function(data) {
  suppressWarnings({
    res <- tryCatch(
      comets::gcm(
        Y = data$Y, X = data$X, Z = data$Z,
        reg_YonZ = "rf", reg_XonZ = "rf",
        alternative = "two.sided"
      ),
      error = function(e) list(p.value = NA)
    )
    list(p.value = res$p.value)
  })
}

# ---- Method 4: Asym PCM-RF ----
asym_pcm_rf_f <- function(data) {
  suppressWarnings({
    res <- tryCatch(
      comets::pcm(
        Y = data$Y, X = data$X, Z = data$Z,
        reg_YonXZ = "rf", reg_YonZ = "rf",
        reg_YhatonZ = "rf", reg_RonZ = "rf"
      ),
      error = function(e) list(p.value = NA)
    )
    list(p.value = res$p.value)
  })
}

# wrap methods in simulatr_function
perm_gcm_rf_spec_f <- simulatr::simulatr_function(f = perm_gcm_rf_f, arg_names = character(0), loop = TRUE)
perm_pcm_rf_spec_f <- simulatr::simulatr_function(f = perm_pcm_rf_f, arg_names = character(0), loop = TRUE)
asym_gcm_rf_spec_f <- simulatr::simulatr_function(f = asym_gcm_rf_f, arg_names = character(0), loop = TRUE)
asym_pcm_rf_spec_f <- simulatr::simulatr_function(f = asym_pcm_rf_f, arg_names = character(0), loop = TRUE)

# assemble method functions
run_method_functions <- list(
  perm_gcm_rf = perm_gcm_rf_spec_f,
  perm_pcm_rf = perm_pcm_rf_spec_f,
  asym_gcm_rf = asym_gcm_rf_spec_f,
  asym_pcm_rf = asym_pcm_rf_spec_f
)

############################
# 5. Assemble specifier
############################
simulatr_spec <- simulatr::simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function,
  run_method_functions
)

# ############################### check results ##################################
# check_results_nofit <- simulatr::check_simulatr_specifier_object(simulatr_spec, B_in = 1)
