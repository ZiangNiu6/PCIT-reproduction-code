#!/usr/bin/env Rscript
###############################################################################
# Local execution of the spline Section 5.2 specifier, through simulatr.
#
# The cluster route is
#   ./run_simulation.sh --sim_name pcm-additive-gam --results_dir DIR
# which drives the same specifier through Nextflow.  Nextflow is not required
# to reproduce these results: simulatr can evaluate the specifier directly, and
# kernel ridge with cached permutation refits is cheap enough to do so on one
# machine.  This script takes that route, so seeding and the output schema are
# simulatr's -- the results object has the same method/grid_id/run_id/output
# structure the pipeline writes and plot-pcm-additive.R consumes.
#
# Usage:
#   Rscript run-pcm-additive-gam-local.R [--B=1000] [--parallel=TRUE]
###############################################################################
suppressPackageStartupMessages({
  library(simulatr)
  library(pcmt)
})

script_dir <- {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
}
cli <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  hit <- grep(paste0("^--", flag, "="), cli, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", flag, "="), "", hit[[1L]])
}
B_reps   <- as.integer(get_opt("B", 1000))
use_par  <- as.logical(get_opt("parallel", "TRUE"))
# --pcm-refit=TRUE forces a genuine refit of the PCM projection learner on
# every permutation instead of reusing the frozen operator.
pcm_refit <- as.logical(get_opt("pcm-refit", "FALSE"))
options(pcmt.cached_refit_pcm = !isTRUE(pcm_refit))

source(file.path(script_dir, "sim-spec-formal", "pcm-additive-gam.R"),
       local = FALSE)

spec <- simulatr::update_B_sim_spec(simulatr_spec, B_reps)
cat(sprintf("pcm-additive-gam via simulatr: %d grid cells x B = %d%s%s\n\n",
            nrow(parameter_grid), B_reps,
            if (isTRUE(use_par)) ", parallel" else "",
            if (isTRUE(pcm_refit)) ", PCM genuine refit" else ""))

t0 <- proc.time()[["elapsed"]]
res <- simulatr::check_simulatr_specifier_object(
  spec, B_in = B_reps, return_data = FALSE, parallel = use_par
)
elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nelapsed: %.0fs\n", elapsed))

out_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
out_dir <- if (nzchar(out_root)) {
  file.path(out_root, "simulation-results", "pcm-additive-gam")
} else {
  file.path(getwd(), "results", "pcm-additive-gam")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(results = res$results, metrics = res$metrics,
             parameter_grid = parameter_grid,
             config = list(B = B_reps, elapsed_sec = elapsed,
                           pcm_genuine_refit = isTRUE(pcm_refit))),
        file.path(out_dir, paste0("pcm-additive-gam",
                                  if (isTRUE(pcm_refit)) "-pcmrefit" else "",
                                  "_results.rds")))
cat("Written to: ", out_dir, "\n", sep = "")
