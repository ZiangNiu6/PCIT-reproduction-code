#!/usr/bin/env Rscript
# Regression check: round-off in derived rates is accepted; altered estimates,
# reordered rows, and even tiny changes to exact design fields are rejected.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "validation-utils.R"))
reference <- data.frame(
  setup_idx = c(1, 2), effect_size = c(0, 0.5),
  method = c("asym_gcm", "perm_gcm"),
  mean_rejection_rate = c(0.0134, 0.05)
)
rounded <- reference
rounded$mean_rejection_rate[[1L]] <- 0.01339999999999999871769
stopifnot(ccle_same_ordered_frame(
  rounded, reference, "mean_rejection_rate"
))
changed <- reference
changed$mean_rejection_rate[[1L]] <- changed$mean_rejection_rate[[1L]] + 0.002
stopifnot(!ccle_same_ordered_frame(
  changed, reference, "mean_rejection_rate"
))
changed <- reference
changed$setup_idx[[1L]] <- changed$setup_idx[[1L]] + 1e-14
stopifnot(!ccle_same_ordered_frame(
  changed, reference, "mean_rejection_rate"
))
changed <- reference
changed$effect_size[[1L]] <- 1e-14
stopifnot(!ccle_same_ordered_frame(
  changed, reference, "mean_rejection_rate"
))
stopifnot(!ccle_same_ordered_frame(
  reference[2:1, ], reference, "mean_rejection_rate"
))
cat("CCLE derived-summary validation smoke test passed.\n")
