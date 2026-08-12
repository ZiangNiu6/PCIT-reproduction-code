#!/usr/bin/env Rscript
###############################################################################
# Merge per-setup results from ccle-typeI-synthY-30targets.R (one job per
# (target, s) pair) into a single combined file.
#
# Inputs and outputs live in $LOCAL_PDRT_DATA_DIR/ccle:
#   ccle-typeI-synthY-30targets-results-{TARGET}-s{S}.rds
#   ccle-typeI-synthY-30targets-results.rds
#   ccle-typeI-synthY-30targets-rejection-rates.csv
#
# To know which setups are expected, this script re-derives the 60-element
# grid the same way the runner does: it samples 10 targets from each frequency
# tercile of the screened mutations using SAMPLE_SEED = 2026.
#
# Usage: Rscript combine-ccle-results.R
###############################################################################

S_VALUES    <- c(5L, 10L)
N_PER_BIN   <- 10L
SAMPLE_SEED <- 2026L

data_root <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
if (!nzchar(data_root)) {
  stop(
    "LOCAL_PDRT_DATA_DIR must point to the external data/results directory.",
    call. = FALSE
  )
}
data_root <- path.expand(data_root)
if (!dir.exists(data_root)) {
  stop("LOCAL_PDRT_DATA_DIR does not exist: ", data_root, call. = FALSE)
}
data_root <- normalizePath(data_root, mustWork = TRUE)
bpath <- file.path(data_root, "ccle")

# ---- Re-derive expected (target, s) grid ----
in_dir <- bpath
input_files <- file.path(in_dir, c("mutation.txt", "response.csv"))
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Missing CCLE input file(s):\n  ",
    paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}

mutations <- read.delim(file.path(in_dir, "mutation.txt"), sep = "\t")
mutations <- mutations[grep("_MUT$", mutations$Name), ]

response <- read.csv(file.path(in_dir, "response.csv"))
response <- response[response$Compound == "PLX4720", ]

idx  <- match(response$CCLE.Cell.Line.Name, colnames(mutations))
resp <- response[!is.na(idx), ]
mut  <- mutations[, c(1, idx[!is.na(idx)])]
rownames(mut) <- mut$Name
mut$Name <- NULL

Y_real <- resp$Amax
cors <- suppressWarnings(apply(mut, 1, function(m) {
  tc <- cor(Y_real, m); ifelse(!is.na(tc), tc, 0)
}))
screen <- cors[abs(cors) > 0.05 & !is.na(cors)]
dat <- data.frame(Y_real = Y_real, t(mut[names(screen), ]))

M_full <- as.matrix(dat[, setdiff(colnames(dat), "Y_real")])
freqs  <- colMeans(M_full)

qs <- quantile(freqs, probs = c(1/3, 2/3), names = FALSE)
pool_low  <- sort(names(freqs)[freqs <= qs[1]])
pool_med  <- sort(names(freqs)[freqs >  qs[1] & freqs <= qs[2]])
pool_high <- sort(names(freqs)[freqs >  qs[2]])

if (any(c(length(pool_low), length(pool_med), length(pool_high)) < N_PER_BIN)) {
  stop(
    "A frequency bin has fewer than ", N_PER_BIN, " mutations: ",
    sprintf(
      "low=%d, med=%d, high=%d",
      length(pool_low), length(pool_med), length(pool_high)
    ),
    call. = FALSE
  )
}

set.seed(SAMPLE_SEED)
targets_low  <- sort(sample(pool_low,  N_PER_BIN))
targets_med  <- sort(sample(pool_med,  N_PER_BIN))
targets_high <- sort(sample(pool_high, N_PER_BIN))

TARGETS_DF <- rbind(
  data.frame(target = targets_low,  tag = "low_sparsity",    stringsAsFactors = FALSE),
  data.frame(target = targets_med,  tag = "medium_sparsity", stringsAsFactors = FALSE),
  data.frame(target = targets_high, tag = "high_sparsity",   stringsAsFactors = FALSE)
)
TARGETS_DF$freq <- freqs[TARGETS_DF$target]

SETUPS <- rbind(
  cbind(TARGETS_DF, s = S_VALUES[1L]),
  cbind(TARGETS_DF, s = S_VALUES[2L])
)
rownames(SETUPS) <- NULL

# ---- Load each per-setup file ----
expected_files <- file.path(bpath,
  sprintf("ccle-typeI-synthY-30targets-results-%s-s%d.rds",
          SETUPS$target, SETUPS$s))
missing <- expected_files[!file.exists(expected_files)]
if (length(missing))
  stop("Missing per-setup files (", length(missing), " of ",
       length(expected_files), "):\n  ",
       paste(head(missing, 20), collapse = "\n  "),
       if (length(missing) > 20) "\n  ..." else "")

results <- list()
expected_methods <- c("asym_gcm", "asym_pcm", "perm_gcm", "perm_pcm", "naive")
for (i in seq_along(expected_files)) {
  obj <- readRDS(expected_files[i])
  if (!is.list(obj) || length(obj) != 1L || !is.list(obj[[1L]])) {
    stop("Malformed per-setup result: ", expected_files[i], call. = FALSE)
  }

  ## Each per-setup file is a length-1 list containing the result record.
  rec <- obj[[1L]]
  required_fields <- c(
    "setup_idx", "target", "sparsity_tag", "freq", "s", "R", "pvals",
    "rejection_rate", "na_rate"
  )
  missing_fields <- setdiff(required_fields, names(rec))
  if (length(missing_fields)) {
    stop(
      "Missing field(s) in ", expected_files[i], ": ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (!identical(as.integer(rec$setup_idx), as.integer(i)) ||
      !identical(as.character(rec$target), as.character(SETUPS$target[i])) ||
      !identical(as.integer(rec$s), as.integer(SETUPS$s[i])) ||
      !identical(as.character(rec$sparsity_tag), as.character(SETUPS$tag[i])) ||
      !isTRUE(all.equal(
        as.numeric(rec$freq), as.numeric(SETUPS$freq[i]), tolerance = 1e-12
      ))) {
    stop("Setup metadata mismatch in ", expected_files[i], call. = FALSE)
  }
  if (!identical(as.integer(rec$R), 500L) ||
      !is.matrix(rec$pvals) ||
      nrow(rec$pvals) != 500L ||
      !identical(colnames(rec$pvals), expected_methods) ||
      !identical(names(rec$rejection_rate), expected_methods) ||
      !identical(names(rec$na_rate), expected_methods)) {
    stop("Method schema mismatch in ", expected_files[i], call. = FALSE)
  }
  if (any(!is.finite(rec$pvals)) ||
      any(rec$pvals < 0 | rec$pvals > 1) ||
      any(!is.finite(rec$rejection_rate)) ||
      any(!is.finite(rec$na_rate)) ||
      any(rec$na_rate != 0)) {
    stop("Invalid or missing p-value in ", expected_files[i], call. = FALSE)
  }
  recomputed_rejection <- colMeans(rec$pvals <= 0.05)
  if (!isTRUE(all.equal(
    unname(rec$rejection_rate), unname(recomputed_rejection),
    tolerance = 1e-12
  ))) {
    stop("Stored rejection rates disagree with p-values in ",
         expected_files[i], call. = FALSE)
  }

  key <- sprintf("%s_s%d", rec$target, rec$s)
  if (key %in% names(results)) {
    stop("Duplicate setup key: ", key, call. = FALSE)
  }
  results[[key]] <- rec
}

saveRDS(results, file.path(bpath, "ccle-typeI-synthY-30targets-results.rds"))

# ---- Long-format rejection table ----
rej_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(target = r$target,
             sparsity_tag = r$sparsity_tag,
             freq = r$freq,
             s = r$s,
             method = names(r$rejection_rate),
             rejection_rate = unname(r$rejection_rate),
             na_rate = unname(r$na_rate[names(r$rejection_rate)]),
             stringsAsFactors = FALSE)
}))
rownames(rej_df) <- NULL
write.csv(rej_df,
          file.path(bpath, "ccle-typeI-synthY-30targets-rejection-rates.csv"),
          row.names = FALSE)

cat("=== Combined 30-targets summary (", length(results), " setups) ===\n", sep = "")
cat(sprintf("Mean rejection rate per (method, s):\n"))
agg <- aggregate(rejection_rate ~ method + s, data = rej_df, FUN = mean)
print(agg)
cat("\nSaved:\n")
cat("  ", file.path(bpath, "ccle-typeI-synthY-30targets-results.rds"), "\n")
cat("  ", file.path(bpath, "ccle-typeI-synthY-30targets-rejection-rates.csv"), "\n")
