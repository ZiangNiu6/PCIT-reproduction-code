#!/usr/bin/env Rscript
# Exercise output verification with isolated valid, malformed and tampered files.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
package_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
checker <- file.path(package_root, "check_outputs.R")
source(file.path(package_root, "real-data-analysis", "method_config.R"))
stopifnot(requireNamespace("jsonlite", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))

run_tests <- function() {
  temporary_root <- tempfile("pcit-output-checks-")
  dir.create(temporary_root)
  on.exit(unlink(temporary_root, recursive = TRUE), add = TRUE)
  figure_dir <- file.path(temporary_root, "figures")
  result_dir <- file.path(temporary_root, "results")
  dir.create(figure_dir)
  dir.create(result_dir)
  figure_names <- c(
    "ccle-typeI-synthY-EIF5A-calibration.pdf",
    "illustrative-example-ci-misspecified.pdf",
    "nb-permutation-signflip-ci-focused-3x3.pdf", "pcm_additive.pdf",
    "ccle-typeI-synthY-30targets-typeI-combined.pdf",
    "frangieh-score5-fdr-power.pdf", "nb-permutation-signflip-mi-focused-3x3.pdf",
    "frangieh-score5-pvalue-histograms.pdf",
    "ccle-typeI-synthY-intro-calibration-power.pdf",
    "ccle-typeI-synthY-30targets-typeI-power-combined.pdf"
  )
  seed_pdf <- file.path(temporary_root, "seed.pdf")
  grDevices::pdf(seed_pdf)
  graphics::plot.new()
  grDevices::dev.off()
  stopifnot(all(file.copy(seed_pdf, file.path(figure_dir, figure_names))))

  run_check <- function(extra = character(), expected_status = 0L,
                        expected_message = NULL, base_arguments = NULL) {
    if (is.null(base_arguments)) {
      base_arguments <- c("--figure-dir", figure_dir, "--results-root", result_dir)
    }
    output <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      shQuote(c("--vanilla", checker, base_arguments, extra)),
      stdout = TRUE, stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    if (status != expected_status ||
        (!is.null(expected_message) && !any(grepl(expected_message, output, fixed = TRUE)))) {
      stop("Unexpected checker response (exit ", status, "):\n",
           paste(output, collapse = "\n"))
    }
  }
  run_check()
  run_check("--help", base_arguments = character(), expected_message = "Usage:")
  run_check("--unknown", expected_status = 1L, expected_message = "Unknown or repeated")
  run_check("--figure-dir", expected_status = 1L, expected_message = "Unknown or repeated")

  target_pdf <- file.path(figure_dir, figure_names[[1L]])
  writeLines("not a PDF file", target_pdf)
  run_check(expected_status = 1L, expected_message = "Invalid PDF header")
  file.copy(seed_pdf, target_pdf, overwrite = TRUE)
  connection <- file(target_pdf, "wb")
  close(connection)
  run_check(expected_status = 1L, expected_message = "Empty or unreadable")
  file.copy(seed_pdf, target_pdf, overwrite = TRUE)
  file.copy(seed_pdf, file.path(figure_dir, "unexpected.pdf"))
  run_check(expected_status = 1L, expected_message = "Unexpected publication figure")
  unlink(file.path(figure_dir, "unexpected.pdf"))
  v4_paths <- file.path(figure_dir, tail(figure_names, 2L))
  unlink(v4_paths)
  run_check(expected_status = 1L, expected_message = "Missing publication figure")
  run_check("--allow-legacy", expected_status = 1L, expected_message = "Unknown or repeated")
  file.copy(seed_pdf, v4_paths[[1L]])
  run_check(expected_status = 1L, expected_message = "Missing publication figure")
  file.copy(seed_pdf, v4_paths[[2L]])
  default_figures <- file.path(result_dir, "figures-and-tables")
  dir.create(default_figures)
  file.copy(file.path(figure_dir, figure_names), default_figures)
  run_check(base_arguments = c("--results-root", result_dir))

  artifact_paths <- c(
    "simulation-results/glm-variance-robustness/illustrative-example-results.rds",
    "simulation-results/nb-permutation-signflip-continuous-dispersion-estimated-theta/nb-permutation-signflip_results.rds",
    "simulation-results/pcm-additive-matched-dcrt/pcm-additive_results.rds",
    "ccle/ccle-typeI-synthY-30targets-matched-v4-results.rds",
    sprintf("frangieh-2021/%s/score_fdr_rep%03d.rds", PERTURB_RESULTS_SUBDIR, 1:20)
  )
  entries <- lapply(seq_along(artifact_paths), function(index) {
    relative_path <- artifact_paths[[index]]
    path <- file.path(result_dir, relative_path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    fixture <- list(fixture = index)
    if (index > 4L) {
      fixture$config <- list(schema_version = PERTURB_SCHEMA_VERSION,
                             q_levels = PERTURB_Q_LEVELS)
    }
    saveRDS(fixture, path)
    list(path = paste0("results/", relative_path), bytes = file.info(path)$size,
         sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
  })
  manifest <- list(schema = "pcit_bundled_results_v1", artifacts = entries)
  manifest_path <- file.path(result_dir, "artifact-manifest.json")
  write_manifest <- function(value) {
    jsonlite::write_json(value, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  }
  write_manifest(manifest)
  run_check("--verify-results", expected_message = "Verified SHA-256 and sizes of 24")
  # Simulate a copied standalone folder after restoring just the two checker
  # dependencies. No parent-project library or R_LIBS_USER is needed there.
  standalone_root <- file.path(temporary_root, "standalone")
  dir.create(standalone_root)
  file.copy(file.path(package_root, c("check_outputs.R", "reproduction-helpers.R", "renv.lock")),
            standalone_root)
  dir.create(file.path(standalone_root, "real-data-analysis"))
  stopifnot(file.copy(file.path(package_root, "real-data-analysis", "method_config.R"),
                     file.path(standalone_root, "real-data-analysis")))
  version <- paste0("R-", R.version$major, ".",
                    strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L])
  standalone_library <- file.path(standalone_root, "renv", "library", "test", version, "test-platform")
  dir.create(standalone_library, recursive = TRUE)
  for (dependency in c("jsonlite", "digest")) {
    stopifnot(file.copy(normalizePath(find.package(dependency)),
                       standalone_library, recursive = TRUE))
  }
  stopifnot(file.copy(result_dir, standalone_root, recursive = TRUE))
  checker <- file.path(standalone_root, "check_outputs.R")
  saved_library_environment <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  saved_override <- Sys.getenv("PCIT_R_LIBRARY", unset = NA_character_)
  saved_results_environment <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = NA_character_)
  saved_figure_environment <- Sys.getenv("PCIT_FIGURE_DIR", unset = NA_character_)
  Sys.unsetenv(c("R_LIBS_USER", "PCIT_R_LIBRARY", "LOCAL_PDRT_DATA_DIR", "PCIT_FIGURE_DIR"))
  run_check("--verify-results", expected_message = "Verified SHA-256 and sizes of 24",
            base_arguments = character())
  checker <- file.path(package_root, "check_outputs.R")
  if (!is.na(saved_library_environment)) Sys.setenv(R_LIBS_USER = saved_library_environment)
  if (!is.na(saved_override)) Sys.setenv(PCIT_R_LIBRARY = saved_override)
  if (!is.na(saved_results_environment)) Sys.setenv(LOCAL_PDRT_DATA_DIR = saved_results_environment)
  if (!is.na(saved_figure_environment)) Sys.setenv(PCIT_FIGURE_DIR = saved_figure_environment)
  # A checksum-valid artifact still has to use the selected single-q schema.
  frangieh_path <- file.path(result_dir, artifact_paths[[5L]])
  original_frangieh <- readRDS(frangieh_path)
  for (field in c("q_levels", "schema_version")) {
    invalid_frangieh <- original_frangieh
    invalid_frangieh$config[[field]] <- if (field == "q_levels") 0.2 else "unexpected_schema"
    saveRDS(invalid_frangieh, frangieh_path)
    bad_manifest <- manifest
    bad_manifest$artifacts[[5L]]$bytes <- file.info(frangieh_path)$size
    bad_manifest$artifacts[[5L]]$sha256 <- digest::digest(
      file = frangieh_path, algo = "sha256", serialize = FALSE
    )
    write_manifest(bad_manifest)
    run_check("--verify-results", expected_status = 1L,
              expected_message = "selected single-q schema")
  }
  saveRDS(original_frangieh, frangieh_path)
  bad_manifest <- manifest
  bad_manifest$artifacts[[1L]]$sha256 <- paste(rep("0", 64L), collapse = "")
  write_manifest(bad_manifest)
  run_check("--verify-results", expected_status = 1L,
            expected_message = "checksum/size mismatch")
  bad_manifest <- manifest
  bad_manifest$artifacts <- bad_manifest$artifacts[-1L]
  write_manifest(bad_manifest)
  run_check("--verify-results", expected_status = 1L,
            expected_message = "exactly the 24 canonical")
  bad_manifest <- manifest
  bad_manifest$artifacts[[1L]]$path <- "results/../outside.rds"
  write_manifest(bad_manifest)
  run_check("--verify-results", expected_status = 1L,
            expected_message = "Malformed bundled-result manifest entry")
  cat("Output-checker regression checks passed.\n")
}
run_tests()
