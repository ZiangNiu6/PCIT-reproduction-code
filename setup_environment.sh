#!/usr/bin/env bash
set -euo pipefail
package_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if (($#)); then
  echo "Usage: $0" >&2
  exit 2
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Install R 4.4.1 and put Rscript on PATH first." >&2
  exit 1
fi
# Verify the exact vendored pcmt source before restoring any packages.
archive="$package_root/renv/cellar/pcmt_b92c1ccb4f8982fb63d35d7b0986cb3fc5af13bd.tar.gz"
expected_sha256="82b6dc0686a376848a44b9baacca9ac3229509902fd5ff277417ba93bbce8d0e"
if [[ ! -f "$archive" ]]; then
  echo "Bundled pcmt source is missing: $archive" >&2
  exit 1
fi
if command -v sha256sum >/dev/null 2>&1; then
  actual_sha256=$(sha256sum "$archive")
elif command -v shasum >/dev/null 2>&1; then
  actual_sha256=$(shasum -a 256 "$archive")
elif command -v openssl >/dev/null 2>&1; then
  actual_sha256=$(openssl dgst -sha256 -r "$archive")
else
  echo "Install sha256sum, shasum, or openssl to verify the bundled source." >&2
  exit 1
fi
actual_sha256=${actual_sha256%% *}
if [[ "$actual_sha256" != "$expected_sha256" ]]; then
  echo "Bundled pcmt source checksum mismatch; restore a clean copy of the archive." >&2
  exit 1
fi
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_PROJECT="$package_root"
export RENV_PATHS_CELLAR="$package_root/renv/cellar"
# renv's cellar lookup uses the active project/cwd, including when restore()
# receives an explicit project argument. Pin both so any caller cwd works.
cd "$package_root"
Rscript --vanilla - "$package_root" <<'RS'
root <- normalizePath(commandArgs(TRUE)[[1L]], mustWork = TRUE)
lock_lines <- readLines(file.path(root, "renv.lock"), warn = FALSE)
record_start <- grep('^[[:space:]]*"renv"[[:space:]]*:[[:space:]]*\\{', lock_lines)
if (length(record_start) != 1L) stop("Cannot identify renv's lockfile record.")
record_lines <- lock_lines[seq.int(record_start, length(lock_lines))]
record_end <- grep('^[[:space:]]*\\}', record_lines)[[1L]]
version_line <- grep('"Version"[[:space:]]*:', record_lines[seq_len(record_end)], value = TRUE)
if (length(version_line) != 1L) stop("Cannot identify renv's pinned version.")
version_match <- regmatches(version_line, regexec('"Version"[[:space:]]*:[[:space:]]*"([0-9.-]+)"', version_line))[[1L]]
if (length(version_match) != 2L) stop("Invalid pinned renv version.")
pinned_version <- version_match[[2L]]
installed_version <- function(lib = NULL) {
  tryCatch(as.character(utils::packageVersion("renv", lib.loc = lib)),
           error = function(e) NA_character_)
}
if (!identical(installed_version(), pinned_version)) {
  bootstrap <- file.path(root, "renv", "bootstrap")
  dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(bootstrap, .libPaths()))
  if (!identical(installed_version(bootstrap), pinned_version)) {
    # The archive URL remains pinned even when CRAN's current release changes.
    urls <- sprintf(c(
      "https://cran.r-project.org/src/contrib/Archive/renv/renv_%s.tar.gz",
      "https://cran.r-project.org/src/contrib/renv_%s.tar.gz"
    ), pinned_version)
    for (url in urls) {
      tryCatch(
        utils::install.packages(url, repos = NULL, type = "source", lib = bootstrap),
        error = function(e) message(conditionMessage(e))
      )
      if (identical(installed_version(bootstrap), pinned_version)) break
    }
    if (!identical(installed_version(bootstrap), pinned_version)) {
      stop("Could not bootstrap renv ", pinned_version,
           " from CRAN into ", bootstrap, ". Check the network and retry.")
    }
  }
}
stopifnot(identical(as.character(utils::packageVersion("renv")), pinned_version))
# Explicitly choose the project library: restore(project=...) alone otherwise
# uses the caller's first library, which can be a system or bootstrap library.
project_library <- renv::paths$library(project = root)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
renv::restore(project = root, library = project_library, prompt = FALSE)
RS
