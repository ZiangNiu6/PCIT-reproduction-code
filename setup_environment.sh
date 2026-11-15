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
# Rhdf5lib supplies HDF5 archives, but macOS still needs OpenSSL's library
# directory when compiling ondisc. Discover it without personal Makevars.
openssl_libdir=""
if [[ "$(uname -s)" == Darwin ]]; then
  openssl_libdir="${PCIT_OPENSSL_LIB_DIR:-}"
  if [[ -n "$openssl_libdir" && ! -f "$openssl_libdir/libcrypto.dylib" && ! -f "$openssl_libdir/libcrypto.a" ]]; then
    echo "PCIT_OPENSSL_LIB_DIR must contain libcrypto.dylib or libcrypto.a." >&2
    exit 1
  fi
  if [[ -z "$openssl_libdir" ]] && command -v pkg-config >/dev/null 2>&1; then
    for module in libcrypto openssl; do
      candidate=$(pkg-config --variable=libdir "$module" 2>/dev/null || true)
      if [[ -f "$candidate/libcrypto.dylib" || -f "$candidate/libcrypto.a" ]]; then
        openssl_libdir="$candidate"; break
      fi
    done
  fi
  if [[ -z "$openssl_libdir" ]] && command -v brew >/dev/null 2>&1; then
    candidate=$(brew --prefix openssl@3 2>/dev/null || true)
    if [[ -f "$candidate/lib/libcrypto.dylib" || -f "$candidate/lib/libcrypto.a" ]]; then
      openssl_libdir="$candidate/lib"
    fi
  fi
  if [[ -n "$openssl_libdir" ]]; then
    openssl_libdir="$(cd "$openssl_libdir" && pwd -P)"
    export LIBRARY_PATH="$openssl_libdir${LIBRARY_PATH:+:$LIBRARY_PATH}"
  fi
fi
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_PROJECT="$package_root"
# Select the active project as well as the explicit restore target so setup
# works from any caller's working directory.
cd "$package_root"
restore_status=0
Rscript --vanilla - "$package_root" <<'RS' || restore_status=$?
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
if ((restore_status)); then
  if [[ "$(uname -s)" == Darwin && -z "$openssl_libdir" ]]; then
    echo "If ondisc could not link libcrypto, install OpenSSL (e.g. brew install openssl@3 pkg-config) or set PCIT_OPENSSL_LIB_DIR to its lib directory." >&2
  fi
  exit "$restore_status"
fi
