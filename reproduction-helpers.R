# Paths shared by the standalone reproduction package. No working-directory
# assumptions or manuscript checkout are needed to use this directory alone.
reproduction_root <- function(script_file = NULL) {
  if (is.null(script_file)) {
    script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    script_file <- if (length(script_arg) == 1L) {
      sub("^--file=", "", script_arg)
    } else getwd()
  }
  current <- normalizePath(script_file, mustWork = TRUE)
  if (!dir.exists(current)) current <- dirname(current)
  repeat {
    if (file.exists(file.path(current, "renv.lock")) &&
        file.exists(file.path(current, "check_outputs.R"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the reproduction-code package root.", call. = FALSE)
    }
    current <- parent
  }
}

reproduction_absolute_path <- function(path, root) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, mustWork = FALSE)
}

reproduction_results_root <- function(root = reproduction_root()) {
  path <- Sys.getenv("LOCAL_PDRT_DATA_DIR", unset = "")
  if (!nzchar(path)) path <- file.path(root, "results")
  reproduction_absolute_path(path, root)
}

reproduction_figure_dir <- function(root = reproduction_root()) {
  path <- Sys.getenv("PCIT_FIGURE_DIR", unset = "")
  if (!nzchar(path)) path <- file.path(reproduction_results_root(root), "figures-and-tables")
  reproduction_absolute_path(path, root)
}

# Discover an already restored library without invoking renv's autoloader or
# installing anything. The parent project library is a development fallback;
# a copied standalone package uses its own renv/library after restoration.
reproduction_library_paths <- function(root = reproduction_root()) {
  override <- Sys.getenv("PCIT_R_LIBRARY", unset = "")
  if (nzchar(override)) {
    override <- reproduction_absolute_path(override, root)
    if (!dir.exists(override)) {
      stop("PCIT_R_LIBRARY does not exist: ", override, call. = FALSE)
    }
  } else override <- character()
  version <- paste0("R-", R.version$major, ".",
                    strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L])
  bases <- c(file.path(root, "renv", "library"),
             file.path(dirname(root), "renv", "library"))
  candidates <- unlist(lapply(bases, function(base) {
    Sys.glob(file.path(base, "*", version, "*"))
  }), use.names = FALSE)
  candidates <- candidates[dir.exists(candidates)]
  unique(c(override, candidates))
}
