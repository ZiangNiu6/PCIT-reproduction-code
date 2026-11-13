#!/usr/bin/env bash
# Shared setup; sourcing this file never installs packages or starts an analysis.

pcit_initialize() {
  local package_dir="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
  export PCIT_REPRODUCTION_ROOT="$package_dir"
  export PCIT_PROJECT_ROOT="$(cd "$package_dir/.." && pwd -P)"

  # Lab configuration is optional. Explicit variables always win over its values.
  local config_path="${PCIT_CONFIG_FILE:-$HOME/.research_config}"
  if [[ "${PCIT_USE_RESEARCH_CONFIG:-0}" == 1 || -n "${PCIT_CONFIG_FILE:-}" ]]; then
    if [[ ! -f "$config_path" ]]; then
      echo "Configuration file not found: $config_path" >&2
      return 1
    fi
    local names=() values=() key
    for key in LOCAL_PDRT_DATA_DIR LOCAL_FRANGIEH_2021_DATA_DIR \
      FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR FRANGIEH_SCORE_PLOTS_DIR \
      PCIT_FIGURE_DIR PCIT_R_LIBRARY; do
      if [[ -n "${!key+x}" ]]; then names+=("$key"); values+=("${!key}"); fi
    done
    set -a
    source "$config_path"
    set +a
    local i
    for ((i=0; i<${#names[@]}; i++)); do export "${names[$i]}=${values[$i]}"; done
  fi

  export LOCAL_PDRT_DATA_DIR="${LOCAL_PDRT_DATA_DIR:-$package_dir/results}"
  [[ "$LOCAL_PDRT_DATA_DIR" == /* ]] || export LOCAL_PDRT_DATA_DIR="$package_dir/$LOCAL_PDRT_DATA_DIR"
  export PCIT_FIGURE_DIR="${PCIT_FIGURE_DIR:-$LOCAL_PDRT_DATA_DIR/figures-and-tables}"
  [[ "$PCIT_FIGURE_DIR" == /* ]] || export PCIT_FIGURE_DIR="$package_dir/$PCIT_FIGURE_DIR"
  export RENV_CONFIG_SANDBOX_ENABLED=FALSE
  export RENV_PROJECT="$package_dir"
  export R_PROFILE_USER="$package_dir/.Rprofile"
}

pcit_require_fresh_directory() {
  local directory="${1:?A fresh output directory is required}"
  shift
  mkdir -p -- "$directory" || return 1
  local entry name allowed_name allowed
  while IFS= read -r -d '' entry; do
    name="${entry##*/}"
    allowed=0
    for allowed_name in "$@"; do
      if [[ "$name" == "$allowed_name" && -f "$entry" ]]; then allowed=1; break; fi
    done
    if ((allowed == 0)); then
      echo "Full computation requires a fresh output directory: $directory" >&2
      echo "Existing entry: $name. Choose a new output directory; use plot-only for saved results." >&2
      return 1
    fi
  done < <(find "$directory" -mindepth 1 -maxdepth 1 -print0)
  # Leave the claim in place after success or failure. Full runs never resume.
  if ! mkdir -- "$directory/.full-run-started"; then
    echo "Full computation has already started in $directory. Choose a new output directory." >&2
    return 1
  fi
}

pcit_check_r() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is required (the pinned R version is 4.4.1)." >&2
    return 1
  fi
  if [[ -n "${PCIT_R_LIBRARY:-}" ]]; then
    [[ "$PCIT_R_LIBRARY" == /* ]] || export PCIT_R_LIBRARY="$PCIT_REPRODUCTION_ROOT/$PCIT_R_LIBRARY"
    if [[ ! -d "$PCIT_R_LIBRARY" ]]; then echo "R library not found: $PCIT_R_LIBRARY" >&2; return 1; fi
    export R_LIBS_USER="$PCIT_R_LIBRARY${R_LIBS_USER:+:$R_LIBS_USER}"
  else
    local project_libraries
    project_libraries=$(Rscript --vanilla - "$PCIT_REPRODUCTION_ROOT" <<'RS'
root <- commandArgs(TRUE)[[1L]]
version <- paste0("R-", R.version$major, ".", strsplit(R.version$minor, ".", fixed=TRUE)[[1L]][1L])
bases <- c(file.path(root, "renv", "library"), file.path(dirname(root), "renv", "library"))
paths <- unlist(lapply(bases, function(base) Sys.glob(file.path(base, "*", version, "*"))))
paths <- paths[dir.exists(paths)]
cat(paste(unique(paths), collapse = .Platform$path.sep))
RS
)
    if [[ -n "$project_libraries" ]]; then
      export R_LIBS_USER="$project_libraries${R_LIBS_USER:+:$R_LIBS_USER}"
    fi
  fi
}

pcit_r() {
  (cd "$PCIT_REPRODUCTION_ROOT" && Rscript --vanilla "$@")
}

pcit_restore() {
  case "${PCIT_RESTORE:-0}" in
    0) return 0 ;;
    1) "$PCIT_REPRODUCTION_ROOT/setup_environment.sh"; pcit_check_r ;;
    *) echo "PCIT_RESTORE must be 0 or 1." >&2; return 2 ;;
  esac
}
