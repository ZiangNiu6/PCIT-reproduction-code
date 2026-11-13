#!/usr/bin/env bash

# Run one fresh formal simulatr specification through the Nextflow pipeline.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd -P)"

sim_spec_dir="$SCRIPT_DIR/sim-spec-formal"
results_dir=""
sim_name=""
profile="local"
B_check=3
B=0
max_gb=4
max_hours=1
benchmark_memory=0
pipeline_revision="${PCIT_SIMULATR_PIPELINE_REVISION:-dd2f1fcdb533a7764f05dd9dd05b175b02b264a8}"

usage() {
  cat <<'EOF'
Usage: run_simulation.sh --sim_name NAME --results_dir DIR [options]

Options:
  --results_dir DIR      New or empty study directory; existing work is refused
  --sim_spec_dir DIR     Directory containing NAME.R (default: sim-spec-formal)
  --profile PROFILE      Nextflow profile: local or standard (default: local)
  --B_check INTEGER      Benchmark replications (default: 3)
  --B INTEGER            Main replications; 0 uses the specifier default
  --max_gb NUMBER        Requested memory in GB (default: 4)
  --max_hours NUMBER     Requested wall time in hours (default: 1)
  --benchmark_memory 0|1 Enable simulatr memory benchmarking (default: 0)
EOF
}

while (($# > 0)); do
  case "$1" in
    --sim_name|--sim-name)
      sim_name="${2:?Missing value for $1}"
      shift 2
      ;;
    --sim_spec_dir|--sim-spec-dir)
      sim_spec_dir="${2:?Missing value for $1}"
      shift 2
      ;;
    --results_dir|--results-dir)
      results_dir="${2:?Missing value for $1}"
      shift 2
      ;;
    --profile)
      profile="${2:?Missing value for $1}"
      shift 2
      ;;
    --B_check|--B-check)
      B_check="${2:?Missing value for $1}"
      shift 2
      ;;
    --B)
      B="${2:?Missing value for $1}"
      shift 2
      ;;
    --max_gb|--max-gb)
      max_gb="${2:?Missing value for $1}"
      shift 2
      ;;
    --max_hours|--max-hours)
      max_hours="${2:?Missing value for $1}"
      shift 2
      ;;
    --benchmark_memory|--benchmark-memory)
      benchmark_memory="${2:?Missing value for $1}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$sim_name" || -z "$results_dir" ]]; then
  echo "Both --sim_name and --results_dir are required." >&2
  usage >&2
  exit 2
fi
if [[ "$profile" != "local" && "$profile" != "standard" ]]; then
  echo "--profile must be local or standard." >&2
  exit 2
fi
if [[ "$sim_name" != nb-permutation-signflip && "$sim_name" != pcm-additive ]]; then
  echo "Supported manuscript studies: nb-permutation-signflip and pcm-additive." >&2
  exit 2
fi
if [[ ! "$B_check" =~ ^[1-9][0-9]*$ || ! "$B" =~ ^(0|[1-9][0-9]*)$ ||
      ! "$max_gb" =~ ^[0-9]+([.][0-9]+)?$ ||
      ! "$max_hours" =~ ^[0-9]+([.][0-9]+)?$ ||
      ! "$benchmark_memory" =~ ^[01]$ ]] ||
    ! awk -v gb="$max_gb" -v hours="$max_hours" 'BEGIN { exit !(gb > 0 && hours > 0) }'; then
  echo "Invalid replication count, memory, wall time, or benchmark-memory option." >&2
  exit 2
fi
if [[ ! -d "$sim_spec_dir" ]]; then
  echo "Simulation-specifier directory not found: $sim_spec_dir" >&2
  exit 2
fi

sim_spec_dir="$(cd "$sim_spec_dir" && pwd -P)"
sim_spec_script="$sim_spec_dir/$sim_name.R"
grid_script="$sim_spec_dir/${sim_name}-grid.R"
if [[ ! -f "$sim_spec_script" ]]; then
  echo "Simulation specifier not found: $sim_spec_script" >&2
  exit 2
fi
if [[ ! -f "$grid_script" ]]; then
  echo "Parameter-grid source not found: $grid_script" >&2
  exit 2
fi

source "$PACKAGE_ROOT/workflow-common.sh"
pcit_initialize "$PACKAGE_ROOT"
pcit_require_fresh_directory "$results_dir"
pcit_check_r

results_dir="$(cd "$results_dir" && pwd -P)"
output_file="$results_dir/${sim_name}_results.rds"
specifier_file="$results_dir/${sim_name}_specifier.rds"

validate_result() {
  if [[ "$sim_name" == nb-permutation-signflip ]]; then
    PCIT_NB_RESULTS_FILE="$output_file" \
      pcit_r "$SCRIPT_DIR/figure-code/plot-nb-permutation-signflip.R" \
        "$output_file" "$grid_script" --validate-only
  else
    PCIT_PCM_RESULTS_FILE="$output_file" \
      pcit_r "$SCRIPT_DIR/figure-code/plot-pcm-additive.R" --validate-only
  fi
}

if ! command -v nextflow >/dev/null 2>&1; then
  echo "Nextflow is required for a fresh simulation run. See reproduction-code/README.md." >&2
  exit 1
fi

export PCIT_SIM_SPEC_SCRIPT="$sim_spec_script"
export PCIT_SIM_SPECIFIER_FILE="$specifier_file"
echo "Building simulatr specifier from current source: $specifier_file"
pcit_r -e '
  suppressPackageStartupMessages(library(simulatr))
  source(Sys.getenv("PCIT_SIM_SPEC_SCRIPT"))
  if (!exists("simulatr_spec", inherits = FALSE)) {
    stop("The specifier source did not create simulatr_spec.")
  }
  saveRDS(simulatr_spec, Sys.getenv("PCIT_SIM_SPECIFIER_FILE"))
'

printf 'env.R_LIBS_USER = "%s"\n' "$R_LIBS_USER" > "$results_dir/nextflow.config"

echo "Running simulation: $sim_name"
(
  cd "$results_dir"
  export NXF_WORK="$results_dir/.nextflow-work"
  nextflow pull katsevich-lab/simulatr-pipeline -r "$pipeline_revision"
  nextflow run katsevich-lab/simulatr-pipeline -r "$pipeline_revision" \
    --simulatr_specifier_fp "$specifier_file" \
    --result_dir "$results_dir" \
    --result_file_name "${sim_name}_results.rds" \
    --B_check "$B_check" \
    --B "$B" \
    --max_gb "$max_gb" \
    --max_hours "$max_hours" \
    --benchmark_memory "$benchmark_memory" \
    -profile "$profile" \
    -with-trace "$results_dir/${sim_name}_trace.txt"
)

if [[ ! -f "$output_file" ]]; then
  echo "Simulation finished without producing $output_file" >&2
  exit 1
fi

validate_result
export PCIT_SIMULATR_PIPELINE_REVISION_USED="$pipeline_revision"
export PCIT_SIM_RESULT_DIR="$results_dir"
pcit_r - <<'RS'
result_dir <- Sys.getenv("PCIT_SIM_RESULT_DIR")
sources <- c(Sys.getenv("PCIT_SIM_SPEC_SCRIPT"),
             sub("[.]R$", "-grid.R", Sys.getenv("PCIT_SIM_SPEC_SCRIPT")))
sources <- sources[file.exists(sources)]
provenance <- list(
  pipeline_revision = Sys.getenv("PCIT_SIMULATR_PIPELINE_REVISION_USED"),
  source_md5 = tools::md5sum(sources),
  lockfile_md5 = tools::md5sum(file.path(Sys.getenv("PCIT_REPRODUCTION_ROOT"), "renv.lock")),
  session_info = utils::sessionInfo()
)
temporary <- tempfile(".provenance-", tmpdir = result_dir)
saveRDS(provenance, temporary)
if (!file.rename(temporary, file.path(result_dir, "reproduction-provenance.rds"))) {
  unlink(temporary)
  stop("Could not publish simulation provenance.")
}
RS

echo "Saved simulation result: $output_file"
