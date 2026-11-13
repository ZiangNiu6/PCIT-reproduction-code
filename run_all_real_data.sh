#!/usr/bin/env bash
# Frangieh score comparison with GB2024 preprocessing at BH q=0.10.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/workflow-common.sh"
usage() {
  cat <<'USAGE'
Usage: ./run_all_real_data.sh [full [n_reps]|plot-only|validate-only|check] [--local|--sge]
  full           Compute q=0.10 replicates in a fresh Frangieh directory (default: 20).
  --local        Execute all selected replicates sequentially on this machine.
  --sge          Submit all selected replicates; use plot-only after completion.
  plot-only      Validate all 20, aggregate, and generate figures.
  validate-only  Validate all 20 saved outputs without writing files.
  check          Run the small numerical, preprocessing and result-contract checks.
All modes use the single-guide GB2024 cohort and BH q=0.10.
Every full invocation requires an empty <results-root>/frangieh-2021 directory;
existing or failed runs are never resumed or overwritten. Choose a new result
root for another run. Restore packages only with PCIT_RESTORE=1.
Lab configuration is opt-in via PCIT_USE_RESEARCH_CONFIG=1.
USAGE
}
mode=full
n_reps=20
execution=sge
execution_set=0
mode_set=0
number_set=0
for option in "$@"; do
  case "$option" in
    full|plot-only|validate-only|check)
      if ((mode_set)); then usage >&2; exit 2; fi
      mode="$option"; mode_set=1 ;;
    --local|--sge)
      if ((execution_set)); then
        echo "Specify exactly one execution selector: --local or --sge." >&2
        exit 2
      fi
      execution="${option#--}"; execution_set=1 ;;
    --help|-h) usage; exit 0 ;;
    *)
      if [[ "$option" =~ ^[1-9][0-9]*$ ]] && ((number_set == 0)); then
        n_reps="$option"; number_set=1
      else usage >&2; exit 2; fi ;;
  esac
done
if ((n_reps > 20)) || [[ "$mode" != full && ( "$number_set" == 1 || "$execution_set" == 1 ) ]]; then
  usage >&2; exit 2
fi
pcit_initialize "$SCRIPT_DIR"
ANALYSIS_DIR="$SCRIPT_DIR/real-data-analysis"
domain_dir="$LOCAL_PDRT_DATA_DIR/frangieh-2021"
export FRANGIEH_INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$domain_dir/intermediate-data-gb2024-single-guide}"
export FRANGIEH_SCORE_RESULTS_DIR="${FRANGIEH_SCORE_RESULTS_DIR:-$domain_dir/score-calibration-gb2024-q10}"
export FRANGIEH_SCORE_PLOTS_DIR="${FRANGIEH_SCORE_PLOTS_DIR:-$domain_dir/plots}"
for key in FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR FRANGIEH_SCORE_PLOTS_DIR; do
  [[ "${!key}" == /* ]] || export "$key=$SCRIPT_DIR/${!key}"
done
if [[ "$mode" == full ]]; then
  for key in FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR FRANGIEH_SCORE_PLOTS_DIR; do
    case "${!key}" in
      "$domain_dir"/*) relative_path="${!key#"$domain_dir"/}" ;;
      *) echo "$key must be inside the fresh Frangieh directory: $domain_dir" >&2; exit 1 ;;
    esac
    case "/$relative_path/" in
      */../*) echo "$key must not traverse outside the fresh Frangieh directory." >&2; exit 1 ;;
    esac
  done
  pcit_require_fresh_directory "$domain_dir"
fi
pcit_check_r
pcit_restore
if [[ "$mode" == check ]]; then
  pcit_r "$ANALYSIS_DIR/smoke_test_run.R"
  pcit_r "$ANALYSIS_DIR/test-workflow-io.R"
  pcit_r "$ANALYSIS_DIR/test-gb2024-preprocessing.R"
  pcit_r "$ANALYSIS_DIR/test-result-contracts.R"
elif [[ "$mode" == validate-only ]]; then
  pcit_r "$ANALYSIS_DIR/data_postprocessing.R" --validate-only
elif [[ "$mode" == plot-only ]]; then
  pcit_r "$ANALYSIS_DIR/data_postprocessing.R"
  pcit_r "$ANALYSIS_DIR/figure-code/plotting_code.R"
else
  if [[ -z "${LOCAL_FRANGIEH_2021_DATA_DIR:-}" ]]; then
    echo "Set LOCAL_FRANGIEH_2021_DATA_DIR to the processed Frangieh input root." >&2; exit 1
  fi
  pcit_r "$ANALYSIS_DIR/smoke_test_run.R"
  pcit_r "$ANALYSIS_DIR/data_preprocessing.R" 1000
  bash "$ANALYSIS_DIR/submit_combined_fdr.sh" "$n_reps" "--$execution"
  if [[ "$execution" == local && "$n_reps" == 20 ]]; then
    pcit_r "$ANALYSIS_DIR/data_postprocessing.R"
    pcit_r "$ANALYSIS_DIR/figure-code/plotting_code.R"
  elif [[ "$execution" == sge ]]; then
    echo "After jobs finish, run: $SCRIPT_DIR/run_all_real_data.sh plot-only"
  fi
fi
