#!/usr/bin/env bash
# Standalone reproduction of the original Frangieh v2 numerical design.
# Current cohort/guide handling is retained by explicit user instruction.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$SCRIPT_DIR/workflow-common.sh"
usage() {
  cat <<'USAGE'
Usage: ./run_all_real_data.sh [full [n_reps]|plot-only|validate-only|check] [--local|--sge] [--force]
  full           Run missing replicates (default: 20, SGE execution).
  --local        Execute missing replicates sequentially on this machine.
  --sge          Submit missing replicates; use plot-only after completion.
  --force        Intentionally replace selected expected replicate outputs.
  plot-only      Validate all 20, aggregate, and generate figures.
  validate-only  Validate all 20 saved outputs without writing files.
  check          Run the synthetic numerical smoke and temporary-file I/O tests.
Restore packages only with PCIT_RESTORE=1. Lab configuration is opt-in via
PCIT_USE_RESEARCH_CONFIG=1; explicitly supplied paths always win.
USAGE
}
mode=full
n_reps=20
execution=sge
execution_set=0
force=0
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
    --force) force=1 ;;
    --help|-h) usage; exit 0 ;;
    *)
      if [[ "$option" =~ ^[1-9][0-9]*$ ]] && ((number_set == 0)); then
        n_reps="$option"; number_set=1
      else usage >&2; exit 2; fi ;;
  esac
done
if ((n_reps > 20)) || [[ "$mode" != full && ( "$force" == 1 || "$number_set" == 1 || "$execution_set" == 1 ) ]]; then
  usage >&2; exit 2
fi
pcit_initialize "$SCRIPT_DIR"
pcit_check_r
pcit_restore
ANALYSIS_DIR="$SCRIPT_DIR/real-data-analysis"
export FRANGIEH_INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/intermediate-data}"
export FRANGIEH_SCORE_RESULTS_DIR="${FRANGIEH_SCORE_RESULTS_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/score-calibration-nonpositive-nulls}"
export FRANGIEH_SCORE_PLOTS_DIR="${FRANGIEH_SCORE_PLOTS_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/plots}"
for key in FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR FRANGIEH_SCORE_PLOTS_DIR; do
  [[ "${!key}" == /* ]] || export "$key=$SCRIPT_DIR/${!key}"
done
if [[ "$mode" == check ]]; then
  pcit_r "$ANALYSIS_DIR/smoke_test_run.R"
  pcit_r "$ANALYSIS_DIR/test-workflow-io.R"
elif [[ "$mode" == validate-only ]]; then
  pcit_r "$ANALYSIS_DIR/data_postprocessing.R" --validate-only
elif [[ "$mode" == plot-only ]]; then
  pcit_r "$ANALYSIS_DIR/data_postprocessing.R"
  pcit_r "$ANALYSIS_DIR/figure-code/plotting_code.R"
else
  options=()
  if ((force)); then options+=(--force); fi
  missing=$(pcit_r "$ANALYSIS_DIR/validate-results.R" --missing "$n_reps" ${options[@]+"${options[@]}"})
  if [[ -z "$missing" ]]; then
    echo "All selected Frangieh replicate outputs are valid; no analysis needed."
  else
    if [[ -z "${LOCAL_FRANGIEH_2021_DATA_DIR:-}" ]]; then
      echo "Set LOCAL_FRANGIEH_2021_DATA_DIR to the processed Frangieh input root." >&2; exit 1
    fi
    pcit_r "$ANALYSIS_DIR/smoke_test_run.R"
    pcit_r "$ANALYSIS_DIR/data_preprocessing.R" 1000
    bash "$ANALYSIS_DIR/submit_combined_fdr.sh" "$n_reps" "--$execution" ${options[@]+"${options[@]}"}
  fi
  if [[ "$execution" == local && "$n_reps" == 20 ]]; then
    pcit_r "$ANALYSIS_DIR/data_postprocessing.R"
    pcit_r "$ANALYSIS_DIR/figure-code/plotting_code.R"
  elif [[ "$execution" == sge && -n "$missing" ]]; then
    echo "After jobs finish, run: $SCRIPT_DIR/run_all_real_data.sh plot-only"
  fi
fi
