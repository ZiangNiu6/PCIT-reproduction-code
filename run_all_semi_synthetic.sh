#!/usr/bin/env bash
# Reproduce the 42 matched-v4 CCLE null/representative-effect settings.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ./reproduction-code/run_all_semi_synthetic.sh [full|plot-only|check|validate-only] [--force]

  full       Smoke-test, resume missing jobs among 42 matched-v4 settings,
             strictly combine, and plot (default).
  plot-only  Plot an existing combined artifact without raw data. If no
             combined v4 exists but v4 setup files exist, strictly combine them.
  check      Run small, data-independent calibration/validation checks.
  validate-only  Validate a saved combined result without writing outputs.
  --force    With full, explicitly replace the 42 expected setup outputs.

CCLE_RESULTS_SCHEMA=auto|v4|v3|v2 selects the saved artifact in plot-only/validate.
Auto prefers v4 > v3 > v2 and never falls back across a partial newer run.
Full requires ccle/mutation.txt and ccle/response.csv under the configured
LOCAL_PDRT_DATA_DIR (default: reproduction-code/results).
Set PCIT_RESTORE=1 to restore renv before running; setup is otherwise read-only.
USAGE
}

mode=full
mode_was_set=0
force_rerun=${CCLE_FORCE_RERUN:-0}
if [[ "$force_rerun" != 0 && "$force_rerun" != 1 ]]; then
  echo "ERROR: CCLE_FORCE_RERUN must be exactly 0 or 1." >&2
  exit 2
fi
for option in "$@"; do
  [[ "$option" != validate-only ]] || option=validate
  case "$option" in
    full|plot-only|check|validate)
      if ((mode_was_set)); then usage >&2; exit 2; fi
      mode=$option
      mode_was_set=1
      ;;
    --force) force_rerun=1 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
done
if [[ "$mode" != full && "$force_rerun" == 1 ]]; then
  echo "ERROR: --force/CCLE_FORCE_RERUN=1 applies only to full mode." >&2
  exit 2
fi
schema_request=${CCLE_RESULTS_SCHEMA:-auto}
case "$schema_request" in
  auto|v4|v3|v2) ;;
  *) echo "ERROR: CCLE_RESULTS_SCHEMA must be auto, v4, v3, or v2." >&2; exit 2 ;;
esac
if [[ "$mode" == full && "$schema_request" != auto && "$schema_request" != v4 ]]; then
  echo "ERROR: Full recomputation produces v4; use plot-only for legacy schemas." >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=workflow-common.sh
source "$script_dir/workflow-common.sh"
pcit_initialize "$script_dir"
pcit_check_r
pcit_restore
cd "$PCIT_REPRODUCTION_ROOT"

analysis_dir="$PCIT_REPRODUCTION_ROOT/semi-synthetic-analysis"
combiner="$analysis_dir/combine-ccle-results.R"
submitter="$analysis_dir/submit-ccle.sh"
plotter="$analysis_dir/figure-code/plot-ccle-results.R"

run_checks() {
  pcit_r "$analysis_dir/smoke-test-summary-validation.R"
  pcit_r "$analysis_dir/smoke-test-matched-calibration.R"
}
if [[ "$mode" == check ]]; then
  run_checks
  exit 0
fi

combined_v4="$LOCAL_PDRT_DATA_DIR/ccle/ccle-typeI-synthY-30targets-matched-v4-results.rds"
combined_v3="$LOCAL_PDRT_DATA_DIR/ccle/ccle-typeI-synthY-30targets-matched-v3-results.rds"
combined_v2="$LOCAL_PDRT_DATA_DIR/ccle/ccle-typeI-synthY-30targets-matched-v2-results.rds"
v4_setup_glob="$LOCAL_PDRT_DATA_DIR/ccle/ccle-typeI-synthY-30targets-matched-v4-results-*-delta*-s10.rds"
v3_setup_glob="$LOCAL_PDRT_DATA_DIR/ccle/ccle-typeI-synthY-30targets-matched-v3-results-*-*-s10.rds"

combine_for_plotting() {
  if [[ "$mode" == validate ]]; then
    echo "ERROR: No combined v4 artifact exists; validation does not write outputs." >&2
    echo "Combine the complete 42-file setup grid explicitly or use plot-only." >&2
    exit 1
  fi
  echo "No combined v4 artifact; strictly combining all 42 setup files ..."
  pcit_r "$combiner"
}

if [[ "$mode" == full ]]; then
  run_checks
  N_SETUPS=42 CCLE_FORCE_RERUN="$force_rerun" bash "$submitter" --local
  pcit_r "$combiner"
  schema_request=v4
elif [[ "$schema_request" == auto ]]; then
  if [[ -f "$combined_v4" ]]; then
    schema_request=v4
  elif compgen -G "$v4_setup_glob" >/dev/null; then
    combine_for_plotting
    schema_request=v4
  elif [[ -f "$combined_v3" ]]; then
    schema_request=v3
  elif compgen -G "$v3_setup_glob" >/dev/null; then
    echo "ERROR: Partial v3 outputs exist; refusing to fall back to v2." >&2
    exit 1
  elif [[ -f "$combined_v2" ]]; then
    schema_request=v2
  else
    echo "ERROR: No matched v4/v3/v2 combined result is available." >&2
    exit 1
  fi
elif [[ "$schema_request" == v4 && ! -f "$combined_v4" ]] &&
     compgen -G "$v4_setup_glob" >/dev/null; then
  combine_for_plotting
fi

# Explicit schema selection is preserved. The plotter validates the complete
# selected artifact and rejects stale schema/figure mixtures without fallback.
export CCLE_RESULTS_SCHEMA="$schema_request"
if [[ "$mode" == validate ]]; then
  pcit_r "$plotter" --validate-only
else
  pcit_r "$plotter"
fi
echo "Semi-synthetic reproduction complete ($mode; $schema_request)."
