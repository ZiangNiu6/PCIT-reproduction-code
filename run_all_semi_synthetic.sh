#!/usr/bin/env bash
# Reproduce the 42 matched-v4 CCLE settings in a fresh output directory.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ./run_all_semi_synthetic.sh [full|plot-only|check|validate-only]
  full           Run all 42 settings in a fresh CCLE directory, combine and plot.
  plot-only      Plot the complete current v4 combined artifact without raw data.
  validate-only  Validate that combined artifact without writing outputs.
  check          Run small, data-independent checks.

Full mode allows only mutation.txt, response.csv and input-provenance.json in
<LOCAL_PDRT_DATA_DIR>/ccle before starting. Use a new directory for every run.
USAGE
}
mode=full
if (($# > 1)); then usage >&2; exit 2; fi
if (($#)); then
  case "$1" in
    full|plot-only|check|validate-only) mode="$1" ;;
    --help|-h) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/workflow-common.sh"
pcit_initialize "$script_dir"
if [[ "$mode" == full ]]; then
  ccle_dir="$LOCAL_PDRT_DATA_DIR/ccle"
  for input_file in mutation.txt response.csv; do
    if [[ ! -f "$ccle_dir/$input_file" ]]; then
      echo "Missing CCLE input: $ccle_dir/$input_file" >&2; exit 1
    fi
  done
  pcit_require_fresh_directory "$ccle_dir" mutation.txt response.csv input-provenance.json
fi
pcit_check_r
pcit_restore
analysis_dir="$PCIT_REPRODUCTION_ROOT/semi-synthetic-analysis"
plotter="$analysis_dir/figure-code/plot-ccle-results.R"
run_checks() {
  pcit_r "$analysis_dir/smoke-test-summary-validation.R"
  pcit_r "$analysis_dir/smoke-test-matched-calibration.R"
}
case "$mode" in
  check) run_checks; bash "$analysis_dir/test-fresh-workflow.sh" ;;
  validate-only) pcit_r "$plotter" --validate-only ;;
  plot-only) pcit_r "$plotter" ;;
  full)
    run_checks
    bash "$analysis_dir/submit-ccle.sh" --local
    pcit_r "$analysis_dir/combine-ccle-results.R"
    pcit_r "$plotter"
    ;;
esac
echo "Semi-synthetic reproduction complete ($mode; v4)."
