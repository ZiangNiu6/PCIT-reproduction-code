#!/usr/bin/env bash
# Reproduce only the three simulation studies used by the current manuscript.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SIMULATION_CODE="$SCRIPT_DIR/simulation-code"
FIGURE_CODE="$SIMULATION_CODE/figure-code"

usage() {
  echo "Usage: $0 [full|check|validate-only|plot-only]"
  echo "full computes all three studies in a fresh simulation-results directory, then plots."
  echo "check runs small data-independent method tests; validate-only checks saved results without writes."
  echo "Use a new LOCAL_PDRT_DATA_DIR for full; existing simulation outputs are never reused or replaced."
}

mode="${1:-full}"
if [[ "$mode" == --help || "$mode" == -h ]]; then usage; exit 0; fi
if (($# > 1)) || [[ "$mode" != full && "$mode" != check &&
  "$mode" != validate-only && "$mode" != plot-only ]]; then
  usage >&2
  exit 2
fi

source "$SCRIPT_DIR/workflow-common.sh"
pcit_initialize "$SCRIPT_DIR"
if [[ "$mode" == full ]]; then
  pcit_require_fresh_directory "$LOCAL_PDRT_DATA_DIR/simulation-results"
fi
pcit_check_r
pcit_restore

if [[ "$mode" == full ]] && ! command -v nextflow >/dev/null 2>&1; then
  echo "Nextflow is required for a full simulation run. See reproduction-code/README.md." >&2
  exit 1
fi

if [[ "$mode" == full || "$mode" == check ]]; then
  pcit_r "$SIMULATION_CODE/smoke-test-pcm-additive-matched.R"
  pcit_r "$SIMULATION_CODE/smoke-test-nb-continuous-dispersion.R"
fi
if [[ "$mode" == check ]]; then
  bash "$SIMULATION_CODE/tests/check-workflow.sh"
  echo "Simulation method checks complete. No production study was run."
  exit 0
fi

if [[ "$mode" == full ]]; then
  pcit_r "$SIMULATION_CODE/glm-variance-robustness.R" "--cores=${PCIT_GLM_CORES:-1}"

  "$SIMULATION_CODE/run_simulation.sh" \
    --sim_name nb-permutation-signflip \
    --results_dir "$LOCAL_PDRT_DATA_DIR/simulation-results/nb-permutation-signflip-continuous-dispersion-estimated-theta" \
    --max_gb "${PCIT_SIM_MAX_GB:-4}" \
    --max_hours "${PCIT_SIM_MAX_HOURS:-1}" \
    --profile "${PCIT_NEXTFLOW_PROFILE:-local}"

  "$SIMULATION_CODE/run_simulation.sh" \
    --sim_name pcm-additive \
    --results_dir "$LOCAL_PDRT_DATA_DIR/simulation-results/pcm-additive-matched-dcrt" \
    --max_gb "${PCIT_SIM_MAX_GB:-4}" \
    --max_hours "${PCIT_PCM_SIM_MAX_HOURS:-2}" \
    --profile "${PCIT_NEXTFLOW_PROFILE:-local}"
fi

if [[ "$mode" == validate-only ]]; then
  pcit_r "$FIGURE_CODE/plot-glm-variance-robustness.R" --validate-only
  pcit_r "$FIGURE_CODE/plot-nb-permutation-signflip.R" --validate-only
  pcit_r "$FIGURE_CODE/plot-pcm-additive.R" --validate-only
else
  # Validate every input before generating any figures.
  pcit_r "$FIGURE_CODE/plot-glm-variance-robustness.R" --validate-only
  pcit_r "$FIGURE_CODE/plot-nb-permutation-signflip.R" --validate-only
  pcit_r "$FIGURE_CODE/plot-pcm-additive.R" --validate-only
  pcit_r "$FIGURE_CODE/plot-glm-variance-robustness.R"
  pcit_r "$FIGURE_CODE/plot-nb-permutation-signflip.R"
  pcit_r "$FIGURE_CODE/plot-pcm-additive.R"
fi

echo "Simulation reproduction complete ($mode)."
