#!/usr/bin/env bash
# Reproduce only the three simulation studies used by the current manuscript.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SIMULATION_CODE="$SCRIPT_DIR/simulation-code"
FIGURE_CODE="$SIMULATION_CODE/figure-code"

usage() {
  echo "Usage: $0 [full|check|validate-only|plot-only]"
  echo "full validates and reuses complete results, computes missing studies, then plots."
  echo "check runs small data-independent method tests; validate-only checks saved results without writes."
  echo "Use an empty LOCAL_PDRT_DATA_DIR for a fresh full rerun. Invalid existing files fail closed."
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
pcit_check_r
pcit_restore

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
  glm_results="$LOCAL_PDRT_DATA_DIR/simulation-results/glm-variance-robustness/illustrative-example-results.rds"
  if [[ -f "$glm_results" ]]; then
    if ! pcit_r "$FIGURE_CODE/plot-glm-variance-robustness.R" --validate-only; then
      echo "Invalid existing GLM result: $glm_results" >&2
      echo "Review and move that file out of this result namespace, or use an empty LOCAL_PDRT_DATA_DIR." >&2
      exit 1
    fi
    echo "Reusing validated GLM result: $glm_results"
  else
    pcit_r "$SIMULATION_CODE/glm-variance-robustness.R" "--cores=${PCIT_GLM_CORES:-1}"
  fi

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
