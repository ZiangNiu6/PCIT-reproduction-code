#!/bin/bash
# Submit independent five-method score-calibration replicates on HPC3.
#
# Usage: bash submit_combined_fdr.sh [n_reps]  # default: 20
#
# Each job runs the five Frangieh methods on one 50-positive/950-null
# replicate. Outputs are written to:
#   results/score-calibration-nonpositive-nulls/score_fdr_repNNN.rds

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
N_REPS="${1:-20}"

# Load path configuration for direct invocation. Required values are embedded
# explicitly in each job script; the full login environment is not forwarded.
if [[ -f "$HOME/.research_config" ]]; then
  # shellcheck disable=SC1090
  source "$HOME/.research_config"
fi

SGE_PROJECT="${PROJECT:-stat_bhaswar_team}"
MEMORY="${MEM:-16G}"

RESULTS_ROOT="${LOCAL_PDRT_DATA_DIR:-$PROJECT_DIR/results}"
if [[ "$RESULTS_ROOT" != /* ]]; then
  RESULTS_ROOT="$PROJECT_DIR/$RESULTS_ROOT"
fi
if [[ "$RESULTS_ROOT" != "/" ]]; then
  RESULTS_ROOT=${RESULTS_ROOT%/}
fi
INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$RESULTS_ROOT/frangieh-2021/intermediate-data}"
SCORE_RESULTS_DIR="${FRANGIEH_SCORE_RESULTS_DIR:-$RESULTS_ROOT/frangieh-2021/score-calibration-nonpositive-nulls}"
if [[ "$INTERMEDIATE_DIR" != /* ]]; then
  INTERMEDIATE_DIR="$PROJECT_DIR/$INTERMEDIATE_DIR"
fi
if [[ "$SCORE_RESULTS_DIR" != /* ]]; then
  SCORE_RESULTS_DIR="$PROJECT_DIR/$SCORE_RESULTS_DIR"
fi
if [[ "$INTERMEDIATE_DIR" != "/" ]]; then
  INTERMEDIATE_DIR=${INTERMEDIATE_DIR%/}
fi
if [[ "$SCORE_RESULTS_DIR" != "/" ]]; then
  SCORE_RESULTS_DIR=${SCORE_RESULTS_DIR%/}
fi

if [[ ! "$N_REPS" =~ ^[1-9][0-9]*$ ]] || ((N_REPS > 20)); then
  echo "n_reps must be an integer between 1 and 20" >&2
  exit 1
fi
if ! command -v qsub >/dev/null 2>&1; then
  echo "qsub is required to submit the Frangieh replicate jobs" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required to submit the Frangieh replicate jobs" >&2
  exit 1
fi
if [[ ! -d "$INTERMEDIATE_DIR" ]]; then
  echo "Frangieh intermediate directory not found: $INTERMEDIATE_DIR" >&2
  echo "Run bash run_all_real_data.sh full before direct submission." >&2
  exit 1
fi
rscript_bin=$(command -v Rscript)

printf 'Project root:                %s\n' "$PROJECT_DIR"
printf 'Rscript:                     %s\n' "$rscript_bin"
printf 'LOCAL_PDRT_DATA_DIR:         %s\n' "$RESULTS_ROOT"
printf 'FRANGIEH_INTERMEDIATE_DIR:   %s\n' "$INTERMEDIATE_DIR"
printf 'FRANGIEH_SCORE_RESULTS_DIR:  %s\n' "$SCORE_RESULTS_DIR"

for ((i = 1; i <= N_REPS; i++)); do
  printf -v padded "%03d" "$i"
  printf '#!/usr/bin/env bash\nset -euo pipefail\nexport LOCAL_PDRT_DATA_DIR=%q\nexport FRANGIEH_INTERMEDIATE_DIR=%q\nexport FRANGIEH_SCORE_RESULTS_DIR=%q\nexport RENV_CONFIG_SANDBOX_ENABLED=FALSE\ncd -- %q\nexec %q %q %d\n' \
    "$RESULTS_ROOT" "$INTERMEDIATE_DIR" "$SCORE_RESULTS_DIR" \
    "$PROJECT_DIR" "$rscript_bin" "$SCRIPT_DIR/run_combined_fdr.R" "$i" | \
    qsub \
      -P "$SGE_PROJECT" \
      -l "m_mem_free=$MEMORY" \
      -cwd \
      -N "score5v2_rep${padded}"
done

printf 'Submitted %d five-method jobs (score5v2_rep001 .. score5v2_rep%03d)\n' \
  "$N_REPS" "$N_REPS"
