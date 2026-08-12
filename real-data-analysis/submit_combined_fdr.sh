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
SGE_PROJECT="${PROJECT:-stat_bhaswar_team}"
MEMORY="${MEM:-16G}"

if [[ ! "$N_REPS" =~ ^[1-9][0-9]*$ ]] || ((N_REPS > 20)); then
  echo "n_reps must be an integer between 1 and 20" >&2
  exit 1
fi
if ! command -v qsub >/dev/null 2>&1; then
  echo "qsub is required to submit the Frangieh replicate jobs" >&2
  exit 1
fi

for ((i = 1; i <= N_REPS; i++)); do
  printf -v padded "%03d" "$i"
  printf '#!/usr/bin/env bash\nset -euo pipefail\ncd %q\nRscript %q %d\n' \
    "$PROJECT_DIR" "$SCRIPT_DIR/run_combined_fdr.R" "$i" | \
    qsub \
      -V \
      -P "$SGE_PROJECT" \
      -l "m_mem_free=$MEMORY" \
      -cwd \
      -N "score5v2_rep${padded}"
done

printf 'Submitted %d five-method jobs (score5v2_rep001 .. score5v2_rep%03d)\n' \
  "$N_REPS" "$N_REPS"
