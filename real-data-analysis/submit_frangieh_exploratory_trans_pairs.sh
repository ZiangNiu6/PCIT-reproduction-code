#!/usr/bin/env bash
# Submit one SGE array task per eligible Frangieh perturbation target.
#
# Required:
#   FRANGIEH_RAW_H5AD=/absolute/path/frangieh_2021_raw.h5ad
#
# Optional:
#   FRANGIEH_INTERMEDIATE_DIR=/absolute/path/intermediate-data
#   FRANGIEH_EXPLORATORY_RESULTS_DIR=/absolute/path/output
#   LOCAL_PDRT_DATA_DIR=/absolute/path/results-root
#   PROJECT=stat_bhaswar_team  MEM=16G  CORES=1  SGE_PE=shared
#   PERMUTATIONS=4999  SEED=20260822  RESPONSE_BLOCK_SIZE=100
#   MIN_TREATED=10  SCAN_MODE=pilot   # SCAN_MODE is pilot or full
#
# CORES defaults to one.  For CORES>1, SGE_PE must name the site's parallel
# environment so the requested cores and the runner's --cores agree.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

if [[ -z "${FRANGIEH_RAW_H5AD:-}" ]]; then
  echo "Set FRANGIEH_RAW_H5AD to the raw scPerturb Frangieh H5AD." >&2
  exit 1
fi
if [[ ! -f "$FRANGIEH_RAW_H5AD" ]]; then
  echo "H5AD not found: $FRANGIEH_RAW_H5AD" >&2
  exit 1
fi

RESULTS_ROOT="${LOCAL_PDRT_DATA_DIR:-$PROJECT_DIR/results}"
INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$RESULTS_ROOT/frangieh-2021/intermediate-data}"
OUTPUT_DIR="${FRANGIEH_EXPLORATORY_RESULTS_DIR:-$RESULTS_ROOT/frangieh-2021/exploratory-trans-pairs}"
SGE_PROJECT="${PROJECT:-stat_bhaswar_team}"
MEMORY="${MEM:-16G}"
N_CORES="${CORES:-1}"
N_PERMUTATIONS="${PERMUTATIONS:-4999}"
BASE_SEED="${SEED:-20260822}"
BLOCK_SIZE="${RESPONSE_BLOCK_SIZE:-100}"
TREATED_FLOOR="${MIN_TREATED:-10}"
SCAN_MODE_VALUE="${SCAN_MODE:-pilot}"

for value_name in N_CORES N_PERMUTATIONS BASE_SEED BLOCK_SIZE TREATED_FLOOR; do
  value="${!value_name}"
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "$value_name must be a positive integer" >&2
    exit 1
  fi
done
if [[ "$TREATED_FLOOR" != "10" && "$TREATED_FLOOR" != "5" ]]; then
  echo "MIN_TREATED must be 10 (primary) or 5 (sensitivity)" >&2
  exit 1
fi
if [[ "$SCAN_MODE_VALUE" != "pilot" && "$SCAN_MODE_VALUE" != "full" ]]; then
  echo "SCAN_MODE must be pilot or full" >&2
  exit 1
fi
if [[ "$SCAN_MODE_VALUE" == "pilot" && "$TREATED_FLOOR" != "10" ]]; then
  echo "The frozen pilot requires MIN_TREATED=10" >&2
  exit 1
fi
if ((N_PERMUTATIONS < 19)); then
  echo "PERMUTATIONS must be at least 19" >&2
  exit 1
fi
if ((N_CORES > 1)) && [[ -z "${SGE_PE:-}" ]]; then
  echo "Set SGE_PE when requesting CORES>1." >&2
  exit 1
fi
if [[ ! -d "$INTERMEDIATE_DIR" ]]; then
  echo "Intermediate directory not found: $INTERMEDIATE_DIR" >&2
  exit 1
fi
if ! command -v qsub >/dev/null 2>&1; then
  echo "qsub is required" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR/logs"
RSCRIPT_BIN="$(command -v Rscript)"
RUNNER="$SCRIPT_DIR/run_frangieh_exploratory_trans_pairs.R"
AUDIT_FILE="$OUTPUT_DIR/manifests/frangieh-exploratory-expression-audit.csv"
if [[ ! -f "$AUDIT_FILE" ]]; then
  echo "Missing locked expression audit: $AUDIT_FILE" >&2
  echo "Run the runner once with --manifest-only and this output directory." >&2
  exit 1
fi

if [[ "$SCAN_MODE_VALUE" == "pilot" ]]; then
  printf -v JOB_BODY '#!/usr/bin/env bash\nset -euo pipefail\ncd -- %q\nexec %q --vanilla %q --h5ad %q --intermediate-dir %q --output-dir %q --target-index "$SGE_TASK_ID" --pilot --min-treated %q --response-block-size %q --permutations %q --seed %q --cores %q\n' \
    "$PROJECT_DIR" "$RSCRIPT_BIN" "$RUNNER" "$FRANGIEH_RAW_H5AD" \
    "$INTERMEDIATE_DIR" "$OUTPUT_DIR" "$TREATED_FLOOR" "$BLOCK_SIZE" \
    "$N_PERMUTATIONS" "$BASE_SEED" "$N_CORES"
else
  printf -v JOB_BODY '#!/usr/bin/env bash\nset -euo pipefail\ncd -- %q\nexec %q --vanilla %q --h5ad %q --intermediate-dir %q --output-dir %q --target-index "$SGE_TASK_ID" --min-treated %q --response-block-size %q --permutations %q --seed %q --cores %q\n' \
    "$PROJECT_DIR" "$RSCRIPT_BIN" "$RUNNER" "$FRANGIEH_RAW_H5AD" \
    "$INTERMEDIATE_DIR" "$OUTPUT_DIR" "$TREATED_FLOOR" "$BLOCK_SIZE" \
    "$N_PERMUTATIONS" "$BASE_SEED" "$N_CORES"
fi

N_TASKS=241
if [[ "$TREATED_FLOOR" == "5" ]]; then
  N_TASKS=244
fi

QSUB_ARGS=(
  -P "$SGE_PROJECT"
  -l "m_mem_free=$MEMORY"
  -cwd
  -t "1-$N_TASKS"
  -N frangieh_trans
  -o "$OUTPUT_DIR/logs"
  -e "$OUTPUT_DIR/logs"
)
if ((N_CORES > 1)); then
  QSUB_ARGS+=(-pe "$SGE_PE" "$N_CORES")
fi

printf '%s' "$JOB_BODY" | qsub "${QSUB_ARGS[@]}"

if [[ "$SCAN_MODE_VALUE" == "pilot" ]]; then
  printf 'Submitted frozen 5,000-pair pilot across 241 balanced tasks\n'
else
  printf 'Submitted Frangieh full exploratory grid: %s tasks x 5,000 responses\n' \
    "$N_TASKS"
fi
printf 'Output: %s\n' "$OUTPUT_DIR"
printf 'mode=%s, min-treated=%s, B=%s, seed=%s, response block=%s, cores/task=%s\n' \
  "$SCAN_MODE_VALUE" "$TREATED_FLOOR" "$N_PERMUTATIONS" "$BASE_SEED" \
  "$BLOCK_SIZE" "$N_CORES"
