#!/usr/bin/env bash
# Submit one SGE array task per eligible Frangieh perturbation target.
#
# Counts are read from the ondisc ODM in the intermediate directory by
# default.  Set FRANGIEH_RAW_H5AD only to run the scPerturb H5AD backend,
# which reconstructs integer counts from log1p(CPM).
#
# Optional:
#   FRANGIEH_RAW_H5AD=/absolute/path/frangieh_2021_raw.h5ad
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

if [[ -n "${FRANGIEH_RAW_H5AD:-}" && ! -f "$FRANGIEH_RAW_H5AD" ]]; then
  echo "H5AD not found: $FRANGIEH_RAW_H5AD" >&2
  exit 1
fi

# Load the lab data-directory config and export it so child R sessions inherit
# it.  The guard keeps the reproduction package usable for anyone who sets the
# variables by hand.  Without this an unset LOCAL_PDRT_DATA_DIR silently fell
# back to $PROJECT_DIR/results, pointing the whole submission at the wrong tree.
if [[ -z "${LOCAL_PDRT_DATA_DIR:-}" && -f "$HOME/.research_config" ]]; then
  set -a; source "$HOME/.research_config"; set +a
fi
if [[ -z "${LOCAL_PDRT_DATA_DIR:-}" ]]; then
  echo "ERROR: LOCAL_PDRT_DATA_DIR is not set and ~/.research_config is absent." >&2
  echo "Set it, or pass FRANGIEH_INTERMEDIATE_DIR and" >&2
  echo "FRANGIEH_EXPLORATORY_RESULTS_DIR explicitly." >&2
  exit 1
fi

RESULTS_ROOT="$LOCAL_PDRT_DATA_DIR"
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

if [[ -n "${FRANGIEH_RAW_H5AD:-}" ]]; then
  printf -v H5AD_FLAG -- '--h5ad %q ' "$FRANGIEH_RAW_H5AD"
  BACKEND_LABEL="scPerturb H5AD"
  BACKEND_ID="scperturb_h5ad"
else
  H5AD_FLAG=""
  BACKEND_LABEL="ondisc ODM"
  BACKEND_ID="ondisc_odm"
fi

# All array tasks share one locked expression audit and one set of manifest
# CSVs.  Build them here on first use rather than requiring a separate
# --manifest-only pass.  The step hashes the counts backing store once so the
# array tasks inherit the MD5 instead of each rehashing several GB.
if [[ ! -f "$AUDIT_FILE" ]]; then
  echo "No locked expression audit in this output directory."
  echo "Building manifests once ($BACKEND_LABEL); this hashes the counts"
  echo "backing store and takes a few minutes ..."
  # shellcheck disable=SC2086
  if ! "$RSCRIPT_BIN" --vanilla "$RUNNER" $H5AD_FLAG \
        --intermediate-dir "$INTERMEDIATE_DIR" --output-dir "$OUTPUT_DIR" \
        --min-treated "$TREATED_FLOOR" --manifest-only; then
    echo "ERROR: manifest build failed; nothing was submitted." >&2
    exit 1
  fi
fi
if [[ ! -f "$AUDIT_FILE" ]]; then
  echo "ERROR: manifest build did not produce $AUDIT_FILE" >&2
  exit 1
fi

# Refuse a stale audit built for the other backend.  Without this the
# submission would succeed and then all 241 array tasks would fail separately
# on the queue when each revalidated the locked counts-source identity.
# Parsed with read.csv rather than awk: the audit row embeds quoted fields that
# themselves contain commas (the rechecked gene indices), so naive splitting on
# "," reads the wrong column.
AUDIT_BACKEND="$(
  "$RSCRIPT_BIN" --vanilla -e \
    'cat(as.character(utils::read.csv(commandArgs(TRUE)[1L],
       stringsAsFactors = FALSE)$backend[[1L]]))' "$AUDIT_FILE" 2>/dev/null
)"
if [[ "$AUDIT_BACKEND" != "$BACKEND_ID" ]]; then
  echo "ERROR: the existing audit was built for backend '$AUDIT_BACKEND'," >&2
  echo "but this submission requests '$BACKEND_ID'." >&2
  echo "Use a fresh FRANGIEH_EXPLORATORY_RESULTS_DIR for the other backend." >&2
  exit 1
fi
if [[ "$SCAN_MODE_VALUE" == "pilot" ]]; then
  PILOT_FLAG="--pilot "
else
  PILOT_FLAG=""
fi
printf -v JOB_BODY '#!/usr/bin/env bash\nset -euo pipefail\ncd -- %q\nexec %q --vanilla %q %s--intermediate-dir %q --output-dir %q --target-index "$SGE_TASK_ID" %s--min-treated %q --response-block-size %q --permutations %q --seed %q --cores %q\n' \
  "$PROJECT_DIR" "$RSCRIPT_BIN" "$RUNNER" "$H5AD_FLAG" \
  "$INTERMEDIATE_DIR" "$OUTPUT_DIR" "$PILOT_FLAG" "$TREATED_FLOOR" \
  "$BLOCK_SIZE" "$N_PERMUTATIONS" "$BASE_SEED" "$N_CORES"

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
printf 'backend=%s, mode=%s, min-treated=%s, B=%s, seed=%s, response block=%s, cores/task=%s\n' \
  "$BACKEND_LABEL" "$SCAN_MODE_VALUE" "$TREATED_FLOOR" "$N_PERMUTATIONS" \
  "$BASE_SEED" "$BLOCK_SIZE" "$N_CORES"
