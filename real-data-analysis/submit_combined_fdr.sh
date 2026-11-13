#!/usr/bin/env bash
# Dispatch every selected Frangieh q=0.10 replicate once, locally or on SGE.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
source "$PACKAGE_DIR/workflow-common.sh"
n_reps=20
execution=sge
execution_set=0
number_set=0
for option in "$@"; do
  case "$option" in
    --local|--sge)
      if ((execution_set)); then
        echo "Specify exactly one execution selector: --local or --sge." >&2
        exit 2
      fi
      execution="${option#--}"; execution_set=1 ;;
    *)
      if [[ "$option" =~ ^[1-9][0-9]*$ ]] && ((number_set == 0)); then n_reps="$option"; number_set=1
      else echo "Usage: $0 [n_reps] [--local|--sge]" >&2; exit 2; fi ;;
  esac
done
if ((n_reps > 20)); then echo "n_reps must be in 1:20" >&2; exit 2; fi
pcit_initialize "$PACKAGE_DIR"
export FRANGIEH_INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/intermediate-data-gb2024-single-guide}"
export FRANGIEH_SCORE_RESULTS_DIR="${FRANGIEH_SCORE_RESULTS_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/score-calibration-gb2024-q10}"
for key in FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR; do
  [[ "${!key}" == /* ]] || export "$key=$PACKAGE_DIR/${!key}"
done
# Dispatch is a single fresh batch; an existing output is never skipped.
if [[ -e "$FRANGIEH_SCORE_RESULTS_DIR" ]] && [[ ! -d "$FRANGIEH_SCORE_RESULTS_DIR" || -n "$(ls -A "$FRANGIEH_SCORE_RESULTS_DIR")" ]]; then
  echo "Frangieh replicate outputs require a fresh empty directory: $FRANGIEH_SCORE_RESULTS_DIR" >&2
  exit 1
fi
mkdir -p "$FRANGIEH_SCORE_RESULTS_DIR"
if ! mkdir "$FRANGIEH_SCORE_RESULTS_DIR/.replicates-dispatched" 2>/dev/null; then
  echo "Frangieh replicates have already been dispatched; choose a fresh result directory." >&2
  exit 1
fi
pcit_check_r
pcit_r "$SCRIPT_DIR/validate-results.R" --intermediates
if [[ "$execution" == sge ]] && ! command -v qsub >/dev/null 2>&1; then echo "qsub is required for SGE mode; use --local for sequential execution." >&2; exit 1; fi
rscript_bin=$(command -v Rscript)
for ((rep=1; rep<=n_reps; rep++)); do
  if [[ "$execution" == local ]]; then
    pcit_r "$SCRIPT_DIR/run_combined_fdr.R" "$rep"
  else
    printf -v padded '%03d' "$rep"
    qsub_args=(-l "m_mem_free=${MEM:-16G}" -cwd -N "score5q10_rep$padded")
    if [[ -n "${PROJECT:-}" ]]; then qsub_args+=(-P "$PROJECT"); fi
    {
      printf '#!/usr/bin/env bash\nset -euo pipefail\n'
      for key in LOCAL_PDRT_DATA_DIR FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR LOCAL_FRANGIEH_2021_DATA_DIR R_LIBS_USER RENV_CONFIG_SANDBOX_ENABLED; do
        if [[ -n "${!key+x}" ]]; then printf 'export %s=%q\n' "$key" "${!key}"; fi
      done
      printf 'cd -- %q\nexec %q --vanilla %q %q\n' "$PACKAGE_DIR" "$rscript_bin" "$SCRIPT_DIR/run_combined_fdr.R" "$rep"
    } | qsub "${qsub_args[@]}"
  fi
done
