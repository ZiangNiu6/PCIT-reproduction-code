#!/usr/bin/env bash
# Resume the original Frangieh v2 experiment locally or on SGE.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
source "$PACKAGE_DIR/workflow-common.sh"
n_reps=20
execution=sge
execution_set=0
force=0
number_set=0
for option in "$@"; do
  case "$option" in
    --local|--sge)
      if ((execution_set)); then
        echo "Specify exactly one execution selector: --local or --sge." >&2
        exit 2
      fi
      execution="${option#--}"; execution_set=1 ;;
    --force) force=1 ;;
    *)
      if [[ "$option" =~ ^[1-9][0-9]*$ ]] && ((number_set == 0)); then n_reps="$option"; number_set=1
      else echo "Usage: $0 [n_reps] [--local|--sge] [--force]" >&2; exit 2; fi ;;
  esac
done
if ((n_reps > 20)); then echo "n_reps must be in 1:20" >&2; exit 2; fi
pcit_initialize "$PACKAGE_DIR"
pcit_check_r
export FRANGIEH_INTERMEDIATE_DIR="${FRANGIEH_INTERMEDIATE_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/intermediate-data}"
export FRANGIEH_SCORE_RESULTS_DIR="${FRANGIEH_SCORE_RESULTS_DIR:-$LOCAL_PDRT_DATA_DIR/frangieh-2021/score-calibration-nonpositive-nulls}"
for key in FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR; do
  [[ "${!key}" == /* ]] || export "$key=$PACKAGE_DIR/${!key}"
done
options=()
if ((force)); then options+=(--force); fi
missing=$(pcit_r "$SCRIPT_DIR/validate-results.R" --missing "$n_reps" ${options[@]+"${options[@]}"})
if [[ -z "$missing" ]]; then echo "All selected replicate outputs are valid; nothing to run."; exit 0; fi
pcit_r "$SCRIPT_DIR/validate-results.R" --intermediates
if [[ "$execution" == sge ]] && ! command -v qsub >/dev/null 2>&1; then echo "qsub is required for SGE mode; use --local for sequential execution." >&2; exit 1; fi
rscript_bin=$(command -v Rscript)
for rep in $missing; do
  if [[ "$execution" == local ]]; then
    pcit_r "$SCRIPT_DIR/run_combined_fdr.R" "$rep" ${options[@]+"${options[@]}"}
  else
    printf -v padded '%03d' "$rep"
    qsub_args=(-l "m_mem_free=${MEM:-16G}" -cwd -N "score5v2_rep$padded")
    if [[ -n "${PROJECT:-}" ]]; then qsub_args+=(-P "$PROJECT"); fi
    {
      printf '#!/usr/bin/env bash\nset -euo pipefail\n'
      for key in LOCAL_PDRT_DATA_DIR FRANGIEH_INTERMEDIATE_DIR FRANGIEH_SCORE_RESULTS_DIR LOCAL_FRANGIEH_2021_DATA_DIR R_LIBS_USER RENV_CONFIG_SANDBOX_ENABLED; do
        if [[ -n "${!key+x}" ]]; then printf 'export %s=%q\n' "$key" "${!key}"; fi
      done
      printf 'cd -- %q\nexec %q --vanilla %q %q' "$PACKAGE_DIR" "$rscript_bin" "$SCRIPT_DIR/run_combined_fdr.R" "$rep"
      if ((force)); then printf ' --force'; fi
      printf '\n'
    } | qsub "${qsub_args[@]}"
  fi
done
