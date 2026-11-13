#!/usr/bin/env bash
# Run all 42 new matched-v4 CCLE settings locally or on SGE.
set -euo pipefail
usage() {
  echo "Usage: submit-ccle.sh [--local|--sge]"
  echo "Dispatches all 42 settings and refuses any existing current result files."
}
mode=--sge
if (($# > 1)); then usage >&2; exit 2; fi
if (($#)); then
  case "$1" in
    --local|--sge) mode="$1" ;;
    --help|-h) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
package_root=$(cd -- "$script_dir/.." && pwd)
source "$package_root/workflow-common.sh"
pcit_initialize "$package_root"
runner="$script_dir/ccle-typeI-synthY-30targets.R"
ccle_dir="$LOCAL_PDRT_DATA_DIR/ccle"
for input_file in mutation.txt response.csv; do
  if [[ ! -f "$ccle_dir/$input_file" ]]; then
    echo "Missing CCLE input: $ccle_dir/$input_file" >&2; exit 1
  fi
done
for output in "$ccle_dir"/ccle-typeI-synthY-30targets-matched-v4-results*.rds \
              "$ccle_dir/ccle-typeI-synthY-30targets-matched-v4-rejection-rates.csv"; do
  if [[ -e "$output" ]]; then
    echo "Existing CCLE output: $output. Use a fresh result directory." >&2; exit 1
  fi
done
pcit_check_r
if [[ "$mode" == --sge ]] && ! command -v qsub >/dev/null 2>&1; then
  echo "qsub is required for SGE mode; use --local for sequential execution." >&2; exit 1
fi
if ! mkdir -- "$ccle_dir/.ccle-dispatched"; then
  echo "CCLE jobs were already dispatched here. Use a fresh result directory." >&2; exit 1
fi
if [[ "$mode" == --local ]]; then
  for setup_idx in {1..42}; do pcit_r "$runner" "$setup_idx"; done
else
  rscript_bin=$(command -v Rscript)
  sge_options=(-l "m_mem_free=${MEM:-8G}" -cwd)
  if [[ -n "${PROJECT:-}" ]]; then sge_options+=(-P "$PROJECT"); fi
  for setup_idx in {1..42}; do
    printf -v job_body \
      '#!/usr/bin/env bash\nset -euo pipefail\nexport PATH=%q:"$PATH"\nexport LOCAL_PDRT_DATA_DIR=%q\nexport PCIT_R_LIBRARY=%q\nsource %q\npcit_initialize %q\npcit_check_r\npcit_r %q %q\n' \
      "$(dirname -- "$rscript_bin")" "$LOCAL_PDRT_DATA_DIR" "${PCIT_R_LIBRARY:-}" \
      "$package_root/workflow-common.sh" "$package_root" "$runner" "$setup_idx"
    printf '%s' "$job_body" | qsub "${sge_options[@]}" -N "ccle_matched_v4_${setup_idx}"
  done
fi
echo "Dispatched all 42 CCLE settings. After SGE jobs finish, combine explicitly:"
echo "Rscript --vanilla $script_dir/combine-ccle-results.R"
echo "Then run: $package_root/run_all_semi_synthetic.sh plot-only"
