#!/usr/bin/env bash
# Run or submit the 60 CCLE semi-synthetic target/sparsity setups.
#
# Usage:
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --local
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --sge
#
# --local runs each setup synchronously. It is used by the repository-level
# full runner so combining and plotting begin only after computation succeeds.
# --sge submits independent Grid Engine jobs and returns after submission; run
# the repository-level plot-only workflow after those jobs have completed.
#
# Optional environment variables for --sge:
#   PROJECT    SGE project (default: stat_bhaswar_team)
#   MEM        m_mem_free request (default: 8G)
#   N_SETUPS   number of setup indices to dispatch (default: 60)

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: submit-ccle.sh {--local|--sge}

  --local  Run all setup indices synchronously in the current process.
  --sge    Submit one Grid Engine job per setup and return immediately.
USAGE
}

if [[ $# -ne 1 ]]; then
  usage >&2
  exit 2
fi

mode=$1
if [[ "$mode" != "--local" && "$mode" != "--sge" ]]; then
  usage >&2
  exit 2
fi

if [[ -z "${LOCAL_PDRT_DATA_DIR:-}" ]]; then
  echo "ERROR: LOCAL_PDRT_DATA_DIR is not set." >&2
  exit 1
fi

data_root=${LOCAL_PDRT_DATA_DIR%/}
ccle_dir="$data_root/ccle"
for input_file in mutation.txt response.csv; do
  if [[ ! -f "$ccle_dir/$input_file" ]]; then
    echo "ERROR: Missing CCLE input: $ccle_dir/$input_file" >&2
    exit 1
  fi
done

n_setups=${N_SETUPS:-60}
if [[ ! "$n_setups" =~ ^[1-9][0-9]*$ ]] || (( n_setups > 60 )); then
  echo "ERROR: N_SETUPS must be an integer from 1 through 60." >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd -- "$script_dir/../.." && pwd)
runner="$script_dir/ccle-typeI-synthY-30targets.R"

if [[ ! -f "$runner" ]]; then
  echo "ERROR: Runner script not found: $runner" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available on PATH." >&2
  exit 1
fi

echo "Project root:        $project_root"
echo "Runner script:       $runner"
echo "Number of setups:    $n_setups"
echo "LOCAL_PDRT_DATA_DIR: $data_root"

if [[ "$mode" == "--local" ]]; then
  cd "$project_root"
  for ((setup_idx = 1; setup_idx <= n_setups; setup_idx++)); do
    echo "Running CCLE setup $setup_idx/$n_setups ..."
    Rscript "$runner" "$setup_idx"
  done
  echo "Completed $n_setups CCLE setup(s)."
  exit 0
fi

if ! command -v qsub >/dev/null 2>&1; then
  echo "ERROR: qsub is not available on PATH." >&2
  exit 1
fi

sge_project=${PROJECT:-stat_bhaswar_team}
memory=${MEM:-8G}
echo "SGE project:         $sge_project"
echo "Memory per job:      $memory"

for ((setup_idx = 1; setup_idx <= n_setups; setup_idx++)); do
  job_name="ccle_synthY_${setup_idx}"
  printf -v job_body \
    '#!/usr/bin/env bash\nset -euo pipefail\ncd -- %q\nexec Rscript %q %q\n' \
    "$project_root" "$runner" "$setup_idx"
  printf '%s' "$job_body" | qsub \
    -V \
    -P "$sge_project" \
    -l "m_mem_free=$memory" \
    -cwd \
    -N "$job_name"
done

echo "Submitted $n_setups CCLE setup job(s)."
echo "After they finish, run: ./run_all_semi_synthetic.sh plot-only"
echo "The root runner will combine the setup files before plotting."
