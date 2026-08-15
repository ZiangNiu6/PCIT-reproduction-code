#!/usr/bin/env bash
# Run or submit 42 matched-v4 CCLE null/representative-effect jobs.
# Power representatives are RIOK3_MUT (high sparsity), PLK1_MUT (medium), and
# ATP11A_MUT (low sparsity); all 30 targets retain their delta = 0 job.
#
# Usage:
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --local
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --sge
#
# With no argument, the script defaults to --sge.
# --local runs each setup synchronously. It is used by the repository-level
# full runner so combining and plotting begin only after computation succeeds.
# --sge submits independent Grid Engine jobs and returns after submission. The
# jobs write only the matched-v4 per-setup namespace; after they finish, use
# the repository-level plot-only workflow to combine and plot those results.
#
# Optional environment variables for --sge:
#   PROJECT    SGE project (default: stat_bhaswar_team)
#   MEM        m_mem_free request (default: 8G)
#   N_SETUPS   number of setup indices to dispatch (default: 42)

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: submit-ccle.sh [--local|--sge]

  no option  Submit one Grid Engine job per setup (same as --sge).
  --local  Run all setup indices synchronously in the current process.
  --sge    Submit one Grid Engine job per setup and return immediately.
USAGE
}

if [[ $# -gt 1 ]]; then
  usage >&2
  exit 2
fi

mode=${1:---sge}
if [[ "$mode" != "--local" && "$mode" != "--sge" ]]; then
  usage >&2
  exit 2
fi

# Load the lab configuration when available. Required values are embedded
# explicitly in each job script; the full login environment is not sent to
# Grid Engine.
if [[ -f "$HOME/.research_config" ]]; then
  # shellcheck disable=SC1090
  source "$HOME/.research_config"
fi

if [[ -z "${LOCAL_PDRT_DATA_DIR:-}" ]]; then
  echo "ERROR: LOCAL_PDRT_DATA_DIR is not set and was not loaded from ~/.research_config." >&2
  exit 1
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd -- "$script_dir/../.." && pwd)
runner="$script_dir/ccle-typeI-synthY-30targets.R"

data_root=$LOCAL_PDRT_DATA_DIR
if [[ "$data_root" != /* ]]; then
  data_root="$project_root/$data_root"
fi
if [[ "$data_root" != "/" ]]; then
  data_root=${data_root%/}
fi
export LOCAL_PDRT_DATA_DIR="$data_root"
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

ccle_dir="$data_root/ccle"
for input_file in mutation.txt response.csv; do
  if [[ ! -f "$ccle_dir/$input_file" ]]; then
    echo "ERROR: Missing CCLE input: $ccle_dir/$input_file" >&2
    exit 1
  fi
done

n_setups=${N_SETUPS:-42}
if [[ ! "$n_setups" =~ ^[1-9][0-9]*$ ]] || (( n_setups > 42 )); then
  echo "ERROR: N_SETUPS must be an integer from 1 through 42." >&2
  exit 2
fi

if [[ ! -f "$runner" ]]; then
  echo "ERROR: Runner script not found: $runner" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available on PATH." >&2
  exit 1
fi
rscript_bin=$(command -v Rscript)

# Jobs run plain `Rscript` after `cd`-ing to the project root, so the renv
# autoloader (project .Rprofile -> renv/activate.R) puts the project library on
# the path automatically. Only the required data path is embedded in the job.
echo "Project root:        $project_root"
echo "Runner script:       $runner"
echo "Result namespace:    matched-v4-representative-effect-grid-s10"
echo "Number of setups:    $n_setups"
echo "LOCAL_PDRT_DATA_DIR: $data_root"
echo "Rscript:             $rscript_bin"

if [[ "$mode" == "--local" ]]; then
  cd "$project_root"
  for ((setup_idx = 1; setup_idx <= n_setups; setup_idx++)); do
    echo "Running CCLE setup $setup_idx/$n_setups ..."
    "$rscript_bin" "$runner" "$setup_idx"
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
  job_name="ccle_matched_v4_${setup_idx}"
  printf -v job_body \
    '#!/usr/bin/env bash\nset -euo pipefail\nexport LOCAL_PDRT_DATA_DIR=%q\nexport RENV_CONFIG_SANDBOX_ENABLED=FALSE\ncd -- %q\nexec %q %q %q\n' \
    "$data_root" "$project_root" "$rscript_bin" "$runner" "$setup_idx"
  printf '%s' "$job_body" | qsub \
    -P "$sge_project" \
    -l "m_mem_free=$memory" \
    -cwd \
    -N "$job_name"
done

echo "Submitted $n_setups CCLE setup job(s)."
echo "After they finish, run: ./run_all_semi_synthetic.sh plot-only"
echo "If the matched-v4 combined result is absent, the root runner will"
echo "strictly validate and combine all 42 matched-v4 setup files before plotting."
