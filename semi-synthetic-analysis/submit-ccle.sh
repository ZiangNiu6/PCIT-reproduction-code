#!/usr/bin/env bash
# Run or submit 42 matched-v4 CCLE null/representative-effect jobs.
# Power representatives are EIF5A_MUT (high sparsity), PLK1_MUT (medium), and
# ATP11A_MUT (low sparsity); all 30 targets retain their delta = 0 job.
#
# Usage:
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --local
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --sge
#   bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh --local --force
#
# With no argument, the script defaults to --sge.
# --local runs each setup synchronously. It is used by the repository-level
# full runner so combining and plotting begin only after computation succeeds.
# --sge submits independent Grid Engine jobs and returns after submission. The
# jobs write only the matched-v4 per-setup namespace; after they finish, use
# the repository-level plot-only workflow to combine and plot those results.
#
# Optional environment variables for --sge:
#   PROJECT    optional SGE project/account
#   MEM        m_mem_free request (default: 8G)
#   N_SETUPS   number of setup indices to dispatch (default: 42)
# For either mode, --force (or CCLE_FORCE_RERUN=1) explicitly reruns every
# selected setup and permits its worker to atomically replace an existing
# expected output. Without it, valid outputs are skipped and malformed outputs
# fail closed without being changed.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: submit-ccle.sh [--local|--sge] [--force]

  no option  Submit one Grid Engine job per setup (same as --sge).
  --local    Run only missing setup indices synchronously.
  --sge      Submit one Grid Engine job per missing setup index and return.
  --force    Rerun every selected setup, replacing reviewed existing outputs.
USAGE
}

if [[ $# -gt 2 ]]; then
  usage >&2
  exit 2
fi

mode=--sge
mode_was_set=0
force_rerun=${CCLE_FORCE_RERUN:-0}
if [[ "$force_rerun" != "0" && "$force_rerun" != "1" ]]; then
  echo "ERROR: CCLE_FORCE_RERUN must be exactly 0 or 1." >&2
  exit 2
fi
for option in "$@"; do
  case "$option" in
    --local|--sge)
      if (( mode_was_set )); then
        usage >&2
        exit 2
      fi
      mode=$option
      mode_was_set=1
      ;;
    --force)
      force_rerun=1
      ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
done

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
package_root=$(cd -- "$script_dir/.." && pwd)
# shellcheck source=../workflow-common.sh
source "$package_root/workflow-common.sh"
pcit_initialize "$package_root"
pcit_check_r
runner="$script_dir/ccle-typeI-synthY-30targets.R"
combiner="$script_dir/combine-ccle-results.R"

data_root=$LOCAL_PDRT_DATA_DIR
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
if [[ ! "$n_setups" =~ ^([1-9]|[1-3][0-9]|4[0-2])$ ]]; then
  echo "ERROR: N_SETUPS must be a canonical integer from 1 through 42." >&2
  exit 2
fi
n_setups_decimal=$((10#$n_setups))
n_setups=$n_setups_decimal

if [[ ! -f "$runner" ]]; then
  echo "ERROR: Runner script not found: $runner" >&2
  exit 1
fi
if [[ ! -f "$combiner" ]]; then
  echo "ERROR: Combiner/status script not found: $combiner" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available on PATH." >&2
  exit 1
fi
rscript_bin=$(command -v Rscript)

# Local and SGE workers use the same explicitly selected package library.
# Only paths needed by each worker are embedded; no full environment export.
echo "Package root:        $package_root"
echo "Runner script:       $runner"
echo "Result namespace:    matched-v4-representative-effect-grid-s10"
echo "Number of setups:    $n_setups"
echo "LOCAL_PDRT_DATA_DIR: $data_root"
echo "Rscript:             $rscript_bin"
echo "Force rerun:         $force_rerun"

status_file=$(mktemp "${TMPDIR:-/tmp}/ccle-v4-output-status.XXXXXX")
cleanup_status_file() {
  rm -f -- "$status_file"
}
trap cleanup_status_file EXIT

cd "$package_root"
pcit_r "$combiner" --output-status "$status_file" "$n_setups" \
  >/dev/null

pending_indices=()
valid_count=0
missing_count=0
invalid_count=0
invalid_detail_limit=10
status_count=0
while IFS=$'\t' read -r setup_idx state filename reason; do
  ((status_count += 1))
  case "$state" in
    valid)
      ((valid_count += 1))
      if (( force_rerun )); then
        pending_indices+=("$setup_idx")
      fi
      ;;
    missing)
      ((missing_count += 1))
      pending_indices+=("$setup_idx")
      ;;
    invalid)
      ((invalid_count += 1))
      if (( force_rerun )); then
        if (( invalid_count <= invalid_detail_limit )); then
          echo "WARNING: setup $setup_idx has an invalid expected output:" >&2
          echo "  $filename" >&2
          echo "  $reason" >&2
        fi
        pending_indices+=("$setup_idx")
      else
        if (( invalid_count <= invalid_detail_limit )); then
          echo "ERROR: setup $setup_idx has an invalid expected output:" >&2
          echo "  $filename" >&2
          echo "  $reason" >&2
        fi
      fi
      ;;
    *)
      echo "ERROR: Unrecognized output-status state for setup $setup_idx: $state" >&2
      exit 1
      ;;
  esac
done < "$status_file"

if (( invalid_count > invalid_detail_limit )); then
  echo "... $((invalid_count - invalid_detail_limit)) more invalid output(s) omitted." >&2
fi

if (( status_count != n_setups )); then
  echo "ERROR: Combiner returned $status_count status rows; expected $n_setups." >&2
  exit 1
fi
if (( invalid_count > 0 && ! force_rerun )); then
  echo "ERROR: Refusing to overwrite $invalid_count invalid output(s)." >&2
  echo "Review them, then use --force or CCLE_FORCE_RERUN=1 to rerun explicitly." >&2
  exit 1
fi

pending_count=${#pending_indices[@]}
echo "Validated existing:  $valid_count"
echo "Missing outputs:     $missing_count"
echo "Invalid outputs:     $invalid_count"
if (( force_rerun )); then
  echo "Setups forced:       $pending_count"
else
  echo "Setups to run:       $pending_count"
fi

if (( pending_count == 0 )); then
  echo "All $n_setups selected CCLE setup outputs are already valid; nothing to run."
  exit 0
fi

export CCLE_FORCE_RERUN="$force_rerun"

if [[ "$mode" == "--local" ]]; then
  pending_position=0
  for setup_idx in "${pending_indices[@]}"; do
    ((pending_position += 1))
    echo "Running CCLE setup $setup_idx ($pending_position/$pending_count pending) ..."
    pcit_r "$runner" "$setup_idx"
  done
  if (( force_rerun )); then
    echo "Completed $pending_count explicitly forced CCLE setup(s)."
  else
    echo "Completed $pending_count CCLE setup(s); skipped $valid_count valid output(s)."
  fi
  exit 0
fi

if ! command -v qsub >/dev/null 2>&1; then
  echo "ERROR: qsub is not available on PATH." >&2
  exit 1
fi

sge_project=${PROJECT:-}
memory=${MEM:-8G}
echo "SGE project:         ${sge_project:-scheduler default}"
echo "Memory per job:      $memory"
sge_options=(-l "m_mem_free=$memory" -cwd)
if [[ -n "$sge_project" ]]; then sge_options+=(-P "$sge_project"); fi

for setup_idx in "${pending_indices[@]}"; do
  job_name="ccle_matched_v4_${setup_idx}"
  printf -v job_body \
    '#!/usr/bin/env bash\nset -euo pipefail\nexport PATH=%q:"$PATH"\nexport LOCAL_PDRT_DATA_DIR=%q\nexport CCLE_FORCE_RERUN=%q\nexport PCIT_R_LIBRARY=%q\nsource %q\npcit_initialize %q\npcit_check_r\npcit_r %q %q\n' \
    "$(dirname -- "$rscript_bin")" "$data_root" "$force_rerun" "${PCIT_R_LIBRARY:-}" \
    "$package_root/workflow-common.sh" "$package_root" "$runner" "$setup_idx"
  printf '%s' "$job_body" | qsub \
    "${sge_options[@]}" \
    -N "$job_name"
done

if (( force_rerun )); then
  echo "Submitted $pending_count explicitly forced CCLE setup job(s)."
else
  echo "Submitted $pending_count CCLE setup job(s); skipped $valid_count valid output(s)."
fi
echo "After they finish, combine explicitly: Rscript semi-synthetic-analysis/combine-ccle-results.R"
echo "Then run: ./run_all_semi_synthetic.sh plot-only"
echo "Plot-only uses an existing combined artifact without raw-data recombination."
