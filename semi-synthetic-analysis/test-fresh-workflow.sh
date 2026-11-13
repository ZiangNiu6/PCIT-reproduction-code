#!/usr/bin/env bash
# Exercise fresh-run/dispatch contracts with mocked R and SGE, never real jobs.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
package_root=$(cd -- "$script_dir/.." && pwd)
sandbox=$(mktemp -d "${TMPDIR:-/tmp}/ccle-fresh-workflow.XXXXXX")
trap 'rm -rf -- "$sandbox"' EXIT
mkdir "$sandbox/bin" "$sandbox/jobs"
cat > "$sandbox/bin/Rscript" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${2:-}" == - ]]; then cat >/dev/null; exit 0; fi
printf '%s\n' "$*" >> "$CCLE_TEST_R_LOG"
MOCK
cat > "$sandbox/bin/qsub" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
name=
while (($#)); do
  if [[ "$1" == -N ]]; then name="$2"; break; fi
  shift
done
[[ -n "$name" ]]
cat > "$CCLE_TEST_JOB_DIR/$name"
MOCK
chmod +x "$sandbox/bin/Rscript" "$sandbox/bin/qsub"
export PATH="$sandbox/bin:$PATH"
export CCLE_TEST_R_LOG="$sandbox/r-calls"
export CCLE_TEST_JOB_DIR="$sandbox/jobs"
export PCIT_USE_RESEARCH_CONFIG=0 PCIT_RESTORE=0
unset PCIT_CONFIG_FILE PCIT_R_LIBRARY R_LIBS_USER
prepare_inputs() {
  mkdir -p "$1/ccle"
  touch "$1/ccle/mutation.txt" "$1/ccle/response.csv"
}
expect_failure() {
  if "$@" > "$sandbox/failure.log" 2>&1; then
    echo "Expected fresh-workflow refusal: $*" >&2; exit 1
  fi
}
# Full runs dispatch every setting once; no resume or forced replacement.
prepare_inputs "$sandbox/local"
LOCAL_PDRT_DATA_DIR="$sandbox/local" bash "$package_root/run_all_semi_synthetic.sh" full > "$sandbox/full.log"
[[ -d "$sandbox/local/ccle/.full-run-started" ]]
[[ -d "$sandbox/local/ccle/.ccle-dispatched" ]]
[[ $(grep -c '/ccle-typeI-synthY-30targets[.]R ' "$CCLE_TEST_R_LOG") == 42 ]]
before_calls=$(wc -l < "$CCLE_TEST_R_LOG")
expect_failure env LOCAL_PDRT_DATA_DIR="$sandbox/local" bash "$package_root/run_all_semi_synthetic.sh" full
[[ $(wc -l < "$CCLE_TEST_R_LOG") == "$before_calls" ]]
expect_failure bash "$package_root/run_all_semi_synthetic.sh" full --force
expect_failure bash "$script_dir/submit-ccle.sh" --force
# Existing files are rejected before computation and remain unchanged.
prepare_inputs "$sandbox/occupied"
printf 'untouched\n' > "$sandbox/occupied/ccle/prior-output.rds"
expect_failure env LOCAL_PDRT_DATA_DIR="$sandbox/occupied" bash "$package_root/run_all_semi_synthetic.sh" full
[[ $(cat "$sandbox/occupied/ccle/prior-output.rds") == untouched ]]
[[ ! -e "$sandbox/occupied/ccle/.full-run-started" ]]
# SGE submission writes42distinct worker templates; immediate redispatch fails.
prepare_inputs "$sandbox/sge"
LOCAL_PDRT_DATA_DIR="$sandbox/sge" bash "$script_dir/submit-ccle.sh" --sge > "$sandbox/sge.log"
for setup_idx in {1..42}; do [[ -s "$sandbox/jobs/ccle_matched_v4_$setup_idx" ]]; done
(( $(find "$sandbox/jobs" -type f | wc -l) == 42 ))
expect_failure env LOCAL_PDRT_DATA_DIR="$sandbox/sge" bash "$script_dir/submit-ccle.sh" --sge
(( $(find "$sandbox/jobs" -type f | wc -l) == 42 ))
# A standalone dispatcher also refuses an already completed result.
prepare_inputs "$sandbox/completed"
printf 'untouched\n' > "$sandbox/completed/ccle/ccle-typeI-synthY-30targets-matched-v4-results-EXAMPLE-delta0-s10.rds"
expect_failure env LOCAL_PDRT_DATA_DIR="$sandbox/completed" bash "$script_dir/submit-ccle.sh" --local
[[ ! -e "$sandbox/completed/ccle/.ccle-dispatched" ]]
echo 'CCLE fresh-only workflow and42-job dispatch checks passed; no analyses executed.'
