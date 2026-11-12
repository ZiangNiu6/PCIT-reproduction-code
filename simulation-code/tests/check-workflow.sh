#!/usr/bin/env bash
# Regression checks for refusing malformed cached results and invalid requests.
set -euo pipefail
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SIM_DIR="$(cd "$TEST_DIR/.." && pwd -P)"
PACKAGE_DIR="$(cd "$SIM_DIR/.." && pwd -P)"
temporary_dir=$(mktemp -d "${TMPDIR:-/tmp}/pcit-simulation-workflow.XXXXXX")
trap 'rm -rf "$temporary_dir"' EXIT

expect_failure() {
  local expected_status="$1"
  shift
  local actual_status=0
  "$@" >"$temporary_dir/command.log" 2>&1 || actual_status=$?
  if [[ "$actual_status" != "$expected_status" ]]; then
    cat "$temporary_dir/command.log" >&2
    echo "Expected exit $expected_status; got $actual_status." >&2
    exit 1
  fi
}

expect_failure 2 bash "$PACKAGE_DIR/run_all_simulations.sh" nonexistent-mode
expect_failure 2 bash "$SIM_DIR/run_simulation.sh" \
  --sim_name pcm-additive --results_dir "$temporary_dir/not-created" --max_hours 0
[[ ! -e "$temporary_dir/not-created" ]]

mkdir "$temporary_dir/nb" "$temporary_dir/pcm"
for study in nb-permutation-signflip pcm-additive; do
  if [[ "$study" == nb-permutation-signflip ]]; then folder=nb; else folder=pcm; fi
  result="$temporary_dir/$folder/${study}_results.rds"
  printf 'malformed result that must not be reused or overwritten\n' > "$result"
  cp "$result" "$temporary_dir/original"
  expect_failure 1 bash "$SIM_DIR/run_simulation.sh" \
    --sim_name "$study" --results_dir "$temporary_dir/$folder"
  cmp "$result" "$temporary_dir/original"
done

echo "Simulation workflow regression checks passed."
