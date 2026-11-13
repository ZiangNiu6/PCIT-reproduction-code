#!/usr/bin/env bash
# Fresh computation must refuse saved/interrupted work before dispatch.
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
expect_fresh_refusal() {
  expect_failure 1 "$@"
  if ! grep -q 'fresh' "$temporary_dir/command.log"; then
    cat "$temporary_dir/command.log" >&2
    echo "Expected a fresh-directory refusal before computation." >&2
    exit 1
  fi
}

expect_failure 2 bash "$PACKAGE_DIR/run_all_simulations.sh" nonexistent-mode
expect_failure 2 bash "$SIM_DIR/run_simulation.sh" \
  --sim_name pcm-additive --results_dir "$temporary_dir/not-created" --max_hours 0
[[ ! -e "$temporary_dir/not-created" ]]

# An R invocation here would prove the refusal happened too late.
mkdir "$temporary_dir/no-compute-bin"
cat > "$temporary_dir/no-compute-bin/Rscript" <<'SH'
#!/usr/bin/env bash
printf 'unexpected R invocation\n' >> "$PCIT_TEST_R_LOG"
exit 99
SH
chmod +x "$temporary_dir/no-compute-bin/Rscript"
export PCIT_TEST_R_LOG="$temporary_dir/unexpected-r.log"

# Parent full mode checks the entire domain before touching any study.
mkdir -p "$temporary_dir/full-busy/simulation-results/previous-study"
printf 'existing output\n' > "$temporary_dir/full-busy/simulation-results/previous-study/result.rds"
expect_fresh_refusal env PATH="$temporary_dir/no-compute-bin:$PATH" \
  LOCAL_PDRT_DATA_DIR="$temporary_dir/full-busy" \
  bash "$PACKAGE_DIR/run_all_simulations.sh" full
[[ ! -e "$PCIT_TEST_R_LOG" ]]
[[ ! -e "$temporary_dir/full-busy/simulation-results/.full-run-started" ]]

# Study-level entry points also reject results, serialized functions and work.
for study in nb-permutation-signflip pcm-additive; do
  for existing in result specifier work claim; do
    directory="$temporary_dir/$study-$existing"
    mkdir "$directory"
    case "$existing" in
      result) entry="$directory/${study}_results.rds" ;;
      specifier) entry="$directory/${study}_specifier.rds" ;;
      work) entry="$directory/.nextflow-work" ;;
      claim) entry="$directory/.full-run-started" ;;
    esac
    printf 'existing data must remain unchanged\n' > "$entry"
    cp "$entry" "$temporary_dir/original"
    expect_fresh_refusal env PATH="$temporary_dir/no-compute-bin:$PATH" \
      bash "$SIM_DIR/run_simulation.sh" --sim_name "$study" --results_dir "$directory"
    cmp "$entry" "$temporary_dir/original"
    [[ ! -e "$PCIT_TEST_R_LOG" ]]
  done
done

source "$PACKAGE_DIR/workflow-common.sh"
pcit_initialize "$PACKAGE_DIR"
pcit_check_r

# Numerical grid definitions are always read from current R source.
expect_failure 1 Rscript --vanilla "$SIM_DIR/figure-code/plot-nb-permutation-signflip.R" \
  "$temporary_dir/unused-results.rds" "$temporary_dir/legacy-grid.rds" --validate-only
grep -q 'must be an R source file' "$temporary_dir/command.log"

# Direct GLM computation refuses both an existing result and prior output work.
for existing in result work claim; do
  directory="$temporary_dir/glm-$existing"
  mkdir "$directory"
  result="$directory/illustrative-example-results.rds"
  case "$existing" in
    result) entry="$result" ;;
    work) entry="$directory/session-info.txt" ;;
    claim) entry="$directory/.full-run-started" ;;
  esac
  printf 'existing GLM work must remain unchanged\n' > "$entry"
  cp "$entry" "$temporary_dir/original"
  expect_fresh_refusal env PCIT_GLM_RESULTS_FILE="$result" \
    Rscript --vanilla "$SIM_DIR/glm-variance-robustness.R" "--output-dir=$directory"
  cmp "$entry" "$temporary_dir/original"
done

# Exercise fresh dispatch without running a numerical study: a tiny spec is
# built from current source, then a stub Nextflow stops with its own exit code.
mkdir "$temporary_dir/specs" "$temporary_dir/dispatch-bin"
cat > "$temporary_dir/specs/nb-permutation-signflip.R" <<'RS'
simulatr_spec <- list(test_marker = "built from current source")
RS
printf '# Test grid placeholder.\n' > "$temporary_dir/specs/nb-permutation-signflip-grid.R"
cat > "$temporary_dir/dispatch-bin/nextflow" <<'SH'
#!/usr/bin/env bash
if [[ "$1" == pull ]]; then exit 0; fi
exit 42
SH
chmod +x "$temporary_dir/dispatch-bin/nextflow"
expect_failure 42 env PATH="$temporary_dir/dispatch-bin:$PATH" \
  bash "$SIM_DIR/run_simulation.sh" --sim_name nb-permutation-signflip \
  --sim_spec_dir "$temporary_dir/specs" --results_dir "$temporary_dir/fresh-study"
[[ -d "$temporary_dir/fresh-study/.full-run-started" ]]
PCIT_TEST_SPEC="$temporary_dir/fresh-study/nb-permutation-signflip_specifier.rds" \
  pcit_r -e 'stopifnot(identical(readRDS(Sys.getenv("PCIT_TEST_SPEC"))$test_marker, "built from current source"))'
expect_fresh_refusal env PATH="$temporary_dir/no-compute-bin:$PATH" \
  bash "$SIM_DIR/run_simulation.sh" --sim_name nb-permutation-signflip \
  --results_dir "$temporary_dir/fresh-study"
[[ ! -e "$PCIT_TEST_R_LOG" ]]

echo "Simulation fresh-run and no-overwrite workflow checks passed."
