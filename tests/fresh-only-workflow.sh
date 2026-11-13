#!/usr/bin/env bash
# Check fresh-only execution without starting any scientific computation.
set -euo pipefail
PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "$PACKAGE_DIR/workflow-common.sh"
fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/pcit-fresh-only.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT

expect_refusal() {
  local status=0
  "$@" >"$fixture_dir/refusal.log" 2>&1 || status=$?
  if [[ "$status" != 1 ]] || ! grep -q 'fresh output directory' "$fixture_dir/refusal.log"; then
    cat "$fixture_dir/refusal.log" >&2
    echo "Expected fresh-output refusal; got exit $status." >&2
    exit 1
  fi
}

# Both successful and interrupted launches permanently claim their directory.
pcit_require_fresh_directory "$fixture_dir/empty"
[[ -d "$fixture_dir/empty/.full-run-started" ]]
expect_refusal pcit_require_fresh_directory "$fixture_dir/empty"

# Preparing raw CCLE inputs must remain compatible with a fresh full run.
mkdir "$fixture_dir/raw"
printf 'raw mutation input\n' > "$fixture_dir/raw/mutation.txt"
printf 'raw response input\n' > "$fixture_dir/raw/response.csv"
pcit_require_fresh_directory "$fixture_dir/raw" mutation.txt response.csv input-provenance.json
[[ -d "$fixture_dir/raw/.full-run-started" ]]
[[ "$(cat "$fixture_dir/raw/mutation.txt")" == 'raw mutation input' ]]

# Relative library overrides resolve against the package from any caller cwd.
mkdir -p "$fixture_dir/package/r-lib"
(
  cd "$fixture_dir"
  PCIT_USE_RESEARCH_CONFIG=0
  unset PCIT_CONFIG_FILE
  pcit_initialize "$fixture_dir/package"
  export PCIT_R_LIBRARY=r-lib
  pcit_check_r
  [[ "$PCIT_R_LIBRARY" == "$fixture_dir/package/r-lib" ]]
)

# A top-level full request must refuse occupied domains before package checks,
# data downloads, model fits or scheduler submissions. Even malformed saved
# output is preserved rather than inspected, reused or overwritten.
for entry in 'simulation-results:run_all_simulations.sh' \
             'ccle:run_all_semi_synthetic.sh' \
             'frangieh-2021:run_all_real_data.sh'; do
  domain="${entry%%:*}"
  runner="${entry#*:}"
  results="$fixture_dir/$domain-root"
  mkdir -p "$results/$domain"
  if [[ "$domain" == ccle ]]; then
    printf 'prepared input fixture\n' > "$results/$domain/mutation.txt"
    printf 'prepared input fixture\n' > "$results/$domain/response.csv"
  fi
  printf 'existing output must remain untouched\n' > "$results/$domain/sentinel.rds"
  cp "$results/$domain/sentinel.rds" "$fixture_dir/expected"
  expect_refusal env -u PCIT_CONFIG_FILE -u FRANGIEH_INTERMEDIATE_DIR \
    -u FRANGIEH_SCORE_RESULTS_DIR -u FRANGIEH_SCORE_PLOTS_DIR \
    LOCAL_PDRT_DATA_DIR="$results" PCIT_USE_RESEARCH_CONFIG=0 \
    PCIT_RESTORE=0 bash "$PACKAGE_DIR/$runner" full
  cmp "$results/$domain/sentinel.rds" "$fixture_dir/expected"
  [[ ! -e "$results/$domain/.full-run-started" ]]
done

echo 'Fresh-only workflow checks passed; no scientific computation was started.'
