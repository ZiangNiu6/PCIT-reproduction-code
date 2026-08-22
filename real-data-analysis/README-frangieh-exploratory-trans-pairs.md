# Frangieh exploratory trans-pair scan

This standalone workflow scans low-expression trans pairs without placing the
positive or negative controls from the manuscript experiment in the candidate
set. It does not modify or feed the production Frangieh FDR workflow.

## Fixed candidate universe

- Cells: the Frangieh Control condition, restricted during testing to cells
  carrying exactly one recognized guide.
- Primary exposures: 241 gene targets with at least ten treated cells, at
  least two treated batches, and no batch containing more than 80% of treated
  cells. ATP1A1, DNAJC9, EIF2S3, PSMA7, RACK1, TUBB, and UBL5 fail the
  ten-cell threshold. `--min-treated 5` retains 244 exposures as a declared
  sensitivity analysis.
- Responses: exclude all 248 perturbation-target genes; require detection in
  at least 577 of all 57,624 local Control-condition cells; rank by the
  precomputed global mean count, breaking ties by gene name; retain the lowest
  5,000.
- Full primary grid: 241 x 5,000 = 1,205,000 trans pairs. No self-target pair or
  non-targeting-guide exposure can occur.

The detection filter and expression ranking do not inspect cell guide
assignments; guide target annotation is used only to exclude all perturbation
targets from the response universe. The panel is therefore unchanged under
every pairwise guide-label permutation. Stored
global means and detection counts were validated against counts reconstructed
from the raw H5AD: all 57,624 cells and 23,712 genes agreed (maximum mean
difference 2.78e-17 and zero nonzero-count differences; maximum inversion
residual 0.00164278), and every shard repeats deterministic raw-count audits.
In the 3,770 exact-one non-targeting cells, the fixed response panel has
20--439 nonzero cells (median 137), with no all-zero response.

The expression ranking is pooled over **all** 57,624 Control-condition cells;
it is not a pure-non-targeting baseline. Pooling without inspecting guide
assignments makes selection label-invariant for this permutation analysis, but
it can still enrich or deplete genes that respond broadly to perturbation.
That biological feature of the exploratory panel should be reported rather
than interpreted as a null-only expression screen.

The tractable primary pilot is a frozen association-blind sample of 5,000
pairs. Every response appears exactly once, and every one of the 241 exposures
receives 20 or 21 responses. A fixed seed assigns the pairs using only the
locked exposure and response manifests—never an association, fold change,
score, or p-value.

## Local commands

Create and audit only the fixed manifests:

```bash
Rscript --vanilla reproduction-code/real-data-analysis/run_frangieh_exploratory_trans_pairs.R \
  --h5ad /absolute/path/frangieh_2021_raw.h5ad \
  --intermediate-dir /absolute/path/frangieh-2021/intermediate-data \
  --output-dir results/frangieh-2021/exploratory-trans-pairs \
  --manifest-only
```

Run this once in the final output directory before launching array tasks. It
performs the across-all-Control raw-count metadata recheck and writes the
locked expression-audit manifest required by every shard. This one-time step
also hashes the complete H5AD. Array shards avoid repeatedly hashing the 6 GB
file: they require the same normalized path, byte size, and UTC modification
time, then inherit the locked full-file MD5 into their provenance. If the file
is copied to a different path or its identity changes, rerun `--manifest-only`
in a new output directory.

Run one perturbation's frozen-pilot shard:

```bash
Rscript --vanilla reproduction-code/real-data-analysis/run_frangieh_exploratory_trans_pairs.R \
  --h5ad /absolute/path/frangieh_2021_raw.h5ad \
  --intermediate-dir /absolute/path/frangieh-2021/intermediate-data \
  --output-dir results/frangieh-2021/exploratory-trans-pairs \
  --target-index 1 --pilot --min-treated 10 --response-block-size 100 \
  --permutations 4999 --seed 20260822 --cores 4
```

`--target GENE` and `--target-range START:END` are alternatives to
`--target-index`. Without `--pilot`, `--response-range START:END` restricts a
full-grid shard to fixed response ranks and is useful for smoke tests or
two-dimensional sharding. The fixed pilot does not accept a response range.
Seeds are a deterministic function of the global exposure index and response
rank, so target order, response block size, and core count do not change the
permutation draws.

## HPCC array

The launcher defaults to the frozen pilot and submits SGE tasks 1--241, one
task per primary eligible perturbation:

```bash
export FRANGIEH_RAW_H5AD=/absolute/path/frangieh_2021_raw.h5ad
export FRANGIEH_INTERMEDIATE_DIR=/absolute/path/frangieh-2021/intermediate-data
export FRANGIEH_EXPLORATORY_RESULTS_DIR=/absolute/path/exploratory-trans-pairs
export PERMUTATIONS=4999
export SEED=20260822
export CORES=1
export MIN_TREATED=10
export SCAN_MODE=pilot
bash reproduction-code/real-data-analysis/submit_frangieh_exploratory_trans_pairs.sh
```

For multiple cores per task, also set `SGE_PE` to the parallel-environment
name configured on the cluster. `SCAN_MODE=full` explicitly requests the
1,205,000-pair grid. `MIN_TREATED=5` is supported only for a full-grid
sensitivity run and produces 244 array tasks.

## Outputs

Every shard writes a CSV and RDS with observed statistics, both one-sided tail
hit counts, the absolute two-sided hit counts, Barry and empirical p-values,
pair seeds, diagnostics, and failure status. The RDS embeds the fixed response
and exposure manifests, data audits, input/script hashes, package versions,
and session information. Byte-identical shared CSV manifests live in the
`manifests/` subdirectory.

The p-values are raw. Do not apply BH within a shard. A complete-shard combiner
must first verify manifest hashes, unique pair IDs, expected coverage, and
provenance, then correct across the full declared pair universe. With 1.205
million hypotheses in the primary full grid, 4,999 permutations have limited
tail resolution. The limitation already applies to the 5,000-pair pilot: with
`B=4999`, the smallest raw p-value is 1/5,000 = 0.0002. BH at FDR 0.05 would
require at least 20 hypotheses tied at that resolution floor before it could
reject any, without resolving their relative tail evidence. Treat that as
inadequate for discovery claims. Use the
fixed scan primarily for exploratory ranking or add a valid adaptive/deeper
permutation stage for multiplicity-controlled discoveries.

## Strict combination

After every frozen-pilot array task finishes, combine only with the strict
coverage gate:

```bash
Rscript --vanilla reproduction-code/real-data-analysis/combine_frangieh_exploratory_trans_pairs.R \
  --input-dir /absolute/path/exploratory-trans-pairs \
  --output-dir /absolute/path/exploratory-trans-pairs/combined \
  --mode pilot --min-treated 10 --permutations 4999 --seed 20260822
```

The combiner requires identical manifest hashes, runner/helper hashes, input
hashes, locked H5AD identity/full MD5, package versions, RNG settings, `B`, and
seed. It row-matches every result to the declared pair, exposure, response,
and frozen-pilot manifests; checks seeds, hit counts, and p-value formulas;
rejects duplicate or non-OK pairs; and refuses output unless the pilot has exactly
5,000 unique responses with 20--21 pairs for each of all 241 exposures.
`--mode full` similarly requires exact coverage of all 1,205,000 primary-grid
pairs (or 1,220,000 pairs with `--min-treated 5`).

Successful combination writes the complete RDS/CSV, a paired raw-p comparison
summary, and a compact plot-ready table. It deliberately does not calculate
BH-adjusted p-values or call discoveries: the combined artifact records the
minimum attainable permutation p-value and an explicit tail-resolution
warning.

By default existing combined files are never replaced. `--overwrite` is
available, but replacement of the four-file bundle is sequential rather than
transactional; for an archival run, prefer a fresh combined-output directory.
