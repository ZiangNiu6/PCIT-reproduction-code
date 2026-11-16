# Frangieh reproduction workflow

The analysis and bundled results use **BH q=0.10**, GB2024 preprocessing,
and schema `perturb_score5_v4_gb2024_q10`.

## Preprocessing contract

The implementation follows the historical
[Barry et al. reproduction script at commit `005c28b`](https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/005c28b9ba30e7d39c3ee18d51111722f49d32ec/R_scripts/uniform_processing_lowmoi.R),
including its pinned `ondisc` behavior:

- Count all original assigned guides before feature eligibility checks;
  retain exactly one guide per cell and apply the first original `passed_qc`
  flag when present. The supplied control metadata has no such flag.
- Use historical stored feature counts. Gene eligibility is original metadata
  `n_nonzero / retained_cell_count >= 0.005`; the numerator is not recomputed
  after cell filtering. Original guide counts of at least 10 select eligible
  target groups and NT exposure guides. All original assignments remain
  available, and all guides within an eligible target group are pooled.
- Normalize feature IDs from underscores to dashes.
- Use `log(n_umis)` and `log(n_nonzero)` from original full-transcriptome
  metadata, without `+1` or batch adjustment.

On the supplied control data this retains 30,486 cells, including 3,770 NT
controls, and 15,449 response genes. There are 215 eligible measured target
groups and 73 eligible NT exposure guides.

PCIT retains its five score calibrations, 50 positive pairs and 950 null pairs,
with two-sided adaptive resampling using `h=15` and at most 25,000 draws per
pair. The positive targets are the 50 eligible genes with the highest original
stored mean expression. The fixed null response-gene pool excludes all
annotated targeted genes and requires mean expression at least 0.005. Within
each replicate, response genes are sampled with replacement and paired with
three distinct random NT guides. The paper's separate effective-sample-size
threshold of seven is not applied. This aligns preprocessing while retaining
our comparison design.

## Checks and bundled results

From the repository root:

```bash
./run_all_real_data.sh check
./run_all_real_data.sh validate-only
./run_all_real_data.sh plot-only
Rscript --vanilla check_outputs.R --verify-results
```

`check` runs small numerical, preprocessing, I/O and result-contract checks.
`validate-only` checks replicate IDs 1–20, their preprocessing/configuration,
p-values, rejection sets, and execution-record structure without writing files. `plot-only`
validates the saved results, aggregates the stored decisions and creates the
q=0.10 figures. Plotting requires no raw inputs or numerical rerun. Adaptive
stopped p-values are specific to q=0.10; they are not passed through BH again.

The bundled numerical results are the selected completed q=0.10 run, adopted
with every numerical result object preserved exactly. Their provenance is
explicitly marked `imported_verified_q10`, with the original result bytes,
execution receipts and source snapshots retained for verification. Historical
code/input metadata is recorded separately; adoption is not reported as a
fresh execution. These fingerprints describe the recorded run; validation
does not require the present source files or raw inputs to match them.
Compact records are in the result folder's `provenance/` directory.

## Complete numerical runs

Obtain original processed inputs with the checksum-verified `prepare_data.py`
route in the [package README](../README.md#frangieh-real-data-analysis).
`LOCAL_FRANGIEH_2021_DATA_DIR` must contain the original
`processed/control/{gene,grna_assignment}` tree before single-guide filtering.
Choose a fresh result root (or a result root whose `frangieh-2021/` directory
is empty):

```bash
export LOCAL_FRANGIEH_2021_DATA_DIR=/absolute/path/frangieh-2021
export LOCAL_PDRT_DATA_DIR=/absolute/path/fresh-results
./run_all_real_data.sh full --local
```

Every `full` invocation requires a fresh `frangieh-2021/` directory, claims it
before computation, and runs every selected replicate from 1 through `n_reps`.
It does not resume, skip or replace previous work. A permanent
`.full-run-started` marker also prevents restarting an interrupted invocation.
Use a new result root for another attempt; use `plot-only` for saved results.

Local mode runs the selected replicates sequentially and aggregates/plots when
all 20 are complete. `full 1 --local` is a smaller execution trial and cannot
produce canonical 20-replicate figures. The default `full` uses SGE; after jobs
finish, run `plot-only`. `PROJECT` selects an optional scheduler project;
`MEM` defaults to `16G` through `m_mem_free`. Paths, R and the selected library
must be accessible on each node.

Paths below are relative to `LOCAL_PDRT_DATA_DIR`:

```text
frangieh-2021/intermediate-data-gb2024-single-guide/
frangieh-2021/score-calibration-gb2024-q10/
frangieh-2021/plots/
figures-and-tables/
```

Explicit path overrides remain supported. During `full`, the intermediate,
replicate-result and diagnostic-plot directories must lie beneath the claimed
fresh `frangieh-2021/` directory, without parent traversal. Read-only validation
and plotting can still select saved results elsewhere. The random-number layout uses four
fixed blocks per replicate, with the selected run's coordinate retained for
q=0.10. Paired permutation calibrations share draws, as do paired sign-flip
calibrations. Removing unused FDR levels does not change these streams.

## Publication and environment

The dispatcher submits all selected replicates once and refuses an existing
result directory containing outputs or its permanent `.replicates-dispatched`
claim. Each worker also refuses an existing target file and publishes its
new result atomically. There are no resume, missing-job or force-replacement
options.

Preprocessing requires a fresh intermediate directory and publishes its single
immutable generation beneath `.generations/`, selected by
`current-generation.rds`. It never reuses or replaces an existing bundle.
The generation records its preprocessing contract, cohort audit, target counts
and provenance. Saved ODMs can be relocated using
`LOCAL_FRANGIEH_2021_DATA_DIR`; the loader still validates physical dataset
identity, original QC/count metadata, normalized IDs, subsets and axis order
before rebinding the backing path.

New replicate files record the preprocessing contract, input/code fingerprints,
package versions, R version and session for inspection. Saved-result validation
checks scientific configuration and decisions, with internal archive integrity
checks for the bundled selected run. It does not compare historical source or
input hashes with today's files, or require execution metadata to match across
replicates.

Run `./setup_environment.sh` once to restore packages. Runners restore only
with `PCIT_RESTORE=1`. Personal configuration is opt-in through
`PCIT_USE_RESEARCH_CONFIG=1` or `PCIT_CONFIG_FILE`; explicit paths take precedence.
