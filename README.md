# Reproduction code

This directory contains the analysis and plotting code for the numerical
results in `project-progress/combine/combine.tex`.  The organization follows a
shallow manuscript-reproduction layout: one directory per analysis domain,
plotting code beside that domain, and three entry points at the repository
root.

## Entry points

From the repository root, run

```bash
./run_all_simulation.sh
./run_all_semi_synthetic.sh
./run_all_real_data.sh
```

With no argument, the simulation and semi-synthetic runners perform the
complete analysis and then generate their figures.  The real-data runner
submits the 20 Perturb-seq replicates to an SGE cluster; after those jobs have
finished, use its `plot-only` mode to validate, aggregate, and plot the
results.  To regenerate figures from existing or downloaded result files, use

```bash
./run_all_simulation.sh plot-only
./run_all_semi_synthetic.sh plot-only
./run_all_real_data.sh plot-only
```

All runners activate the repository-level `renv` environment.  Numerical
results are stored outside the source tree under `LOCAL_PDRT_DATA_DIR` (or the
documented local fallback), while final PDFs are written to
`project-progress/combine/figures-and-tables`.

## Contents

- `simulation-code/`: GLM variance robustness, negative-binomial score
  calibration, continuous-regression GCM/PCM comparisons, and regularization
  sensitivity experiments.
- `semi-synthetic-analysis/`: the 30-target CCLE synthetic-response analysis.
- `real-data-analysis/`: the Frangieh Perturb-seq FDR and power analysis.
- `plotting-style.R`: the shared visual encoding used by the GCM/PCM and CCLE
  figures.
- `check_outputs.R`: verifies that all 11 active manuscript PDFs exist and are
  nonempty.

Exploratory analyses not used by the active manuscript are intentionally
excluded.

## Data and results

The simulations require no external raw data.  The CCLE workflow expects
`mutation.txt` and `response.csv` under
`$LOCAL_PDRT_DATA_DIR/ccle`.  The Frangieh workflow expects the processed
control-condition `ondisc` inputs documented in
`real-data-analysis/data_preprocessing.R`. Set
`LOCAL_FRANGIEH_2021_DATA_DIR` to the directory containing the `processed/`
input tree; generated intermediate files and results remain under
`LOCAL_PDRT_DATA_DIR`.

The Frangieh input tree must contain these four files:

```text
processed/control/gene/gene_expression_matrix.odm
processed/control/gene/gene_expression_metadata.rds
processed/control/grna_assignment/grna_assignments_ungrouped.odm
processed/control/grna_assignment/grna_assignments_ungrouped_metadata.rds
```

The complete simulation workflow requires Nextflow and pins the
`katsevich-lab/simulatr-pipeline` workflow to commit
`dd2f1fcdb533a7764f05dd9dd05b175b02b264a8`. It uses Nextflow's `standard`
profile by default; set `PCIT_NEXTFLOW_PROFILE=local` for an interactive run.
The complete CCLE workflow runs locally by default. It first runs a small,
data-independent matched-calibration smoke test, then runs all 60 setups,
strictly combines them, and regenerates the figures. The domain-level
`submit-ccle.sh` supports asynchronous SGE submission for the full rerun; see
the CCLE instructions below. The complete Frangieh workflow requires an SGE
cluster and the processed `ondisc` data. The `plot-only` workflows require
only the corresponding saved result files. Supplying fewer than 20 Frangieh
replicates is useful for a cluster test, but the postprocessor deliberately
requires replicate IDs 1--20 before it will overwrite manuscript figures.

## CCLE matched-calibration v2

The corrected CCLE analysis uses the result namespace
`ccle-typeI-synthY-30targets-matched-v2`. Within each GCM or PCM pair, normal
and permutation calibration are applied to the same observed statistic. They
therefore share the same synthetic response, nuisance fits, and, for PCM, the
same sample split and learned direction; the pair differs only in its
reference distribution. This is a calibration comparison, not a comparison
of two independently fitted test statistics.

The runner assigns deterministic, purpose-specific random-number streams to
the synthetic response, PCM split, permutation stream, and learner roles.
Consequently, adaptive permutation stopping or changing method execution
order does not alter later responses, splits, or nuisance fits. Each per-setup
RDS records the schema and provenance needed by the combiner to verify this
contract.

The v2 filenames deliberately differ from the legacy unmatched filenames.
The v2 combiner and plotter accept only the matched-v2 schema, so old results
cannot be silently mixed into the corrected figures. Per-setup results and
the combined result are written under `$LOCAL_PDRT_DATA_DIR/ccle`; the
combined file is
`ccle-typeI-synthY-30targets-matched-v2-results.rds`.

For a local full rerun, use

```bash
./run_all_semi_synthetic.sh full
```

For an SGE rerun, first restore the repository environment, run the smoke
test, and then submit one job per setup:

```bash
RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript -e 'renv::restore(prompt = FALSE)'
Rscript reproduction-code/semi-synthetic-analysis/smoke-test-matched-calibration.R
bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh
```

The submission script automatically sources `~/.research_config`, resolves
the absolute `Rscript` path, and embeds only that executable and
`LOCAL_PDRT_DATA_DIR` in each job script. Each worker starts in the project
root, where the project's `.Rprofile` activates `renv`. The launcher
deliberately avoids `qsub -V` and does not forward the full login environment.
The configured paths must be visible at the same location on every compute
node. `PROJECT` and `MEM` optionally override the SGE project and per-job
memory request. After all 60 jobs have finished successfully, run

```bash
./run_all_semi_synthetic.sh plot-only
```

If the matched-v2 combined RDS is absent, `plot-only` strictly validates and
combines all 60 per-setup files before plotting. If it is already present,
the plotter validates and uses it. Legacy combined results are ignored.

The corrected Frangieh workflow uses schema
`perturb_score5_v2_nonpositive_nulls`; legacy v1 results are not accepted as
inputs to the new postprocessing code.

## Verification

After generating all figures, run

```bash
Rscript --vanilla reproduction-code/check_outputs.R
```

The development `pcmt` dependency is pinned in `renv.lock` to commit
`b92c1ccb4f8982fb63d35d7b0986cb3fc5af13bd` of
`ZiangNiu6/pcmt`. That repository must be publicly readable (or made
available to reviewers) for a clean `renv::restore()`.
