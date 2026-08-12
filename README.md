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
The complete CCLE
workflow runs locally by default; the domain-level `submit-ccle.sh` also
supports SGE submission.  The complete Frangieh workflow requires an SGE
cluster and the processed `ondisc` data.  The `plot-only` workflows require
only the corresponding saved result files. Supplying fewer than 20 Frangieh
replicates is useful for a cluster test, but the postprocessor deliberately
requires replicate IDs 1--20 before it will overwrite manuscript figures.

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
