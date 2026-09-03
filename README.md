# Reproduction code

This directory contains the analysis and plotting code for the numerical
results in `manuscript/combine.tex`.  The organization follows a
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
documented local fallback), while the ten manuscript-facing PDFs are written
to `manuscript/figures-and-tables`. Additional diagnostic plots remain beside
their analysis results rather than being copied into the manuscript bundle.

## Contents

- `simulation-code/`: GLM variance robustness, negative-binomial score
  calibration, continuous-regression GCM/PCM comparisons, and regularization
  sensitivity experiments.
- `semi-synthetic-analysis/`: the 30-target CCLE synthetic-response analysis.
- `real-data-analysis/`: the Frangieh Perturb-seq FDR and power analysis.
- `plotting-style.R`: the shared visual encoding used by the GCM/PCM and CCLE
  figures.
- `check_outputs.R`: verifies the eight baseline/fallback manuscript PDFs and,
  when v4 CCLE composites are present, requires both current composite figures;
  every checked PDF must exist and be nonempty, and unreferenced image files
  are rejected from the manuscript figure directory.

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

The Section 5.1 negative-binomial study writes to
`simulation-results/nb-permutation-signflip-continuous-dispersion-estimated-theta`.
In its variance-misspecification scenario, `Z` remains standard normal and the
true NB size is `theta(Z) = 46 exp(-2.1 Z)`, while `glm.nb` estimates one common
working size from each dataset. This common-size working model remains
misspecified because the true size varies with `Z`. The plotter strictly checks
this contract and requires matching
`continuous_nb_dispersion_estimated_theta` schema markers from both the grid
and each method output, so fixed-size or copied legacy results fail validation.
It does not read the former `nb-permutation-signflip-continuous-dispersion` or
`nb-permutation-signflip` result directories. Output PDF filenames are
unchanged so the manuscript includes remain stable. The five
signal magnitudes for this row are `0, 0.055, 0.11, 0.165, 0.22` under MI and
`0, 0.0675, 0.135, 0.2025, 0.27` under CI.

The Section 5.2 continuous-regression study writes to the isolated
`simulation-results/pcm-additive-matched-dcrt` namespace. Its formal simulatr
specifier registers two outer methods: `matched_gcm_pcm_rf` and `dcrt_rf`.
Within each GCM or PCM pair, the normal and permutation p-values are computed
from the same observed statistic, random-forest fits, and (for PCM) sample
split. The `dcrt_rf` entry is specifically the in-sample learned-law dCRT
analyzed by Niu et al., using a Gaussian working distribution with RF mean
fits. It uses the same 500-tree nuisance regressions as GCM and the
homoskedastic working law
`X | Z ~ N(mu_X_hat(Z), sigma_X_hat^2)`, with 499 fixed conditional draws and
an independently generated draw stream for each dataset. It reports the
plus-one, absolute-statistic two-sided p-value; 499 matches the permutation cap
used by the other methods and is a study-specific computational choice. The
Gaussian law is correctly shaped in the first two rows but intentionally
misspecified in the heteroscedastic, shifted-exponential row. The
full runner first executes
`simulation-code/smoke-test-pcm-additive-matched.R`; the result validator then
requires 27,000 complete outputs from each outer method. The manuscript figure
displays the four matched GCM/PCM procedures together with dCRT as a fifth
curve. Its legends identify dCRT as the test and conditional resampling as its
calibration. dCRT remains separately identifiable in the result object, and
the plot workflow also writes its cellwise rejection rates to
`pcm-additive-dcrt-rejection-summary.csv` beside that result. Pass
`--summary-only` to `simulation-code/plotting-code/plot-pcm-additive.R` to
validate both outer methods and write this table without regenerating the PDF.
Legacy `pcm-additive` and pre-dCRT `pcm-additive-matched` results are never
accepted by this plotter. The full runner budgets two hours per Section 5.2
job by default; set `PCIT_PCM_SIM_MAX_HOURS` to override that request.

The complete CCLE workflow runs locally by default. It first runs a small,
data-independent matched-calibration smoke test, then resumes the 30 null and
12 nonnull representative-target settings at `s = 10` (42 jobs total),
strictly combines them, and regenerates the figures. The domain-level
`submit-ccle.sh` supports the same missing-only resume behavior for asynchronous
SGE submission; see the CCLE instructions below. The complete Frangieh workflow requires an SGE
cluster and the processed `ondisc` data. The `plot-only` workflows require
only the corresponding saved result files. Supplying fewer than 20 Frangieh
replicates is useful for a cluster test, but the postprocessor deliberately
requires replicate IDs 1--20 before it will overwrite manuscript figures.

## CCLE representative effect-grid matched calibration (v4)

The CCLE analysis uses the result namespace
`ccle-typeI-synthY-30targets-matched-v4`. Within each GCM or PCM pair, normal
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

The v4 design fixes the synthetic-response active-set size at `s = 10` and
retains one `delta = 0` null setting for each of the 30 mutation targets used
in the paper's calibration analysis. Nonzero effects are run only for the
three conditional-effect representatives: `EIF5A_MUT` (high sparsity/rare),
`PLK1_MUT` (medium sparsity), and `ATP11A_MUT` (low sparsity/common). Their
paired grids are

```text
Y^(delta) = mu(Z) + epsilon + delta * sigma * X_j
EIF5A_MUT:  delta in {0, 1, 2, 3, 4}
PLK1_MUT:   delta in {0, 0.5, 1, 1.5, 2}
ATP11A_MUT: delta in {0, 0.5, 1, 1.5, 2}
```

The four-carrier target `EIF5A_MUT` is the high-sparsity representative in
both the paper's null-calibration and p-value-quantile figures and the v4
conditional-effect diagnostic.

Thus each nonzero effect is a `delta`-noise-SD per-carrier shift and differs
from that target's paired null only through the added target effect: the active
set, coefficients, Gaussian noise, PCM split, permutation streams, and learner
streams are shared across its five grid points. The binary target is not
variance-standardized, avoiding implausibly large effects for one-carrier
mutations. The 30 `delta = 0` jobs deliberately reuse the old matched-v2
`s = 10` random-number streams, so existing and new null results can be checked
exactly. The complete contract therefore contains 42 records: 30 null records
plus four additional effect records for each of the three representatives.

The v4 filenames deliberately differ from legacy unmatched, matched-v2, and
matched-v3 filenames. Per-setup files use
`ccle-typeI-synthY-30targets-matched-v4-results-{TARGET}-{EFFECT_ID}-s10.rds`,
where `EFFECT_ID` is the target-appropriate member of `delta0`, `delta0p5`,
`delta1`, `delta1p5`, `delta2`, `delta3`, or `delta4`. The v4 combiner accepts
only the exact 30-null/12-alternative paired schema, so older results cannot be
silently mixed into the effect-grid analysis. Per-setup and combined results
are written under `$LOCAL_PDRT_DATA_DIR/ccle`; the combined files are
`ccle-typeI-synthY-30targets-matched-v4-results.rds` and
`ccle-typeI-synthY-30targets-matched-v4-rejection-rates.csv`.

For figure-only reproduction during the transition, set
`CCLE_RESULTS_SCHEMA=auto|v4|v3|v2`. The default `auto` uses the first available
artifact in the strict priority `v4 > v3 > v2`; an explicit value requires that
exact artifact. V3 produces its single-effect transitional diagnostic, while
v2 retains only its 30 `s = 10` null records and skips the effect diagnostic.
The plotter never accepts unmatched legacy results.

For a local full rerun, use

```bash
./run_all_semi_synthetic.sh full
```

This command is resumable: before dispatch, the combiner's strict per-record
validator checks every expected output. Valid setup files are skipped and only
missing setup indices are run. Malformed, stale, heterogeneous, or wrong-setup
files fail closed and are left unchanged; obsolete unexpected matched-v4 files
must be moved out of the result namespace. To deliberately replace all 42
reviewed expected outputs, use

```bash
./run_all_semi_synthetic.sh full --force
```

For an SGE rerun, first restore the repository environment, run the smoke
test, and then submit jobs for the missing setup indices:

```bash
RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript -e 'renv::restore(prompt = FALSE)'
Rscript reproduction-code/semi-synthetic-analysis/smoke-test-matched-calibration.R
bash reproduction-code/semi-synthetic-analysis/submit-ccle.sh
```

The same command can be rerun safely after an interrupted array: strictly valid
outputs are skipped. `--force` (or `CCLE_FORCE_RERUN=1`) explicitly submits all
selected indices and permits their workers to replace reviewed expected files;
it does not authorize an obsolete unexpected filename, which must first be
moved out of the namespace.

The submission script sources `~/.research_config` only when
`LOCAL_PDRT_DATA_DIR` is not already set, resolves the absolute `Rscript` path,
and embeds only that executable and
`LOCAL_PDRT_DATA_DIR` in each job script. Each worker starts in the project
root, where the project's `.Rprofile` activates `renv`. The launcher
deliberately avoids `qsub -V` and does not forward the full login environment.
The configured paths must be visible at the same location on every compute
node. `PROJECT` and `MEM` optionally override the SGE project and per-job
memory request. If all selected outputs are already valid, SGE mode exits
successfully without requiring `qsub`. After all 42 jobs have finished
successfully, run

```bash
./run_all_semi_synthetic.sh plot-only
```

If the matched-v4 combined RDS is absent but any matched-v4 per-setup output
exists, `plot-only` strictly requires and combines all 42 setup files; it never
silently falls back after a partial v4 run. With a complete v4 artifact, the
plotter validates the exact 30 nulls, the three target-specific five-point
grids, and within-target pairing. It uses only `delta = 0` for the unchanged
null-calibration figures and writes
`ccle-alternative-synthY-30targets-rejection.pdf`. That GCM-only diagnostic
shows separate EIF5A, PLK1, and ATP11A curves on target-specific x-axis ranges:
0 to 4 for EIF5A and 0 to 2 for PLK1 and ATP11A. The panel labels state each
target's actual grid. It displays
raw rejection at the nominal `0.05` cutoff over `delta`, including the paired
null. Because the methods can have different null rejection rates, this is not
a size-adjusted power curve. V4 also writes two manuscript composites:
`ccle-typeI-synthY-intro-calibration-power.pdf` places the EIF5A null bar plot
beside its conditional-effect curve, and
`ccle-typeI-synthY-30targets-typeI-power-combined.pdf` adds the full
three-target conditional-effect plot as row (c) beneath the two null rows.
The manuscript automatically uses these composites when they exist and keeps
the null-only legacy figures as a transition fallback. To prevent stale panels,
the plotter fails closed if a v2/v3 result is selected while either v4
manuscript composite remains in the figure directory. Unmatched legacy
combined results are ignored.

## CCLE asymptotic-GCM effect-grid pilot

The standalone effect-grid pilot is for choosing a production alternative;
its estimates are selection-only and are not accepted by the matched-v4
combiner or plotter. It keeps the 30 paper targets and `s = 10`, uses the
independent seed base `30260000`, and evaluates the paired per-carrier effects
`0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2`. Within every target and replicate it
generates the outcome components once, caches one `X | Z` random forest across
all eight effects, and fits only `Y | Z` separately for each effect. Each target
artifact includes an exact comparison of the cached public-`ranger` statistic
against `pcmt::perm_gcm(B = 1)` from its initial run.

The eight-effect runner uses schema `ccle_asym_gcm_effect_pilot_v3`. If its
target-specific RDS/CSV pair is a strictly valid prefix artifact—either the
six-effect v1 grid ending at `1.25` or the seven-effect v2 grid ending at
`1.5`—it retains those rows, regenerates and checks the paired outcome
components and cached `X | Z` fit for every replicate, strictly revalidates the
original preflight, computes only the missing suffix, and transactionally
upgrades the same files. A valid
eight-effect v3 pair is skipped. Missing, mismatched, or malformed artifact
pairs fail closed and are not overwritten; a target with no artifacts runs all
eight effects.

The runner accepts one target, comma-separated targets, inclusive `:` or `-`
ranges, or `all`. For example, two long-lived local processes avoid reloading
the CCLE inputs for every target:

```bash
source ~/.research_config
export PILOT_R=25 PILOT_RF_TREES=500 PILOT_RF_THREADS=4
Rscript reproduction-code/semi-synthetic-analysis/pilot-ccle-asym-gcm-effect-grid.R 1:15 &
Rscript reproduction-code/semi-synthetic-analysis/pilot-ccle-asym-gcm-effect-grid.R 16:30 &
wait
Rscript reproduction-code/semi-synthetic-analysis/combine-pilot-ccle-asym-gcm-effect-grid.R
```

`PILOT_OUTPUT_DIR` optionally overrides the default
`$LOCAL_PDRT_DATA_DIR/ccle/pilot-asym-gcm-grid`. The runner writes one atomic
RDS/CSV pair per target, with filenames recording target index, mutation,
replicate count, and tree count. After all 30 matching target files exist, the
combiner transactionally publishes a combined RDS, an equal-target summary
CSV, and a target-level summary CSV named from `R` and the tree count. The
summaries contain nominal rejection rates, paired gains over effect zero, and Monte
Carlo standard errors; the combined RDS also retains validated raw rows and
full method, seed, preflight, and package provenance.

The pilot considered all eight effect sizes independently of the production
replicates. For target screening, the four-carrier `EIF5A_MUT` was also checked
on `delta = 0, 1, 2, 3, 4` with `R = 25`. Its asymptotic-GCM nominal rejection
rates were `0, 0.08, 0.16, 0.24, 0.40`; a local adaptive-permutation GCM check
gave `0.04, 0.12, 0.32, 0.48, 0.56`. These estimates have wide Monte Carlo
uncertainty and are used only to screen the illustrative target and grid. The
choice must be confirmed by the independent `R = 500` production run, and no
pilot replicates are pooled with it. The production grids are `0, 1, 2, 3, 4`
for `EIF5A_MUT` and `0, 0.5, 1, 1.5, 2` for `PLK1_MUT` and `ATP11A_MUT`;
nonzero effects are not run for the other 27 targets. Production plots show
both asymptotic and permutation GCM; neither screening check is presented as
confirmatory power evidence.

For auditability, the `R = 25` selection-only equal-target asymptotic-GCM
nominal rejection rates were:

| Displayed sparsity group | 0 | 0.25 | 0.5 | 0.75 | 1 | 1.25 | 1.5 | 2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| High sparsity (rare) | .004 | .008 | .028 | .032 | .044 | .052 | .072 | .120 |
| Medium sparsity | .040 | .080 | .172 | .284 | .388 | .572 | .704 | .864 |
| Low sparsity (common) | .044 | .104 | .252 | .468 | .668 | .812 | .928 | .976 |

The stored pilot tags describe mutation-frequency terciles, so the plot and
table intentionally invert the extreme labels: internal `low_sparsity` means
rare mutations and is displayed as high sparsity, while internal
`high_sparsity` means common mutations and is displayed as low sparsity.

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
