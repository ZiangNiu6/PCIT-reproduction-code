# Reproduction code

This standalone folder reproduces the numerical figures for **Conditional
independence testing through marginal permutation**. Like
[SKSD-reproduction-code](https://github.com/ZiangNiu6/SKSD-reproduction-code), it
contains its environment, entry points, saved results, and analysis/figure code.
Commands below run from this directory; prefix paths with `reproduction-code/`
when running from the parent repository. Scripts also work from another working
directory.

## Environment

Use **R 4.4.1**, a Unix-like shell, and an internet connection for initial package
installation. Source-package installation can require a compiler and HDF5
development libraries for `ondisc`/`Rhdf5lib`. Python 3 is needed only for the
optional data-preparation script. Nextflow and Java are needed only to compute
missing NB/continuous-regression simulation results; SGE is optional.

```bash
./setup_environment.sh
```

`renv.lock` pins the R/Bioconductor packages; `R-packages.txt` is a readable
inventory. The exact `pcmt` source revision is included in `renv/cellar/` so it
does not require reviewer GitHub credentials. Runners never restore packages or
contact package repositories unless explicitly invoked with `PCIT_RESTORE=1`.

## Quick reproduction from bundled results

The 24 canonical result files under `results/` are copies of the existing
publication results, including the combined CCLE v4 artifact and all 20
Frangieh replicate files. `results/artifact-manifest.json` records byte counts
and SHA-256 hashes. These results predate the new provenance fields; they have
not been relabeled as new computations.

```bash
./run_all_simulations.sh plot-only
./run_all_semi_synthetic.sh plot-only
./run_all_real_data.sh plot-only
Rscript --vanilla check_outputs.R --verify-results
```

These commands require no raw data, Nextflow, or cluster. They create ten PDFs
in `results/figures-and-tables`. Diagnostic plots and tables remain beside
their domain results. Figure checking validates PDF headers, the exact figure
set, and nonempty files. `--verify-results` additionally checks the bundled
numerical artifacts against the manifest; omit it for fresh numerical reruns,
whose bytes may differ across systems. Use `--figure-dir DIR` or
`--results-root DIR` for another output location.

Small, data-independent checks do not start publication-scale analyses:

```bash
./run_all_simulations.sh check
./run_all_semi_synthetic.sh check
./run_all_real_data.sh check
```

To validate saved numerical results without writing figures or summaries, use
`validate-only` on any of the three entry points.

## Complete numerical reruns

`full` validates existing outputs and computes only missing work. For a fresh
run, choose a new `LOCAL_PDRT_DATA_DIR`; keep the bundled archive intact.
Set `PCIT_FIGURE_DIR` to override the default `<results-root>/figures-and-tables`.
Relative runner environment paths are resolved against this package, not the
calling directory. Explicit Python/checker CLI paths use the calling directory.

### Simulations

```bash
LOCAL_PDRT_DATA_DIR=/absolute/path/fresh-results ./run_all_simulations.sh full
```

Only the active studies run: the illustrative GLM variance example, the
continuous-dispersion NB study with estimated working size, and the matched
GCM/PCM comparison with dCRT. The Nextflow workflow is pinned to
`katsevich-lab/simulatr-pipeline@dd2f1fcdb533a7764f05dd9dd05b175b02b264a8`.
The default profile is `local`; `PCIT_NEXTFLOW_PROFILE=standard` selects the
pipeline's cluster profile. Nextflow/Java installation is a separate system
dependency. Complete cached results are validated and can be reused without
Nextflow installed.

The GLM study has 8,000 replicate rows. The NB validator requires 30 grid cells
× 5 methods × 2,000 replicates (300,000 outputs). The continuous-regression
validator requires 27 cells × 1,000 replicates from each of two outer methods
(matched GCM/PCM and dCRT), expanding to five displayed procedures. Normal and
permutation members of each matched pair share the observed statistic and fits.

The default resource requests are 4 GB and 1 hour per NB job, and 4 GB and
2 hours per continuous-regression job; these are per-job limits, not total
runtime estimates. Override `PCIT_SIM_MAX_GB`, `PCIT_SIM_MAX_HOURS`,
`PCIT_PCM_SIM_MAX_HOURS`, and `PCIT_GLM_CORES` as needed. Publication-scale work
can take many hours or longer; total wall time depends on hardware and dispatch.
Malformed/stale cached artifacts stop the run with recovery instructions.

### CCLE semi-synthetic analysis

The two public source files have been checked against the existing analysis
inputs: the response CSV is identical, and the mutation input is identical
after removing exactly two GCT header lines. Download and verify them with:

```bash
python3 prepare_data.py ccle --destination /absolute/path/fresh-results/ccle
python3 prepare_data.py ccle --destination /absolute/path/fresh-results/ccle --check
LOCAL_PDRT_DATA_DIR=/absolute/path/fresh-results ./run_all_semi_synthetic.sh full
```

The inputs come from Broad's
[binary mutation/CNA calls](https://data.broadinstitute.org/ccle_legacy_data/binary_calls_for_copy_number_and_mutation_data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct)
and [2015-02-24 drug-response table](https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv).
`prepare_data.py` pins the original download SHA-256 values, checks the GCT
version/dimensions, records the transformation and prepared-file hashes, and
refuses to replace incomplete or modified inputs. `--source-dir DIR` prepares
already downloaded original files offline. Allow roughly 200 MB for downloading
and staging these approximately 100 MB of inputs.

The production v4 design has 30 null settings plus four nonzero effects for
each of EIF5A, PLK1, and ATP11A: 42 settings × 500 replicates, with 500-tree
forests and up to 499 adaptive permutation draws. EIF5A uses effects
`0,1,2,3,4`; PLK1 and ATP11A use `0,0.5,1,1.5,2`, at fixed active-set size 10.
Local execution is sequential. For parallel SGE submission use
`semi-synthetic-analysis/submit-ccle.sh --sge`, then `plot-only` after completion.
Optional `PROJECT` and `MEM` override scheduler settings. Expect a substantial
publication-scale workload; no new end-to-end runtime benchmark is claimed.

Valid setup files are skipped; wrong-design, mixed-runtime, or malformed files
fail closed. `full --force` explicitly replaces the expected setup outputs.
New outputs record input/code hashes and package/session provenance. Plot-only
uses a complete combined artifact directly, without requiring raw inputs or
recombination. `CCLE_RESULTS_SCHEMA=auto|v4|v3|v2` is respected: auto prefers
v4 > v3 > v2, and partial newer namespaces cannot trigger silent fallback.
Legacy schemas are a diagnostic compatibility option, not the current
ten-figure bundle; the current manuscript uses v4.

The current main CCLE composite has two sparsity rows and three columns:
all-ten-target calibration, representative-target calibration, and
conditional-effect rejection. Its displayed
power representatives are EIF5A and ATP11A; the separate effect-grid diagnostic
also includes PLK1. The introduction composite shows EIF5A calibration and
conditional-effect rejection. These are nominal rejection curves, not
size-adjusted power estimates.

### Frangieh real-data analysis

**The existing Frangieh cohort and guide-handling rules are deliberately
retained in this revision.** No single-guide restriction or correction to the
historical `pure_nt` definition has been applied. The saved results therefore
retain the previously discussed multi-guide/control-cohort limitation. Workflow
repairs do not establish that those scientific results would be unchanged after
a future cohort correction.

The processed inputs are published in
[Barry et al.'s Zenodo data archive, version v1](https://zenodo.org/records/10976417).
The Frangieh ZIP has published MD5 `1ad2adc53fde94fa17cad9ac6ebca006`.
This helper downloads the approximately 2.8 GB ZIP, verifies it, and extracts
only the four control-condition input files (about 4.3 GB); allow at least
8 GB of free space during preparation:

```bash
python3 prepare_data.py frangieh --destination /absolute/path/frangieh-2021
export LOCAL_FRANGIEH_2021_DATA_DIR=/absolute/path/frangieh-2021
export LOCAL_PDRT_DATA_DIR=/absolute/path/fresh-results
./run_all_real_data.sh full --local
```

`LOCAL_FRANGIEH_2021_DATA_DIR` must contain:

```text
processed/control/gene/gene_expression_matrix.odm
processed/control/gene/gene_expression_metadata.rds
processed/control/grna_assignment/grna_assignments_ungrouped.odm
processed/control/grna_assignment/grna_assignments_ungrouped_metadata.rds
```

These are the author's processed inputs before this package's cohort selection;
do not substitute SCEPTRE's single-guide-filtered analysis intermediates.
The original raw-to-ondisc import recipe is pinned at
[import-frangieh-2021/de111dbc](https://github.com/Katsevich-Lab/import-frangieh-2021/tree/de111dbc442d63f7537217fa7dcdca8500afeb3c).
The archive route avoids rebuilding the original raw conversion. The complete
multi-GB download was not repeated during this maintenance task; its identity
and expected archive layout were checked against the public record.

Local mode runs missing replicates sequentially and aggregates/plots when all
20 are complete. `full 1 --local` is a small execution trial, but cannot produce
canonical manuscript figures. The default `full` mode uses SGE; after jobs
finish, run `plot-only`. Scheduler memory defaults to 16 GB via `m_mem_free`;
`PROJECT` is optional and defaults to the site's scheduler account. Each
replicate analyzes 1,000 pairs at two FDR levels with up to 25,000 adaptive draws
per pair. See [real-data workflow details](real-data-analysis/README-workflow.md)
for resume/force behavior, atomic intermediate generations, ODM path relocation,
and new-result provenance.

## Layout and configuration

```text
reproduction-code/
  README.md, R-packages.txt, renv.lock, renv/, .Rprofile
  setup_environment.sh, prepare_data.py, check_outputs.R
  run_all_simulations.sh, run_all_semi_synthetic.sh, run_all_real_data.sh
  simulation-code/           active specifications, method checks, figure-code/
  semi-synthetic-analysis/   active CCLE analysis, validation, figure-code/
  real-data-analysis/        active Frangieh analysis, validation, figure-code/
  results/                   canonical numerical artifacts and generated figures
```

The package does not load personal configuration by default.
`PCIT_USE_RESEARCH_CONFIG=1` opts into `~/.research_config`, or use
`PCIT_CONFIG_FILE=/absolute/path/config`. Explicit path variables take
precedence. `PCIT_R_LIBRARY` can select an already installed library. Within
the full research checkout, the legacy parent `renv` library is also recognized;
an exported standalone folder uses its own restored library.

Removed regularization figures, exploratory searches, pilots, and alternative
data analyses live under the parent repository's `development-code/`. They are
not part of this standalone publication package. The repository-root runners
are compatibility wrappers; new workflows should use the entry points here.
