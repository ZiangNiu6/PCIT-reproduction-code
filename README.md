# Reproduction code

This standalone repository reproduces the numerical figures for **Conditional
independence testing through marginal permutation**. Like
[SKSD-reproduction-code](https://github.com/ZiangNiu6/SKSD-reproduction-code), it
contains its environment, entry points, saved results, and analysis/figure code.
Commands below run from the repository root. Scripts also work from another
working directory.

## Environment

Use **R 4.4.1**, a Unix-like shell, and an internet connection for initial package
installation. Source-package installation can require a compiler and HDF5
development libraries for `ondisc`/`Rhdf5lib`. Python 3 is needed only for the
optional data-preparation script. Nextflow and Java are needed only to compute
missing NB/continuous-regression simulation results; SGE is optional.

```bash
./setup_environment.sh
```

On macOS, source builds of `ondisc` also need OpenSSL. With Homebrew, install
`openssl@3` and `pkg-config` before setup. The setup script discovers the
library directory through `pkg-config` or Homebrew; set `PCIT_OPENSSL_LIB_DIR`
to the directory containing `libcrypto` for a custom installation. Personal
R Makevars files are not required for this library-path configuration.

`renv.lock` pins the R/Bioconductor packages; `R-packages.txt` is a readable
inventory. Setup installs [`pcmt`](https://github.com/ZiangNiu6/pcmt) from its
public GitHub repository at commit
`92c17ecc2359cfd4cc02c306bef793e1fbd4c77e`; no GitHub credentials are required.
Runners never restore packages or contact package repositories unless
explicitly invoked with `PCIT_RESTORE=1`.

## Quick reproduction from bundled results

The 24 canonical result files under `results/` bundle the publication
results, including the combined CCLE v4 artifact and all 20
single-guide Frangieh replicate files at BH **q=0.10**.
`results/artifact-manifest.json` records byte counts and SHA-256 hashes.
The Frangieh files preserve the selected completed run's numerical results
exactly and record their verified adoption separately from a fresh execution.

```bash
./run_all_simulations.sh plot-only
./run_all_semi_synthetic.sh plot-only
./run_all_real_data.sh plot-only
Rscript --vanilla check_outputs.R --verify-results
```

These commands require no raw data, Nextflow, or cluster. They create ten PDFs
in `results/figures-and-tables`. Only the 24 canonical numerical artifacts,
their manifest, compact Frangieh execution provenance, and those ten PDFs are
bundled under `results/`. Derived summaries, diagnostic plots, and previews are
rebuilt beside their domain results and ignored by Git. Figure checking
validates PDF headers, the exact figure set, and nonempty files.
`--verify-results` additionally checks the bundled
numerical artifacts against the manifest; omit it for fresh numerical reruns,
whose bytes may differ across systems. Use `--figure-dir DIR` or
`--results-root DIR` for another output location.

Small, data-independent checks do not start publication-scale analyses:

```bash
./run_all_simulations.sh check
./run_all_semi_synthetic.sh check
./run_all_real_data.sh check
bash tests/fresh-only-workflow.sh
```

To validate saved numerical results without writing figures or summaries, use
`validate-only` on any of the three entry points.

## Complete numerical reruns

`full` always computes a new study. Choose a fresh `LOCAL_PDRT_DATA_DIR`;
the selected study's output directory must contain no previous run files.
CCLE's prepared raw-input files are allowed before its first run.
The bundled manuscript results remain available through `plot-only` and
`validate-only` and are never reused by `full`.

A `.full-run-started/` marker permanently claims each numerical output
directory. Completed, concurrent and interrupted runs all require a different
directory for another `full` invocation. There is no automatic resume or
`--force` replacement. Partial outputs remain available for inspection.
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
dependency. `plot-only` and `validate-only` use saved results without Nextflow.

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
Existing simulation outputs stop `full` before numerical work starts.
Each simulation specification is rebuilt from its current source for a new run.

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
Local execution is sequential. For parallel SGE submission to a fresh CCLE
output directory use `semi-synthetic-analysis/submit-ccle.sh --sge`, explicitly
combine the complete 42-file grid, then run `plot-only`.
Optional `PROJECT` and `MEM` override scheduler settings. Expect a substantial
publication-scale workload; no new end-to-end runtime benchmark is claimed.

All 42 settings are computed afresh. Existing setup or combined outputs are
refused rather than skipped or replaced. New outputs retain input/code and
package/session metadata for inspection; source hashes are not used to resume
older work. Plot-only reads the complete current v4 combined artifact directly,
without raw inputs or automatic recombination. Older v2/v3 schemas and their
selection/fallback paths have been removed.

The current main CCLE composite has two sparsity rows and three columns:
all-ten-target calibration, representative-target calibration, and
conditional-effect rejection. Its displayed
power representatives are EIF5A and ATP11A; the separate effect-grid diagnostic
also includes PLK1. The introduction composite shows EIF5A calibration and
conditional-effect rejection. These are nominal rejection curves, not
size-adjusted power estimates.

### Frangieh real-data analysis

**The bundled results and new runs use the literal preprocessing of Barry et al.'s 2024 Genome
Biology implementation**, pinned to manuscript commit
[`005c28b`](https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/005c28b9ba30e7d39c3ee18d51111722f49d32ec/R_scripts/uniform_processing_lowmoi.R).
They retain exactly one original assigned guide per cell before feature QC,
apply the original author's `passed_qc` flag when present, normalize feature
IDs from underscores to dashes, and use `log(n_umis)` and `log(n_nonzero)`
from the original full-transcriptome cell metadata as covariates. There is no
`+1` transformation or batch covariate.

Feature QC deliberately follows the historical stored-count implementation:
gene eligibility is original metadata `n_nonzero / retained_cell_count >= 0.005`,
without recomputing the numerator after cell filtering. The original guide
count threshold of 10 selects eligible target groups and NT exposure guides;
all guide assignments are retained for the initial cell QC, and all guides
within an eligible target group remain pooled. On the local control data this
gives 30,486 retained cells, 3,770 NT controls, 15,449 eligible genes, 215 eligible
measured target groups, and 73 eligible NT guides.

This aligns preprocessing, while retaining **PCIT's five-method, 50-positive /
950-null comparison**, its seeds and resampling settings. The 50 positive
targets are ranked by original stored mean expression after historical
eligibility checks; null candidates exclude all annotated targeted genes.
The paper's separate effective-sample-size threshold of seven is not applied,
and this is not a reproduction of its complete benchmarking protocol.

The active analysis uses only **BH q=0.10**. Superseded result bundles,
their plots, and the historical schema selector have been removed. The
selected completed q=0.10 run is retained exactly, with a compact archive of
its original execution receipts and source snapshots. Validation does not
depend on the removed result files or on access to raw data for plotting.

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

Local mode runs new replicates sequentially and aggregates/plots when all
20 are complete. `full 1 --local` is a small execution trial, but cannot produce
canonical manuscript figures. The default `full` mode uses SGE; after jobs
finish, run `plot-only`. Scheduler memory defaults to 16 GB via `m_mem_free`;
`PROJECT` is optional and defaults to the site's scheduler account. Each
replicate analyzes 1,000 pairs at `q = 0.10`, retaining `h = 15`
and at most 25,000 adaptive draws per pair. All tables and headline figures
use that level and the saved adaptive rejection sets. Across the 20 bundled
replicates, the empirically studentized permutation method has mean power
75.3% and mean false discovery proportion 9.61%.
See [real-data workflow details](real-data-analysis/README-workflow.md)
for fresh-run behavior, intermediate publication, ODM path relocation,
and execution metadata.

Under the selected result root, the workflow uses these paths:

```text
frangieh-2021/intermediate-data-gb2024-single-guide/
frangieh-2021/score-calibration-gb2024-q10/
frangieh-2021/plots/
figures-and-tables/
```

`check_outputs.R --verify-results` checks the ten bundled figure files and
all 24 canonical numerical artifacts, including the selected Frangieh run.

## Layout and configuration

```text
PCIT-reproduction-code/
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
precedence. The repository discovers its restored `renv` library;
`PCIT_R_LIBRARY` can select another already installed library.

Development analyses and historical verification receipts remain in
[PCIT-project](https://github.com/ZiangNiu6/PCIT-project). This repository
preserves the history of its `reproduction-code/` folder through project commit
[`892167f`](https://github.com/ZiangNiu6/PCIT-project/commit/892167fd2159c5861c895ffe43386586f095cc1e).
The research repository's compatibility runners use this sibling checkout.
