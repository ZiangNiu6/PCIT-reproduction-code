# Frangieh reproduction workflow

The standalone package runner uses the existing five-method v2 experiment.
The user explicitly retained its current scientific cohort and guide handling:
there is no single-guide filter, and the historical `pure_nt` mask excludes
only guides whose targets occur in the expression feature table. This differs
from Barry et al.'s single-guide analysis. These workflow changes do not correct
that scientific issue or change the saved numerical results, seeds, positive
and null pair construction, or calibration methods.

From the `reproduction-code` directory:

```bash
./run_all_real_data.sh check
./run_all_real_data.sh validate-only
./run_all_real_data.sh plot-only
```

`validate-only` checks every canonical replicate ID 1–20, configuration, p-value
and decision structure, and provenance when present, without writing files.
`plot-only` aggregates those saved rejection sets and writes figures beneath
`results/figures-and-tables`. Neither command restores packages automatically.
Run `./setup_environment.sh` first when dependencies are missing, or explicitly
set `PCIT_RESTORE=1` on a runner invocation.

For numerical reruns, obtain the processed Frangieh data from the original
SCEPTRE data archive at [Zenodo](https://doi.org/10.5281/zenodo.10976417), with
the import route documented in the original
[import-frangieh-2021 repository](https://github.com/Katsevich-Lab/import-frangieh-2021).
Set `LOCAL_FRANGIEH_2021_DATA_DIR` to the directory containing:

```text
processed/control/gene/gene_expression_matrix.odm
processed/control/gene/gene_expression_metadata.rds
processed/control/grna_assignment/grna_assignments_ungrouped.odm
processed/control/grna_assignment/grna_assignments_ungrouped_metadata.rds
```

Use a separate result directory for a fresh run:

```bash
export LOCAL_FRANGIEH_2021_DATA_DIR=/absolute/path/frangieh-2021
export LOCAL_PDRT_DATA_DIR=/absolute/path/new-results
./run_all_real_data.sh full --local
```

Local mode runs missing replicates sequentially and, for the complete 20-replicate
run, aggregates and plots on completion. A smaller `full 1 --local` run is a
cluster-independent execution check; it cannot overwrite canonical figures.
The default `full` mode uses SGE. Set `PROJECT` and `MEM` for the site's project
and memory resource. With `PROJECT` unset, SGE uses the site default project;
the memory default is `16G` using `m_mem_free`. Paths and the selected R library must be accessible on each node.
After SGE jobs finish, run `plot-only`.

Both modes validate all existing expected files before dispatch and run only
missing replicates. Existing complete v2 artifacts remain accepted even though
they predate the new provenance fields. Invalid, noncanonical, or stale outputs
fail closed; unexpected replicate filenames are always rejected. To deliberately
replace the selected expected replicate outputs, add `--force`. That flag does
not remove unexpected files or silently replace invalid intermediate bundles.

Complete flat intermediate directories from the original workflow remain
readable. Newly generated intermediates are written as immutable generations
under `.generations/`; an atomic `current-generation.rds` pointer publishes the
whole bundle. The production and exploratory readers resolve that pointer.
`data_preprocessing.R 1000 --force` explicitly publishes a replacement generation
while preserving the old one. Saved ODMs copied from another machine can be
relocated by supplying `LOCAL_FRANGIEH_2021_DATA_DIR`: the loader verifies the
saved ODM identity and preserves all saved subsets and metadata while rebinding
only the backing path.

New replicate files are published atomically and record intermediate hashes,
analysis-code hashes, package versions, R version, and session information.
Existing v2 artifacts are not rewritten to manufacture missing provenance.
`PCIT_USE_RESEARCH_CONFIG=1` opts into `~/.research_config`; `PCIT_CONFIG_FILE`
selects another config. Explicit environment paths take precedence.

Run the data-independent numerical smoke check and the temporary-file I/O checks:

```bash
./run_all_real_data.sh check
```
