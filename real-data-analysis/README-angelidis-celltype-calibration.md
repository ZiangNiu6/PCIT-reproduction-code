# Angelidis cell-type calibration

This workflow assesses conditional random-label calibration on the observed
mouse-level cell-composition table from Angelidis et al. (Nature
Communications, 2019; GEO GSE124872). It holds the counts and retained-cell
totals fixed and enumerates all `choose(15, 7) = 6,435` allocations of seven
old labels among the 15 mice.

The primary PCIT statistic uses the beta-binomial proportional-variance
quasi-score. For cell type `j` and mouse `i`,

```text
Var(Y_ij | M_i) = M_i mu_j (1-mu_j) {1 + (M_i-1) phi_j}
e_ij = (Y_ij - M_i mu_j) / {1 + (M_i-1) phi_j}.
```

The null-only `phi_j` is estimated without age, and `mu_j` solves the
corresponding intercept quasi-score equation. For every allocation, the test
multiplies `e_ij` by `X_i - 7/15`, empirically studentizes the contributions,
and obtains a two-sided p-value by inclusive rank over the complete orbit.
There are no nuisance covariates in this application, so the conditional
exposure law is exactly the fixed-allocation law and does not require an
estimated propensity model.

DCATS remains on its ordinary chi-square reference; it is not
permutation-calibrated. The controlled primary comparator uses a pooled
beta-binomial dispersion estimated once under a label-blind null model and
then frozen. A cell-specific null-dispersion LRT is reported as a sensitivity
analysis. These are no-similarity, total-denominator DCATS-style likelihood
ratios; the optional DCATS classification-similarity correction is outside
this calibration experiment.

The paper-defined 30 biological identities are primary. The three metadata
labels excluded by the original paper (`low_quality_cells`,
`red_blood_cells`, and `Gamma-Delta_T_cells`) are restored only in the
33-label sensitivity analysis. Denominators are recomputed after this
source-defined filtering: two cells have missing annotations, 14,811 cells
have one of the 33 metadata labels, and 14,124 cells enter the primary
analysis.

## Run

Download the public metadata file:

```bash
curl -L \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124872/suppl/GSE124872_Angelidis_2018_metadata.csv.gz \
  -o /tmp/GSE124872_Angelidis_2018_metadata.csv.gz
```

Then, from the repository root, run:

```bash
Rscript --vanilla \
  reproduction-code/real-data-analysis/angelidis-celltype-calibration.R \
  --metadata=/tmp/GSE124872_Angelidis_2018_metadata.csv.gz \
  --cores=6
```

The implementation uses base R for beta-binomial profile likelihoods. If the
CRAN package `aod` is available, a deterministic subset of fixed-dispersion
likelihood ratios is checked against the implementation used by DCATS.

## Outputs

- `data/angelidis/angelidis-celltype-permutation-calibration.csv`: per-type
  rejection fractions at alpha 0.01, 0.025, 0.05, and 0.10.
- `data/angelidis/angelidis-celltype-permutation-calibration-alpha05.csv`:
  compact alpha-0.05 table.
- `data/angelidis/angelidis-celltype-observed-discoveries.csv`: corrected
  observed-label effects, score ESS, and raw/BH-adjusted p-values for all four
  methods.
- `data/angelidis/angelidis-celltype-bh-calibration.csv`: complete-orbit
  family-level BH calibration at q = 0.05.
- `data/angelidis/angelidis-celltype-permutation-calibration.rds`: complete
  assignment-level statistics, p-value matrices, diagnostics, corrected age
  metadata, input hash, and session information.
- `project-progress/combine/figures-and-tables/angelidis-celltype-calibration.pdf`
  and `.png`: paper-30 per-type comparison.
- `project-progress/combine/figures-and-tables/angelidis-observed-discoveries.pdf`
  and `.png`: signed observed BH evidence for the four methods.

## Interpretation

This is a conditional permutation-orbit diagnostic. Under the hypothetical
law that the seven old labels are exchangeable among these 15 fixed mice, it
measures how often each method rejects over the complete assignment space.
It is not a repeated-data beta-binomial simulation, and age was not randomized
in the original study. Consequently, it does not establish causal effects or
future-study Type-I error under arbitrary sampling and confounding mechanisms.

## Current full-orbit result

For the paper-defined 30 types at nominal alpha 0.05:

| Method | Mean rate | Per-type range | Mean absolute error from 0.05 |
|---|---:|---:|---:|
| PCIT permuted proportional-variance score | 0.04988 | 0.04973--0.04988 | 0.00012 |
| Same score, normal reference | 0.08649 | 0.00699--0.12914 | 0.04011 |
| DCATS-style LRT, pooled null phi | 0.05524 | 0--0.36503 | 0.06435 |
| DCATS-style LRT, cell-specific null phi | 0.04714 | 0.01772--0.07801 | 0.01105 |

The pooled-phi DCATS average is close to 0.05 only because severe liberal and
conservative errors cancel across cell types: its median per-type rate is
0.01740. The largest rates are 0.36503 for Type-2 pneumocytes and 0.31313 for
alveolar macrophages, while ten types never reject over the orbit. The
cell-specific-phi sensitivity is substantially better but remains visibly
less uniform than permutation calibration.

At the family level, applying BH at q = 0.05 to each of the 6,435 complete-null
label assignments gives:

| Method | Probability of any BH rejection | Mean number rejected |
|---|---:|---:|
| PCIT permuted proportional-variance score | 0.04351 | 0.07988 |
| Same score, normal reference | 0.14981 | 0.30179 |
| DCATS-style LRT, pooled null phi | 0.24973 | 0.40171 |
| DCATS-style LRT, cell-specific null phi | 0.01026 | 0.01880 |

For the corrected observed age labels, the four methods make 11, 11, 12, and
9 BH discoveries, respectively. Eight are shared by all four: B cells, CD4 T
cells, classical monocytes, eosinophils, megakaryocytes, non-classical
monocytes, plasma cells, and Type-2 pneumocytes. PCIT additionally detects CD8
T cells, Club cells, and interstitial macrophages. The normal-reference score
swaps CD8 T cells for vascular endothelial cells. The pooled-phi LRT adds
alveolar macrophages and vascular endothelial cells but drops interstitial
macrophages; the cell-specific-phi LRT is a nine-type subset of PCIT.

These disagreements align with the orbit diagnostics but do not identify
biological truth. For example, the pooled-phi LRT's null rejection fraction is
31.3% for alveolar macrophages, whereas its fraction for interstitial
macrophages is 0.062%. The normal-reference score's null rejection fraction is
9.23% for vascular endothelial cells.

The custom fixed-phi likelihood implementation was checked against 35
deterministic `aod::betabin` fits. The maximum absolute likelihood-ratio
difference was `6.21e-6` and the mean absolute difference was `5.87e-7`.
