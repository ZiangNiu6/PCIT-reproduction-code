# Ren2021 covariate-adjusted calibration

`ren2021-covariate-calibration.R` analyzes the DCATS-packaged major-cell-type
counts from Ren et al. (2021), using the control versus severe/critical
progression contrast (20 versus 48 PBMC libraries).

The analysis follows the requested PCIT specification:

- the null outcome model adjusts for age, sex, and fresh/frozen PBMC status;
- those covariates enter only the `Y | Z` fit;
- disease labels are randomized marginally with the 20/48 allocation fixed;
- there is no `X | Z` model, propensity fit, exposure residualization, or
  restricted/stratified permutation;
- the statistic uses the beta-binomial proportional-variance residual
  `(Y - M * mu) / {1 + (M - 1) * phi}` and empirical studentization;
- DCATS-style pooled- and cell-specific-dispersion LRTs retain their native
  one-degree-of-freedom chi-square p-values. They are not permutation
  calibrated.

The script uses a label-blind null dispersion for the controlled DCATS
comparators. This differs from the published Figure 4 pipeline, which estimated
a pooled dispersion from the full design (including state) and applied a
synthetic 10% similarity correction. The present comparison isolates
small-effective-sample calibration on the same raw annotated counts.

## Input

The compact `Ren2021.rda` object is distributed with DCATS. Pass its path or
install DCATS. For a DCATS source checkout at `/tmp/DCATS`, run:

```bash
Rscript --vanilla \
  reproduction-code/real-data-analysis/ren2021-covariate-calibration.R \
  --ren-data=/tmp/DCATS/data/Ren2021.rda \
  --outer-draws=1000 \
  --reference-draws=19999 \
  --cores=6
```

The production seed is `20260822`. The script validates a deterministic subset
of its direct fixed-dispersion LRTs against `aod` when that package is
available. Use `--validate-aod=false` to skip only that numerical audit.

## Production result

The production run used 1,000 outer marginal allocations and 19,999 independent
reference allocations. At nominal alpha 0.05, the results across the ten cell
types were:

| Method | Mean rejection | Cell-type range | Mean absolute error | Maximum absolute error |
|---|---:|---:|---:|---:|
| PCIT permuted proportional-variance score | 0.0516 | 0.042–0.059 | 0.0044 | 0.009 |
| Same score, normal reference | 0.0619 | 0.052–0.080 | 0.0119 | 0.030 |
| DCATS LRT, pooled null phi | 0.0511 | 0.000–0.226 | 0.0545 | 0.176 |
| DCATS LRT, cell-specific null phi | 0.0460 | 0.013–0.072 | 0.0142 | 0.037 |

The pooled-phi DCATS average is therefore misleading: monocytes reject in
22.6% of draws and CD8 T cells in 15.1%, while dendritic cells and
megakaryocytes never reject. Cell-specific dispersion substantially improves
DCATS, but remains strongly conservative for megakaryocytes (1.3%) and plasma
cells (1.6%). The PCIT range is compatible with Monte Carlo variation from
1,000 outer draws (the binomial standard error near 0.05 is about 0.0069).

The score-information effective sample size ranges from 2.69 to 32.26 despite
68 nominal libraries. It is lowest for megakaryocytes (2.69), macrophages
(3.18), neutrophils (3.27), and plasma cells (6.28). The fitted pooled null
dispersion is 0.07052. All fits returned finite p-values. Nine direct LRTs were
checked against `aod`; the maximum absolute likelihood-ratio difference was
6.85e-5.

## Interpretation

This is a conditional marginal-randomization diagnostic on a fixed response
and covariate table, not a repeated-cohort simulation. Age was not randomized,
and marginal label permutations deliberately break the observed association
between disease state and the covariates. The CI justification therefore comes
from the correctly specified, outcome-only `Y | Z` score plus studentization,
not from design exchangeability or a conditional randomization law.

DCATS treats PBMC library rows as independent. In the source metadata, the 68
libraries in this contrast correspond to 59 patients; repeated libraries are
not represented in the packaged design and are not modeled here. A donor-aware
analysis using the original sample workbook should be reported as a separate
sensitivity analysis.

## Outputs

- `data/ren2021/ren2021-covariate-calibration.csv`: per-cell-type calibration
  at four nominal thresholds.
- `data/ren2021/ren2021-covariate-calibration-alpha05.csv`: the primary 0.05
  table.
- `data/ren2021/ren2021-covariate-observed.csv`: observed-label association
  p-values and BH-adjusted q-values.
- `data/ren2021/ren2021-covariate-calibration.rds`: allocation matrices,
  fitted null objects, statistics, p-value matrices, diagnostics, and session
  metadata.
- `project-progress/combine/figures-and-tables/ren2021-covariate-calibration.*`:
  calibration plot.
