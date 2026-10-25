# Constant-size sparse NB Type-I search

The requested `manuscript-update/` directory is absent in this checkout. This
search uses the sparse conditional-independence row in
`manuscript/combine.tex` and its formal implementation in
`sim-spec-formal/nb-permutation-signflip.R`.

## Fixed model and comparisons

For the original null datasets, \(n=400\), \(Z\sim N(0,1)\),
\(X\mid Z\sim\mathrm{Bernoulli}\{\operatorname{expit}(a+\delta Z)\}\), and
\(Y\mid X,Z\sim\mathrm{NB}\{10\exp(\beta_Z Z),\theta\}\). The logistic
intercept \(a\) is recalibrated at each \(\delta\) to keep marginal
\(\Pr(X=1)=0.025\). The exposure effect on the outcome is always
\(\gamma=0\), so conditional independence holds. The NB size \(\theta\)
is **constant across observations**; its inverse \(1/\theta\) is the
dispersion parameter. For Barry and our tests, the working null model is
`MASS::glm.nb(Y ~ Z)`, with a common size estimated from each dataset.

Because the NB working variance family is correctly specified for Barry and
our tests, the manuscript's large-sample result for information
studentization applies to that comparison. Any failure of Barry's NB-based
test found here would be a finite-sample effect of the sparse design.

The initial search evaluated Barry et al.'s information-studentized
permutation test and our empirically studentized permutation test. They use
the same dataset, fitted null model, and raw-label permutations. Both
p-values use the two-sided plus-one rule, and rejection means
\(p\leq0.05\). A failed fit or calibration counts as nonrejection in the
all-dataset rate; outputs also give rates among successful fits.

## Search and independent confirmation

The broad screen covered all 144 combinations of
\(\theta\in\{0.05,0.1,0.25,0.5,1,2\}\),
\(\beta_Z\in\{-1.5,-0.75,0,0.5,1,1.5\}\), and
\(\delta\in\{0.25,0.5,1,2\}\). An additional 80-cell screen used
\(\theta\in\{0.1,0.25,0.5,1,2\}\),
\(\beta_Z\in\{-3,-2,2,3\}\), and \(\delta\in\{1,2,3,4\}\).
Both screens used 100 datasets per cell and 399 permutations. Their largest
individual rejection rates are exploratory: the cells were selected from
224 noisy estimates. In the extreme screen, 45 of 80 cells had more than 5%
fitting/calibration failures.

Five promising cells were rerun with **new seeds**, 2,000 datasets, and
1,999 permutations, matching the manuscript's Monte Carlo settings. The
original sparse setting was rerun as a reference. Exact binomial 95% Monte
Carlo intervals are shown for the all-dataset rejection rates.

| NB size \(\theta\) | Outcome \(\beta_Z\) | Exposure \(\delta\) | Barry rejections / 2,000 (95% CI) | Ours rejections / 2,000 (95% CI) | Failures |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 0.5 | 0.5 | 88; 4.40% (3.54–5.39%) | 94; 4.70% (3.82–5.72%) | 0 |
| 0.1 | 1 | 0.5 | 108; 5.40% (4.45–6.48%) | 63; 3.15% (2.43–4.01%) | 11 |
| 0.25 | −1.5 | 1 | 95; 4.75% (3.86–5.78%) | 72; 3.60% (2.83–4.51%) | 2 |
| 1 | 1 | 2 | 108; 5.40% (4.45–6.48%) | 91; 4.55% (3.68–5.56%) | 0 |
| 0.5 | 2 | 1 | 103; 5.15% (4.22–6.21%) | 70; 3.50% (2.74–4.40%) | 33 |
| 0.5 | −2 | 2 | 83; 4.15% (3.32–5.12%) | 89; 4.45% (3.59–5.45%) | 20 |

**Conclusion:** The constant-size NB search did not establish a setting in
which Barry et al. fails to control Type-I error at 5%. Its largest
independently confirmed rate was 5.4%, and every Barry interval includes
5%. The gap between the methods in some cells reflects our test's
conservatism as well as any possible Barry inflation, so the gap alone does
not show failure of nominal control. These finite-grid results do not prove
uniform control over all parameters.

## Poisson-only SCEPTRE fitting step

The search script now includes `poisson_only` as a third test. It exactly
matches the first response-fitting step in SCEPTRE's
[`perform_response_precomputation`](https://github.com/Katsevich-Lab/sceptre/blob/main/R/precomputation_functions.R):
`stats::glm.fit(Y ~ Z, family = poisson())` on the design matrix. It then
uses a Poisson working score, with \(W=\hat\mu\), score contributions
\(\tilde X(Y-\hat\mu)\), and Barry's information-studentized, two-sided
plus-one raw-label permutation calibration. **It never estimates or uses NB
size.** This is a Poisson-only variant of the test, not SCEPTRE's full NB
method; the latter also estimates size from the Poisson fit.

The Poisson fit and p-value have separate success flags, so they remain
available when `MASS::glm.nb` fails. The same stored permutation seed pairs
all three tests. At the small, constant size \(\theta=0.1\),
\(\beta_Z=1\), and \(\delta=0.5\), two 2,000-dataset runs with 1,999
permutations yielded:

| Seed | Barry | Ours | Poisson-only | NB failures | Poisson failures |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20261011 | 108/2,000 (5.40%) | 63/2,000 (3.15%) | 129/2,000 (6.45%) | 11 | 0 |
| 20261020 | 106/2,000 (5.30%) | 63/2,000 (3.15%) | 137/2,000 (6.85%) | 11 | 0 |

Pooling the two independent runs gives 266/4,000 = **6.65%** for the
Poisson-only test, with an exact binomial 95% Monte Carlo interval of
**5.90–7.47%**. Thus this Poisson working-score variant overrejects at this
selected setting. It does not establish a Type-I failure of Barry's NB-size
test. The first run reproduces the previous Barry and our p-values exactly.

## Pitfall of a standard-normal reference for our score

The search script now also reports `ours_normal`. This method uses **the
identical fitted NB model and observed empirical-studentized score** as our
permutation test, but assigns the two-sided p-value
\(p_{\mathrm{normal}}=2\Phi(-|T_{\mathrm{ours}}|)\) instead of comparing
\(T_{\mathrm{ours}}\) to studentized label permutations. It is distinct from
the Poisson-only method above. Degenerate scores get \(p=1\); failed NB fits
remain nonrejections in the all-dataset rates.

Applying this deterministic transformation to the saved null statistics gives:

| Constant-size setting | Barry permutation | Our permutation | Our normal reference |
| --- | ---: | ---: | ---: |
| Baseline: \(\theta=0.5,\ \beta_Z=0.5,\ \delta=0.5\); 2,000 datasets | 88 (4.40%) | 94 (4.70%) | 237 (11.85%) |
| Matched effects: \(\theta=0.5,\ \beta_Z=1,\ \delta=0.5\); 2,000 fresh datasets | 89 (4.45%) | 94 (4.70%) | 222 (11.10%) |
| Small size: \(\theta=0.1,\ \beta_Z=1,\ \delta=0.5\); 4,000 datasets over two seeds | 214 (5.35%) | 126 (3.15%) | 1,073 (26.83%) |

Holding both covariate effects fixed at \((1,0.5)\), the normal reference
overrejects even at size \(0.5\): 11.10% (exact 95% Monte Carlo interval
9.76–12.56%), versus 4.70% for our permutation test and 4.45% for Barry.
The problem becomes more severe at size \(0.1\).

For the small-size setting, the exact binomial 95% Monte Carlo interval for
the normal-reference rate is **25.46–28.23%**. Of its 1,073 rejections,
1,071 have a negative observed score. The observed score's 2.5th and 97.5th
percentiles are **−2.93 and 1.22**, far from the symmetric standard-normal
cutoffs \(\pm1.96\). The normal approximation therefore treats common
negative deviations as significant. The permutation comparison rejected
only 126 datasets; 951 were rejected by the normal reference alone.
These are finite-sample findings under sparse exposure, not a failure of the
score's asymptotic normal limit.

The diagnostic plot and exact counts are saved as
`results/simulation-results/nb-sparse-ci-search/normal-reference-diagnostic.png`
and `normal-reference-confirmations.csv`. The plot holds both covariate
effects fixed while varying only NB size. Recreate both files with
`plot-nb-normal-reference.R`.

## Larger sample with 0.5% exposure

To check a rarer-exposure design at larger \(n\), we fixed the constant NB
size at \(\theta=0.5\), the outcome covariate effect at \(\beta_Z=1\), and
the exposure covariate effect at \(\delta=0.5\). We changed only
\(n=10{,}000\) and the calibrated marginal exposure prevalence to 0.005.
The outcome exposure effect remained zero. Each of two independent seeds
generated 2,000 null datasets. Barry, our permutation test, and Poisson-only
used 1,999 two-sided plus-one permutations per dataset; the normal-reference
p-value used the standard-normal tail. All NB and Poisson fits succeeded.

| Seed | Barry | Ours, permutation | Ours, normal reference | Poisson-only |
| ---: | ---: | ---: | ---: | ---: |
| 20261026 | 125/2,000 (6.25%) | 95 (4.75%) | 151 (7.55%) | 226 (11.30%) |
| 20261027 | 93/2,000 (4.65%) | 100 (5.00%) | 145 (7.25%) | 217 (10.85%) |
| **Pooled** | **218/4,000 (5.45%)** | **195 (4.88%)** | **296 (7.40%)** | **443 (11.08%)** |

The exact 95% Monte Carlo intervals for the pooled rates are **4.77–6.20%**
for Barry, **4.23–5.59%** for our permutation test, **6.61–8.26%** for our
normal-reference test, and **10.12–12.09%** for Poisson-only. The first
Barry run's 6.25% rate was above nominal, but the independent run returned
4.65%; the pooled interval includes 5%. This design therefore does **not**
provide convincing evidence that Barry fails to control Type-I error.
Across all 4,000 datasets, Barry alone rejected 153 and ours alone rejected
130; their paired rejection-rate difference (Barry minus ours) was 0.58
percentage points (approximate 95% interval −0.25 to 1.40 points; exact
discordance \(p=0.191\)).

The average exposed count was 49.77 (SD 6.91), versus about 10 at
\(n=400\) and 2.5% exposure. The average number of exposed observations with
a nonzero outcome was 40.38. The observed empirical-score effective sample
size averaged only 8.43 (median 7.33). The empirical-studentized score's
2.5th and 97.5th percentiles were −2.79 and 1.64. Of the 296
normal-reference rejections, 273 had a negative score. Thus the normal
reference remains anti-conservative even
with roughly 50 exposed observations. The Poisson-only working score is
more anti-conservative still. Because \(n\) and prevalence changed together
from the earlier run, neither change can be assigned the entire difference
in rejection rates.

The two dataset-level RDS bundles and CSV summaries are saved under
`results/simulation-results/nb-sparse-ci-search/` in the directories
`confirm-n10k-p0p005-theta0p5-beta1-delta0p5` and
`replicate-n10k-p0p005-theta0p5-beta1-delta0p5`.
For a new run of this design, use
`--nb-sizes=0.5 --beta-z=1 --x-z-slopes=0.5 --n=10000
--x-prevalence=0.005 --reps=2000 --permutations=1999
--permutation-engine=sparse` with a fresh `--seed` and `--output-dir`.

## Fitting-time comparison at small NB size

In the installed SCEPTRE implementation, response precomputation fits a
Poisson GLM with `stats::glm.fit`, estimates NB size from the fitted Poisson
means in compiled code, and retains the Poisson coefficients. It does not
alternate the mean and size fits as `MASS::glm.nb` does. The benchmark uses
the same NB datasets at true constant size \(\theta=0.1\), with
\(Z\sim N(0,1)\) and \(\mu=10\exp(0.5Z)\). Times below are mean milliseconds
per response fit on this machine; they exclude permutation tests.

| \(n\) | SCEPTRE Poisson-only | SCEPTRE Poisson + NB size | `MASS::glm.nb` |
| ---: | ---: | ---: | ---: |
| 400 | 1.39 | 1.77 | 6.83 |
| 10,000 | 25.02 | 26.70 | 105.58 |
| 100,000 | 127.55 | 147.80 | 665.50 |

The full six-size timing output is in
`results/simulation-results/nb-fit-scaling/benchmark-sceptre-theta0p1.csv`,
and the reproducible benchmark is `benchmark-sceptre-nb-fit.R`. The SCEPTRE
Poisson-only timing is its exact first `glm.fit` step. SCEPTRE's internal
routine may change with package versions; this run used installed version
0.10.2.

## Reproduction

Run `search-nb-sparse-ci-null.R` from the repository root using
`Rscript --vanilla` (plain `Rscript` can wait on this project's `renv`
sandbox lock). The script accepts comma-separated `--nb-sizes`, `--beta-z`,
and `--x-z-slopes`, plus `--n`, `--x-prevalence`, `--reps`, `--permutations`,
`--permutation-engine`, `--cores`, `--seed`, and `--output-dir`. The `sparse`
engine samples exposed subsets directly and preserves the exact binary-label
permutation distribution. Its projection formulas and rejection counts were
checked against the original dense implementation on identical sampled
subsets; the largest statistic difference was \(8.0\times10^{-15}\), and all
three methods had identical permutation hit counts. It writes the grid,
dataset-level p-values, method summaries,
Barry-versus-ours paired summaries, and an RDS bundle. The method summaries
include method-specific failure counts for the Poisson-only and normal-reference
tests. The saved
outputs are in
`results/simulation-results/nb-sparse-ci-search/`; that directory is ignored
by Git. The paired permutation code was checked against both formal method
functions on 16 smoke datasets, with exact p-value agreement.
