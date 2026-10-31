# caustests 1.1.7

## Breaking changes

* `xtpcaus()`, PFTY test: `zbar` and `zbar_pv` are `NA` unless `kmax = 0` and `dmax = 0` (plain VAR(p) specification). For that specification they are the Dumitrescu and Hurlin (2012) Z-bar of eq. (9) computed from the individual Wald statistics with the unbiased residual variance, and `zbar_pv` is one-sided (upper tail). The 1.1.6 Z-bar (asymptotic moments applied to the Fourier Toda-Yamamoto Wald statistics, two-sided p-value) is no longer computed.
* `xtpcaus()`, PFTY test: lag orders whose unrestricted regression leaves fewer than 5 residual degrees of freedom are not considered; a warning reports the bound, and the function stops when no lag order is admissible. Selected lags and Wald statistics change on short panels where version 1.1.6 selected such a lag order.
* `xtpcaus()`, PFTY test: bootstrap p-values change for the same seed. The bootstrap now resamples the residuals of the restricted model (tested lags of x removed), draws the same time indices for all units (cross-sectional block of size N, Dumitrescu and Hurlin 2012, Section 6.2), generates each unit's bootstrap series recursively from the restricted model and repeats the lag and frequency selection on it; version 1.1.6 resampled the unrestricted residuals independently per unit on a fixed design with the selected lag and frequency held fixed.
* `xtpcaus()`, PFTY test: a unit whose regressor matrix is singular has `ind_wald`, `ind_pval_b` set to `NA` with a warning and is excluded from the Fisher, W-bar and Z statistics (`fisher_df = 2 * N_used`); version 1.1.6 set its Wald statistic to 0 with p-value 1.
* `xtpcaus()`, PFTY test: `kmax = 0` (no Fourier terms) is now a specification of its own. In 1.1.6 `kmax = 0` was accepted without error, but the frequency loop never ran, so every unit was fitted with the lag order p = 1 and the frequency f = 1 (a Fourier regression with no lag selection). Results with `kmax = 0` change accordingly; for `grunfeld_panel()` with `pmax = 2`, `dmax = 0` the selected lag orders change from 1, 1, 1, 1, 1 to 1, 2, 2, 2, 2, the selected frequencies from 1 to 0, and `W-bar` from 1.3603 to 4.2288.

## Other changes

* `xtpcaus()`, PFTY test: new elements `ztilde`, `ztilde_pv` (Dumitrescu and Hurlin 2012, eq. 26 with the moments of eq. 21 and 22, one-sided p-value), `ztilde_boot_pv` (one-sided cross-sectional block bootstrap p-value), `wbar_dh`, `ind_wald_dh`, `ind_df`, `N_used`. The Z-tilde elements are `NA` unless `kmax = 0` and `dmax = 0`.
* `xtpcaus()`: `dmax = 0` omits the Toda-Yamamoto augmentation.
* `xtpcaus()`: missing values in `y`, `x`, `panel_id` or `time_id` stop with the count and row positions.
* `xtpcaus()`: every unit must be observed at the same set of time values exactly once; duplicated or misaligned rows stop with an error (version 1.1.6 compared only the row count with N times T).
* `xtpcaus()`: `seed` is applied in a local RNG scope; the global `.Random.seed` is restored on exit.
* `print.xtpcaus()`: PFTY output shows the degrees of freedom per unit, the Dumitrescu and Hurlin statistics for the plain specification, and a note instead of Z-bar for the other specifications; the arrow symbol in the PQC labels replaced by "to".
* Documentation: Dumitrescu and Hurlin (2012) added to the references of `xtpcaus()`; Details section describing the specification, the degrees-of-freedom bound, the Z statistics and the bootstrap; "and" instead of "&" in the printed test names of `caustests()` and `xtpcmg()`.
* Tests: hand recomputation of W_i, W-bar, Z-bar and Z-tilde with `lm()` for the plain and the Fourier specifications, the degrees-of-freedom error, the missing-value error, the misaligned-panel error, the singular-unit warning and the RNG scope.
* Validation ledger added as `inst/validation/ledger_xtpcaus.md`.

# caustests 1.1.6

* `xtpcaus()`, PFTY test: the Z-bar statistic is documented as asymptotic in T (sequential limit, T before N). On short panels (for example T = 17) it rejects at p-values near zero even when the bootstrap Fisher statistic does not; the printed output now carries a note recommending the bootstrap Fisher p-value in that case. No change to the computation.

# caustests 1.1.5

* Tests 6 and 7 (quantile causality): `quantreg::summary.rq(se = "nid")` warns "k non-positive fis" when the sparsity estimate is non-positive at some observations. The package called it once for the observed sample and once per bootstrap replication at every quantile, so the same warning was repeated up to (nboot + 1) times per quantile and direction (for example 965 times for test 7 on `caustests_data` with 5 quantiles and `nboot = 99`), which users mistook for an error. These warnings are now counted inside the loops and reported once per call, with the directions, quantiles and number of fits affected. Other warnings are not suppressed.
* New element `quantreg_warnings` in the returned object (tests 6-7): a data frame with the number of affected quantile regression fits per direction and quantile. `print()` adds one line when any occurred.
* No numerical results change: Wald statistics and bootstrap p-values are identical to version 1.1.4 for the same seed.

# caustests 1.1.4

* `xtpcmg()`: the one-sided long-run covariance used in the fully modified bias correction was transposed (it estimated the sum of E(u_t v_{t+j}) instead of E(v_t u_{t+j})), and the quadratic spectral and Daniell kernels did not use the one-sided weights of the authors' code; group-mean and pooled FM-OLS estimates were therefore biased. The long-run covariance now follows the authors' lr_varmod.m.
* `xtpcmg()`, pooled model: the covariance matrix is now the asymptotic covariance of de Jong and Wagner (2022) for one-way and two-way effects (as in the authors' PanelEKC code), with a heteroskedasticity-robust sandwich for the controls; the previous version used sigma^2 (X'X)^-1 from the FM residuals, and a unit matrix when X'X was singular.
* `xtpcmg()`, cross-section robust covariance (`corr_rob = TRUE`): uses the conditional long-run covariance between units instead of the covariance of u alone.
* All `xtpcmg()` estimates and standard errors now reproduce the Stata module xtpcmg 1.0.2; tests added.

# caustests 1.1.3

* Bug fix: in the Toda-Yamamoto type tests (1, 3, 5, 6 and 7) the lag order was selected on the VAR augmented with the extra `dmax` lags, so the selected p was that of the augmented model. The lag order p is now selected on the VAR in levels with p lags, and the test VAR then has p + dmax lags, as in Toda and Yamamoto (1995).
* Bug fix: the Fourier terms were evaluated on a time index that restarted at 1 after the lags were dropped, so sin(2 pi k t / T) was shifted and scaled by the effective sample size. They now use the time index t = 1, ..., T of the full sample, in the selection step as in the test regression.
* With these changes the Wald statistics, lags and frequencies of tests 1 to 5 agree with the Stata command caustests (SSC) on the same data (for example 9.326 for test 1 and 21.944 for test 2).

# caustests 1.1.2

* Corrected the DOI of Wang and Nguyen (2022) to 10.1080/1331677X.2021.1948436.
* Removed a DOI attached to de Jong and Wagner (2022) that could not be verified in CrossRef; the citation text is unchanged.
* No changes to code.

# caustests 1.0.0

## Initial CRAN Release

* Implemented 7 Granger causality tests:
  - Test 1: Toda-Yamamoto (1995)
  - Test 2: Single Fourier Granger (Enders and Jones, 2016)
  - Test 3: Single Fourier Toda-Yamamoto (Nazlioglu et al., 2016)
  - Test 4: Cumulative Fourier Granger (Enders and Jones, 2019)
  - Test 5: Cumulative Fourier Toda-Yamamoto (Nazlioglu et al., 2019)
  - Test 6: Quantile Toda-Yamamoto (Cai et al., 2023)
  - Test 7: Bootstrap Fourier Granger Causality in Quantiles (Cheng et al., 2021)

* Features:
  - Automatic lag order selection via AIC or BIC
  - Optimal Fourier frequency selection
  - Bootstrap inference for robust p-values
  - Support for multivariate systems (all pairwise directions)
  - Quantile causality testing across distribution

* S3 methods: `print()`, `summary()`, `plot()`
* Example dataset: `caustests_data`
