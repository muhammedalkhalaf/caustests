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
