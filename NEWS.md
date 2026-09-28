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
  - Test 2: Single Fourier Granger (Enders & Jones, 2016)
  - Test 3: Single Fourier Toda-Yamamoto (Nazlioglu et al., 2016)
  - Test 4: Cumulative Fourier Granger (Enders & Jones, 2019)
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
