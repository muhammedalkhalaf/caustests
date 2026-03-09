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
