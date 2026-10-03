## caustests 1.1.5

This release changes how a quantreg warning is reported; no numerical results change. The version on CRAN is 1.1.4.

* Tests 6 and 7 (quantile causality): `quantreg::summary.rq(se = "nid")` warns "k non-positive fis" when the sparsity estimate is non-positive at some observations. The package called it once for the observed sample and once per bootstrap replication at every quantile, so the same warning was repeated up to (nboot + 1) times per quantile and direction (for example 965 times for test 7 on `caustests_data` with 5 quantiles and `nboot = 99`), which users mistook for an error. These warnings are now counted inside the loops and reported once per call, with the directions, quantiles and number of fits affected. Other warnings are not suppressed.
* New element `quantreg_warnings` in the returned object (tests 6-7): a data frame with the number of affected quantile regression fits per direction and quantile. `print()` adds one line when any occurred.
* No numerical results change: Wald statistics and bootstrap p-values are identical to version 1.1.4 for the same seed.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
