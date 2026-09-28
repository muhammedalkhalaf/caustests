## caustests 1.1.3

This release corrects the computations below; the 1.1.2 submission (reference metadata only) should be discarded in favour of this one.

* Bug fix: in the Toda-Yamamoto type tests (1, 3, 5, 6 and 7) the lag order was selected on the VAR augmented with the extra `dmax` lags, so the selected p was that of the augmented model. The lag order p is now selected on the VAR in levels with p lags, and the test VAR then has p + dmax lags, as in Toda and Yamamoto (1995).
* Bug fix: the Fourier terms were evaluated on a time index that restarted at 1 after the lags were dropped, so sin(2 pi k t / T) was shifted and scaled by the effective sample size. They now use the time index t = 1, ..., T of the full sample, in the selection step as in the test regression.
* With these changes the Wald statistics, lags and frequencies of tests 1 to 5 agree with the Stata command caustests (SSC) on the same data (for example 9.326 for test 1 and 21.944 for test 2).

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
