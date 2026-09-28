## caustests 1.1.4

This release corrects the computations below and replaces the 1.1.3 submission, which should be discarded.

* `xtpcmg()`: the one-sided long-run covariance used in the fully modified bias correction was transposed (it estimated the sum of E(u_t v_{t+j}) instead of E(v_t u_{t+j})), and the quadratic spectral and Daniell kernels did not use the one-sided weights of the authors' code; group-mean and pooled FM-OLS estimates were therefore biased. The long-run covariance now follows the authors' lr_varmod.m.
* `xtpcmg()`, pooled model: the covariance matrix is now the asymptotic covariance of de Jong and Wagner (2022) for one-way and two-way effects (as in the authors' PanelEKC code), with a heteroskedasticity-robust sandwich for the controls; the previous version used sigma^2 (X'X)^-1 from the FM residuals, and a unit matrix when X'X was singular.
* `xtpcmg()`, cross-section robust covariance (`corr_rob = TRUE`): uses the conditional long-run covariance between units instead of the covariance of u alone.
* All `xtpcmg()` estimates and standard errors now reproduce the Stata module xtpcmg 1.0.2; tests added.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
