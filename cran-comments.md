## caustests 1.1.7

CRAN version: 1.1.6.

### Changes

* xtpcaus(test = "pfty"): lag orders leaving fewer than 5 residual degrees of freedom are excluded (warning), error when none is admissible.
* xtpcaus(test = "pfty"): Dumitrescu and Hurlin (2012) Z-bar (eq. 9) and Z-tilde (eq. 26, moments of eq. 21 and 22) with one-sided p-values, reported for the plain VAR(p) specification (kmax = 0, dmax = 0) only; NA for the Fourier or Toda-Yamamoto augmented specification.
* xtpcaus(test = "pfty"): bootstrap resamples restricted-model residuals with the same time indices for all units, regenerates the series recursively and repeats the lag and frequency selection; bootstrap p-value of Z-tilde added.
* xtpcaus(test = "pfty"): singular unit regressions give NA with a warning and are excluded from the panel statistics.
* xtpcaus(): missing values and duplicated or misaligned (panel_id, time_id) rows stop with an informative error.
* xtpcaus(): seed applied in a local RNG scope; the global .Random.seed is restored.
* kmax = 0 and dmax = 0 accepted by xtpcaus().
* Printed labels: "to" instead of the arrow symbol, "and" instead of the ampersand between author names.
* Reference to Dumitrescu and Hurlin (2012) added; Details section added to ?xtpcaus.
* New tests (tests/testthat/test-xtpcaus-dh.R); validation ledger in inst/validation.

### Test environments

* Linux (x86_64), R 4.3.3. R-devel was not available in the build environment; the package was checked on R 4.3.3 only.

### R CMD check results

R CMD check --as-cran --no-manual (R 4.3.3, _R_CHECK_SYSTEM_CLOCK_=0): 0 errors, 0 warnings, 1 NOTE.

* NOTE (environment): checking CRAN incoming feasibility reports the URL https://github.com/muhammedalkhalaf/caustests as "Status: 403 Forbidden". The check ran behind a proxy that blocks the request; the repository is public.

Examples run with --run-donttest (28 s). Tests: tests/testthat (all pass).
