---
title: 'R Packages for Cointegration Testing with Structural Breaks, Fourier Approximations, and Causality Analysis'
tags:
  - R
  - econometrics
  - cointegration
  - structural breaks
  - Fourier approximation
  - panel data
  - Granger causality
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 18 March 2026
bibliography: paper.bib
---

# Summary

We present eight R packages for cointegration testing and Granger causality analysis: `cointsmall`, `xtcadfcoint`, `xtbreakcoint`, `hatemicoint`, `xtdhcoint`, `fcoint`, `makicoint`, and `caustests`. These packages implement a wide range of cointegration tests—from small-sample methods with as few as twelve observations to large panel tests accommodating cross-sectional dependence, structural breaks, and smooth Fourier approximations—as well as comprehensive Granger causality tests including Fourier, quantile, and bootstrap variants. All packages are open-source under the GPL-3 license.

# Statement of Need

Cointegration analysis is fundamental to empirical economics, underpinning the estimation of long-run equilibrium relationships between non-stationary variables. While the `urca` package [@Pfaff2008] provides classical Engle-Granger and Johansen tests, and the `cointReg` package offers cointegrating regression estimators, the econometric literature has advanced substantially in several directions that remain underserved in R:

1. **Small-sample cointegration**: Standard tests perform poorly with fewer than 50 observations. The @Trinh2022 method handles samples as small as 12 observations with structural breaks.
2. **Panel cointegration with breaks and common factors**: Tests by @Banerjee2015 and @Banerjee2024 accommodate structural breaks, cross-sectional dependence via common factors, and heterogeneous break locations.
3. **Structural break cointegration**: The @HatemiJ2008 and @Maki2012 tests allow for multiple unknown regime shifts in the cointegrating relationship.
4. **Fourier-based cointegration**: Tests by @Banerjee2017fadl and @Tsong2016 use trigonometric terms to capture smooth structural changes without specifying break dates.
5. **Durbin-Hausman panel tests**: @Westerlund2008 provides tests robust to cross-sectional dependence through common factor extraction.
6. **Advanced causality testing**: Beyond standard Granger tests, Fourier [@Enders2016; @Nazlioglu2019], quantile [@Cai2023], and bootstrap Fourier quantile [@Cheng2021] causality tests offer richer causal inference.

# Packages

## cointsmall

Implements the @Trinh2022 cointegration test for very small samples (as few as 12 observations). Applies residual-based ADF-type tests with up to two endogenous structural breaks in the constant or both constant and slope. Includes automatic lag selection via BIC and simulation-based critical values.

```r
library(cointsmall)
result <- cointsmall(y, x, model = "cs", max_breaks = 2)
summary(result)
```

## xtcadfcoint

Implements the @Banerjee2024 panel cointegration test with structural instabilities. Tests the null of no cointegration using cross-sectionally augmented Dickey-Fuller (CADF) regressions with Common Correlated Effects (CCE) estimation, allowing for structural breaks in cointegrating vectors, factor loadings, and deterministic components.

```r
library(xtcadfcoint)
result <- xtcadfcoint(y ~ x1 + x2, data = panel_data,
                      id = "country", time = "year", model = "break")
summary(result)
```

## xtbreakcoint

Implements panel cointegration tests with structural breaks and cross-section dependence following @Banerjee2015. Provides iterative factor-break estimation, individual ADF tests on defactored residuals, standardized panel statistics, and the @BaiNg2004 MQ test for identifying common stochastic trends. Five model specifications with varying deterministic components.

```r
library(xtbreakcoint)
result <- xtbreakcoint(y ~ x1 + x2, data = panel_data,
                       id = "country", time = "year",
                       model = 3, max_breaks = 2)
summary(result)
```

## hatemicoint

Implements the @HatemiJ2008 cointegration test with two unknown regime shifts. Provides ADF*, Zt*, and Zα* test statistics with endogenously determined break dates. Critical values are based on the simulation tables of Hatemi-J (2008).

```r
library(hatemicoint)
result <- hatemicoint(y, x, model = "cs")
summary(result)
```

## xtdhcoint

Implements the @Westerlund2008 Durbin-Hausman panel cointegration tests, which are robust to cross-sectional dependence through common factor extraction via principal components. Provides both group-mean (DHg) and panel (DHp) test statistics with automatic factor number selection.

```r
library(xtdhcoint)
result <- xtdhcoint(y ~ x1 + x2, data = panel_data,
                    id = "country", time = "year")
summary(result)
```

## fcoint

Implements four Fourier-based cointegration tests for smooth structural breaks: the Fourier ADL test (FADL) of @Banerjee2017fadl, the Fourier Engle-Granger tests (FEG and FEG2), and the @Tsong2016 KPSS-type test. Lag and frequency selection via AIC or BIC.

```r
library(fcoint)
result <- fcoint(y, x, test = "fadl", max_freq = 3, ic = "aic")
summary(result)
```

## makicoint

Implements the @Maki2012 cointegration test allowing for an unknown number of structural breaks (up to five). Extends the Gregory-Hansen framework to multiple breaks with four model specifications: level shift, level shift with trend, regime shift, and regime shift with trend.

```r
library(makicoint)
result <- makicoint(y, x, model = "regime_shift",
                    max_breaks = 3, trim = 0.15)
summary(result)
```

## caustests

A comprehensive suite of Granger causality tests: standard Toda-Yamamoto [@TodaYamamoto1995], Fourier single-frequency [@Enders2016] and cumulative-frequency [@Nazlioglu2019] tests, quantile causality [@Cai2023], and bootstrap Fourier Granger causality in quantiles [@Cheng2021]. All tests include bootstrap inference for robust p-values.

```r
library(caustests)
# Fourier Toda-Yamamoto test
result <- ftycaus(y, x, max_lag = 4, max_freq = 3, nboot = 1000)
summary(result)

# Quantile causality test
result <- qcaus(y, x, taus = c(0.25, 0.5, 0.75), nboot = 500)
summary(result)
```

# Acknowledgements

The author acknowledges Anindya Banerjee and Josep Lluís Carrion-i-Silvestre for making their original GAUSS code available, and Joakim Westerlund for the original GAUSS implementation of the Durbin-Hausman test.

# References
