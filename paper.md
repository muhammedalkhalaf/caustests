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

Eight R packages implement cointegration and Granger causality methods absent from the R ecosystem: `cointsmall` for small-sample cointegration, `xtcadfcoint` and `xtbreakcoint` for panel cointegration with breaks and common factors, `hatemicoint` and `makicoint` for regime-shift cointegration, `xtdhcoint` for Durbin–Hausman panel tests, `fcoint` for Fourier-based cointegration, and `caustests` for Fourier/quantile/bootstrap causality. All are open-source under GPL-3.

# Statement of Need

The `urca` package [@Pfaff2008] provides classical Engle–Granger and Johansen tests, yet several important advances remain unavailable in R: small-sample methods [@Trinh2022], panel tests with breaks and common factors [@Banerjee2015; @Banerjee2024], regime-shift cointegration [@HatemiJ2008; @Maki2012], Fourier-based tests [@Banerjee2017fadl; @Tsong2016], and advanced causality tests [@Enders2016; @Cai2023; @Cheng2021]. These packages fill those gaps.

# Packages

## cointsmall

Implements the @Trinh2022 cointegration test for very small samples (as few as 12 observations) with up to two endogenous structural breaks.

```r
library(cointsmall)
result <- cointsmall(y, x, breaks = 2, model = "cs",
                     criterion = "adf", trim = 0.15)
summary(result)
```

## xtcadfcoint

Implements the @Banerjee2024 panel cointegration test with structural instabilities using cross-sectionally augmented Dickey–Fuller regressions and Common Correlated Effects estimation.

```r
library(xtcadfcoint)
result <- xtcadfcoint(y ~ x1 + x2, data = panel_data,
                      id = "country", time = "year",
                      model = 1, breaks = 1)
summary(result)
```

## xtbreakcoint

Panel cointegration with structural breaks and cross-section dependence following @Banerjee2015, with iterative factor-break estimation and the @BaiNg2004 MQ test.

```r
library(xtbreakcoint)
result <- xtbreakcoint(y ~ x1 + x2, data = panel_data,
                       id = "country", time = "year",
                       model = "trendshift", max_factors = 5)
summary(result)
```

## hatemicoint

Implements the @HatemiJ2008 cointegration test with two unknown regime shifts, providing ADF*, Zt*, and Zα* statistics with endogenous break dates.

```r
library(hatemicoint)
result <- hatemicoint(y, x, maxlags = 8,
                      lag_selection = "aic", trimming = 0.15)
summary(result)
```

## xtdhcoint

Implements the @Westerlund2008 Durbin–Hausman panel cointegration tests, robust to cross-sectional dependence via principal components. Provides group-mean (DHg) and panel (DHp) statistics.

```r
library(xtdhcoint)
result <- xtdhcoint(y ~ x1 + x2, data = panel_data,
                    id = "country", time = "year",
                    criterion = "bic")
summary(result)
```

## fcoint

Four Fourier-based cointegration tests: FADL [@Banerjee2017fadl], FEG, FEG2, and the @Tsong2016 KPSS-type test.

```r
library(fcoint)
result <- fcoint(y, x, test = "fadl", max_freq = 3,
                 criterion = "aic")
print(result)
```

## makicoint

Implements the @Maki2012 test allowing an unknown number of breaks (up to five) under four model specifications: level shift, level shift with trend, regime shift, and regime shift with trend.

```r
library(makicoint)
result <- makicoint(y, x, model = 2, max_breaks = 3,
                    trimming = 0.10)
summary(result)
```

## caustests

Granger causality suite: Toda–Yamamoto [@TodaYamamoto1995], Fourier single-frequency [@Enders2016] and cumulative-frequency, quantile causality [@Cai2023], and bootstrap Fourier Granger causality in quantiles [@Cheng2021]. The `test` argument (1–7) selects the method.

```r
library(caustests)
# Fourier Toda-Yamamoto (test = 2)
result <- caustests(data = cbind(y, x), test = 2,
                    pmax = 4, kmax = 3, nboot = 1000)
summary(result)

# Quantile causality (test = 6)
result <- caustests(data = cbind(y, x), test = 6,
                    quantiles = c(0.25, 0.5, 0.75), nboot = 500)
summary(result)
```

# Acknowledgements

The author acknowledges Anindya Banerjee and Josep Lluís Carrion-i-Silvestre for their original GAUSS code, and Joakim Westerlund for the Durbin–Hausman GAUSS implementation.

# References
