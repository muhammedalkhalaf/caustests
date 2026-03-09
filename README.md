# caustests: Multiple Granger Causality Tests

Comprehensive suite of Granger causality tests for time series analysis, including:

1. **Toda-Yamamoto (1995)** - Standard Granger causality robust to integration order
2. **Single Fourier Granger** (Enders & Jones, 2016) - Captures smooth structural breaks
3. **Single Fourier Toda-Yamamoto** (Nazlioglu et al., 2016) - Combines TY with Fourier
4. **Cumulative Fourier Granger** (Enders & Jones, 2019) - Multiple Fourier frequencies
5. **Cumulative Fourier Toda-Yamamoto** (Nazlioglu et al., 2019)
6. **Quantile Toda-Yamamoto** (Cai et al., 2023) - Causality across quantiles
7. **Bootstrap Fourier Granger in Quantiles** (Cheng et al., 2021) - BFGC-Q

All tests include bootstrap inference for robust p-values.

## Installation

```r
# Install from CRAN (when available)
install.packages("caustests")

# Or install development version from GitHub
# install.packages("devtools")
devtools::install_github("muhammedalkhalaf/caustests")
```

## Usage

```r
library(caustests)

# Load example data
data(caustests_data)

# Test 1: Toda-Yamamoto test
result1 <- caustests(caustests_data, test = 1, nboot = 999)
print(result1)

# Test 3: Single Fourier Toda-Yamamoto
result3 <- caustests(caustests_data, test = 3, kmax = 3, nboot = 999)
summary(result3)

# Test 6: Quantile causality
result6 <- caustests(caustests_data, test = 6, 
                     quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9),
                     nboot = 999)
print(result6)
plot(result6)
```

## References

- Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. *Journal of Econometrics*, 66(1-2), 225-250.
- Enders, W., & Jones, P. (2016). Grain prices, oil prices, and multiple smooth breaks in a VAR. *Studies in Nonlinear Dynamics & Econometrics*, 20(4), 399-419.
- Nazlioglu, S., Gormus, N. A., & Soytas, U. (2016). Oil prices and real estate investment trusts (REITs). *Energy Economics*, 60, 168-175.
- Nazlioglu, S., Soytas, U., & Gormus, N. A. (2019). Oil prices and monetary policy in emerging markets. *Emerging Markets Finance and Trade*, 55(1), 105-117.
- Cai, Y., Chang, T., Xiang, Y., & Chang, H. L. (2023). Testing Granger causality in quantiles. *Finance Research Letters*, 58, 104327.
- Cheng, S. C., et al. (2021). Bootstrap Fourier Granger causality test in quantiles. *Letters in Spatial and Resource Sciences*, 14, 31-49.

## Author

Dr. Merwan Roudane (merwanroudane920@gmail.com)

## License

GPL-3
