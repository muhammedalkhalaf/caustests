#' caustests: Multiple Granger Causality Tests for Time Series and Panel Data
#'
#' Comprehensive suite of Granger causality tests for time series and panel
#' data including standard Toda-Yamamoto, Fourier-based tests, quantile
#' causality tests, panel Fourier Toda-Yamamoto, panel quantile causality,
#' and Group-Mean/Pooled FM-OLS for panel cointegrating polynomial regressions.
#'
#' @section Time Series Functions:
#' \itemize{
#'   \item \code{\link{caustests}}: Perform Granger causality tests (time series)
#' }
#'
#' @section Panel Data Functions:
#' \itemize{
#'   \item \code{\link{xtpcaus}}: Panel Granger causality tests (PFTY and PQC)
#'   \item \code{\link{xtpcmg}}: Panel cointegrating polynomial regressions (FM-OLS)
#' }
#'
#' @section Available Tests (Time Series):
#' \enumerate{
#'   \item Toda-Yamamoto (1995) - Robust to integration order
#'   \item Single Fourier Granger (Enders & Jones, 2016)
#'   \item Single Fourier Toda-Yamamoto (Nazlioglu et al., 2016)
#'   \item Cumulative Fourier Granger (Enders & Jones, 2019)
#'   \item Cumulative Fourier Toda-Yamamoto (Nazlioglu et al., 2019)
#'   \item Quantile Toda-Yamamoto (Cai et al., 2023)
#'   \item Bootstrap Fourier Granger in Quantiles (Cheng et al., 2021)
#' }
#'
#' @section Available Tests (Panel Data):
#' \enumerate{
#'   \item Panel Fourier Toda-Yamamoto (Yilanci and Gorus, 2020)
#'   \item Panel Quantile Causality (Wang and Nguyen, 2022)
#'   \item Group-Mean FM-OLS (Wagner and Reichold, 2023)
#'   \item Pooled FM-OLS (de Jong and Wagner, 2022)
#' }
#'
#' @section Data:
#' \itemize{
#'   \item \code{\link{caustests_data}}: Example time series dataset
#'   \item \code{\link{grunfeld_panel}}: Example panel dataset for xtpcaus
#'   \item \code{\link{grunfeld_cmg}}: Example panel dataset for xtpcmg
#' }
#'
#' @references
#' Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector 
#' autoregressions with possibly integrated processes. \emph{Journal of 
#' Econometrics}, 66(1-2), 225-250. \doi{10.1016/0304-4076(94)01616-8}
#'
#' Enders, W., & Jones, P. (2016). Grain prices, oil prices, and multiple 
#' smooth breaks in a VAR. \emph{Studies in Nonlinear Dynamics & Econometrics},
#' 20(4), 399-419. \doi{10.1515/snde-2014-0101}
#'
#' Nazlioglu, S., Gormus, N. A., & Soytas, U. (2016). Oil prices and real 
#' estate investment trusts (REITs): Gradual-shift causality and volatility 
#' transmission analysis. \emph{Energy Economics}, 60, 168-175. 
#' \doi{10.1016/j.eneco.2016.09.009}
#'
#' Nazlioglu, S., Soytas, U., & Gormus, N. A. (2019). Oil prices and monetary 
#' policy in emerging markets: Structural shifts in causal linkages. 
#' \emph{Emerging Markets Finance and Trade}, 55(1), 105-117. 
#' \doi{10.1080/1540496X.2018.1434072}
#'
#' Cai, Y., Chang, T., Xiang, Y., & Chang, H. L. (2023). Testing Granger 
#' causality in quantiles between the stock and the foreign exchange markets 
#' of Japan. \emph{Finance Research Letters}, 58, 104327. 
#' \doi{10.1016/j.frl.2023.104327}
#'
#' Cheng, S. C., Hsueh, H. P., Ranjbar, O., Wang, M. C., & Chang, T. (2021). 
#' Bootstrap Fourier Granger causality test in quantiles and the asymmetric 
#' causal relationship between CO2 emissions and economic growth. 
#' \emph{Letters in Spatial and Resource Sciences}, 14, 31-49. 
#' \doi{10.1007/s12076-020-00263-0}
#'
#' Yilanci, V. and Gorus, M.S. (2020). Does economic globalization have 
#' predictive power for ecological footprint in MENA counties? A panel 
#' causality test with a Fourier function. \emph{Environmental Science and 
#' Pollution Research}, 27, 40552-40562. \doi{10.1007/s11356-020-10092-9}
#'
#' Wang, K.M. and Nguyen, T.B. (2022). A quantile panel-type analysis of 
#' income inequality and healthcare expenditure. \emph{Economic Research}, 
#' 35(1), 873-893. \doi{10.1080/1331677X.2021.1952089}
#'
#' Wagner, M. and Reichold, K. (2023). Panel cointegrating polynomial 
#' regressions: group-mean fully modified OLS estimation and inference. 
#' \emph{Econometric Reviews}, 42(4), 358-392. 
#' \doi{10.1080/07474938.2023.2178141}
#'
#' de Jong, R.M. and Wagner, M. (2022). Panel cointegrating polynomial 
#' regressions. \emph{Annals of Applied Statistics}, 16(1), 416-442. 
#' \doi{10.1214/21-AOAS1536}
#'
#' @docType package
#' @name caustests-package
#' @aliases caustests-package
#' @keywords internal
"_PACKAGE"
