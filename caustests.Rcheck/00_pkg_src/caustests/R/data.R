#' Example Time Series Dataset for Causality Tests
#'
#' A simulated dataset containing three time series variables for demonstrating
#' Granger causality tests. The data includes one dependent variable (Y) and
#' two potential causal variables (X1, X2) with known causal relationships.
#'
#' @format A data frame with 200 observations and 3 variables:
#' \describe{
#'   \item{Y}{Dependent variable, generated as AR(2) plus causal effects from X1}
#'   \item{X1}{First explanatory variable, AR(1) process}
#'   \item{X2}{Second explanatory variable, independent AR(1) process}
#' }
#'
#' @details
#' The data generating process is:
#' \itemize{
#'   \item X1 and X2 are independent AR(1) processes
#'   \item Y depends on its own lags plus lagged values of X1 (but not X2)
#'   \item This creates a true causal relationship from X1 to Y
#'   \item There is no true causality from X2 to Y or from Y to X1/X2
#' }
#'
#' This allows users to verify that the causality tests correctly identify
#' the causal direction X1 => Y while finding no significant causality
#' in other directions (with appropriate sample sizes and test settings).
#'
#' @examples
#' data(caustests_data)
#' head(caustests_data)
#' summary(caustests_data)
#'
#' # Check correlations
#' cor(caustests_data)
#'
#' @source Simulated data for package demonstration
#' @keywords datasets
"caustests_data"
