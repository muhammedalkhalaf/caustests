#' @importFrom stats as.formula coef complete.cases fitted pchisq qchisq
#' @importFrom graphics abline barplot par plot
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "caustests ", utils::packageVersion("caustests"), 
    " - Multiple Granger Causality Tests\n",
    "Type 'citation(\"caustests\")' for citation information."
  )
}
