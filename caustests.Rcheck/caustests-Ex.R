pkgname <- "caustests"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "caustests-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('caustests')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("caustests")
### * caustests

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: caustests
### Title: Multiple Granger Causality Tests
### Aliases: caustests

### ** Examples

# Load example data
data(caustests_data)

## No test: 
# Test 1: Toda-Yamamoto test
result1 <- caustests(caustests_data, test = 1, nboot = 199)
print(result1)
summary(result1)

# Test 3: Single Fourier Toda-Yamamoto
result3 <- caustests(caustests_data, test = 3, kmax = 2, nboot = 199)
print(result3)

# Test 6: Quantile causality (fewer quantiles for speed)
result6 <- caustests(caustests_data, test = 6, 
                     quantiles = c(0.25, 0.50, 0.75), nboot = 199)
print(result6)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("caustests", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("caustests_data")
### * caustests_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: caustests_data
### Title: Example Time Series Dataset for Causality Tests
### Aliases: caustests_data
### Keywords: datasets

### ** Examples

data(caustests_data)
head(caustests_data)
summary(caustests_data)

# Check correlations
cor(caustests_data)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("caustests_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.caustests")
### * plot.caustests

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.caustests
### Title: Plot Quantile Causality Results
### Aliases: plot.caustests

### ** Examples

## No test: 
data(caustests_data)
result <- caustests(caustests_data, test = 6,
                    quantiles = c(0.25, 0.50, 0.75), nboot = 199)
plot(result)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.caustests", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
