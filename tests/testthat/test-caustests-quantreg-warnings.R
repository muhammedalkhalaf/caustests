test_that("repeated quantreg sparsity warnings are counted and muffled", {
  # Fake computation that warns three times (twice the quantreg sparsity
  # warning, once something else); only the sparsity warnings must be muffled.
  fake_core <- function(y, p_opt, k_opt, dmax, tau, nboot, test) {
    warning("3 non-positive fis")
    warning("7 non-positive fis")
    warning("some other warning")
    list(wald = 1.5, pval_boot = 0.2)
  }
  local_mocked_bindings(.quantile_wald_boot_core = fake_core)
  res <- NULL
  w <- capture_warnings(
    res <- caustests:::.quantile_wald_boot(NULL, 1L, 0L, 1L, 0.5, 99L, 6L)
  )
  expect_identical(w, "some other warning")
  expect_identical(res$n_fis_warnings, 2L)
  expect_identical(res$n_fits, 100L)
  expect_equal(res$wald, 1.5)
  expect_equal(res$pval_boot, 0.2)
})

test_that("caustests() emits exactly one summary warning and records the counts", {
  skip_if_not_installed("quantreg")
  data(caustests_data, envir = environment())
  d <- as.matrix(caustests_data[, 1:2])
  set.seed(1)
  n_fis <- 0L
  other <- character()
  res <- withCallingHandlers(
    caustests(d, test = 6, quantiles = c(0.1, 0.5), nboot = 99,
              pmax = 2, verbose = FALSE),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("non-positive fis", msg, fixed = TRUE)) {
        n_fis <<- n_fis + 1L
      } else {
        other <<- c(other, msg)
      }
      invokeRestart("muffleWarning")
    })
  # One warning in total that mentions the quantreg message, none of the
  # individual quantreg warnings
  expect_identical(n_fis, 1L)
  expect_identical(other, character())
  fl <- res$quantreg_warnings
  expect_s3_class(fl, "data.frame")
  expect_named(fl, c("direction", "quantile", "n_warnings", "n_fits"))
  expect_equal(nrow(fl), 2L * 2L)  # two directions, two quantiles
  expect_true(all(fl$n_fits == 100L))
  expect_true(sum(fl$n_warnings) > 0)
  expect_true(all(fl$n_warnings <= fl$n_fits))
  # print mentions it in one line
  out <- capture.output(print(res))
  expect_length(grep("quantreg reported non-positive fis", out), 1L)
})

test_that("no summary warning when quantreg did not warn", {
  fl <- data.frame(direction = "a => b", quantile = c(0.25, 0.75),
                   n_warnings = c(0L, 0L), n_fits = c(100L, 100L),
                   stringsAsFactors = FALSE)
  expect_no_warning(caustests:::.warn_quantreg_fis(fl))
  fl$n_warnings[2] <- 5L
  expect_warning(caustests:::.warn_quantreg_fis(fl),
                 "1 of 2 quantiles \\(tau = 0.75; 5 of 200 quantile regression fits\\)")
})
