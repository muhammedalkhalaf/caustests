#' Multiple Granger Causality Tests
#'
#' Performs various Granger causality tests including Toda-Yamamoto,
#' Fourier-based tests (single and cumulative frequency), and quantile
#' causality tests with bootstrap inference.
#'
#' @param data A data frame or matrix with time series variables (columns).
#' @param test Integer 1-7 specifying the test type:
#'   \itemize{
#'     \item 1: Toda-Yamamoto (1995)
#'     \item 2: Single Fourier Granger (Enders & Jones, 2016)
#'     \item 3: Single Fourier Toda-Yamamoto (Nazlioglu et al., 2016)
#'     \item 4: Cumulative Fourier Granger (Enders & Jones, 2019)
#'     \item 5: Cumulative Fourier Toda-Yamamoto (Nazlioglu et al., 2019)
#'     \item 6: Quantile Toda-Yamamoto (Cai et al., 2023)
#'     \item 7: Bootstrap Fourier Granger Causality in Quantiles (Cheng et al., 2021)
#'   }
#' @param pmax Maximum lag order for model selection (default: 8).
#' @param ic Information criterion: 1 for AIC, 2 for SBC/BIC (default: 1).
#' @param nboot Number of bootstrap replications (default: 1000).
#' @param kmax Maximum Fourier frequency (default: 3, used for tests 2-5, 7).
#' @param dmax Extra lags for Toda-Yamamoto augmentation. If NULL, automatically
#'   set to 0 for tests 2, 4 (differences) and 1 for tests 1, 3, 5, 6, 7 (levels).
#' @param quantiles Numeric vector of quantiles for tests 6-7
#'   (default: seq(0.1, 0.9, 0.1)).
#' @param verbose Logical; print progress messages (default: TRUE).
#'
#' @return An object of class \code{"caustests"} containing:
#'   \item{results}{Data frame with test results for each direction}
#'   \item{test}{Test number used}
#'   \item{test_name}{Name of the test}
#'   \item{pmax}{Maximum lag considered}
#'   \item{ic}{Information criterion used}
#'   \item{nboot}{Number of bootstrap replications}
#'   \item{kmax}{Maximum Fourier frequency}
#'   \item{dmax}{Augmentation lags}
#'   \item{quantiles}{Quantiles used (for tests 6-7)}
#'   \item{quantile_results}{Detailed quantile results (for tests 6-7)}
#'
#' @details
#' The package implements seven Granger causality tests:
#'
#' \strong{Test 1: Toda-Yamamoto (1995)}
#' Standard Granger causality in levels using VAR with extra lags equal to 
#' the maximum integration order (dmax). This approach is robust to unknown 
#' integration and cointegration properties.
#'
#' \strong{Tests 2-3: Single Fourier Frequency}
#' Incorporate a single Fourier frequency to capture smooth structural breaks.
#' Test 2 uses first differences, Test 3 uses levels (Toda-Yamamoto style).
#'
#' \strong{Tests 4-5: Cumulative Fourier Frequency}
#' Use cumulative Fourier frequencies (1 to k) for more flexible break patterns.
#' Test 4 uses first differences, Test 5 uses levels.
#'
#' \strong{Test 6: Quantile Toda-Yamamoto}
#' Extends Toda-Yamamoto to quantile regression, allowing causality analysis
#' across different quantiles of the conditional distribution.
#'
#' \strong{Test 7: Bootstrap Fourier Granger Causality in Quantiles (BFGC-Q)}
#' Combines Fourier flexibility with quantile regression for robust inference
#' under structural breaks and across quantiles.
#'
#' @references
#' Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector 
#' autoregressions with possibly integrated processes. \emph{Journal of 
#' Econometrics}, 66(1-2), 225-250. \doi{10.1016/0304-4076(94)01616-8}
#'
#' Enders, W., & Jones, P. (2016). Grain prices, oil prices, and multiple 
#' smooth breaks in a VAR. \emph{Studies in Nonlinear Dynamics & Econometrics},
#' 20(4), 399-419. \doi{10.1515/snde-2015-0004}
#'
#' Nazlioglu, S., Gormus, N. A., & Soytas, U. (2016). Oil prices and real 
#' estate investment trusts (REITs): Gradual-shift causality and volatility 
#' transmission analysis. \emph{Energy Economics}, 60, 168-175. 
#' \doi{10.1016/j.eneco.2016.09.009}
#'
#' Nazlioglu, S., Soytas, U., & Gormus, N. A. (2019). Oil prices and monetary 
#' policy in emerging markets: Structural shifts in causal linkages. 
#' \emph{Emerging Markets Finance and Trade}, 55(1), 105-117. 
#' \doi{10.1080/1540496X.2017.1402843}
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
#' \doi{10.1007/s12076-020-00256-3}
#'
#' @examples
#' # Load example data
#' data(caustests_data)
#'
#' \donttest{
#' # Test 1: Toda-Yamamoto test
#' result1 <- caustests(caustests_data, test = 1, nboot = 199)
#' print(result1)
#' summary(result1)
#'
#' # Test 3: Single Fourier Toda-Yamamoto
#' result3 <- caustests(caustests_data, test = 3, kmax = 2, nboot = 199)
#' print(result3)
#'
#' # Test 6: Quantile causality (fewer quantiles for speed)
#' result6 <- caustests(caustests_data, test = 6, 
#'                      quantiles = c(0.25, 0.50, 0.75), nboot = 199)
#' print(result6)
#' }
#'
#' @export
caustests <- function(data, test, pmax = 8, ic = 1, nboot = 1000,
                      kmax = 3, dmax = NULL, 
                      quantiles = seq(0.1, 0.9, 0.1),
                      verbose = TRUE) {
  
  # Input validation
  if (!test %in% 1:7) {
    stop("test must be an integer between 1 and 7")
  }
  if (!ic %in% 1:2) {
    stop("ic must be 1 (AIC) or 2 (SBC/BIC)")
  }
  if (nboot < 99) {
    stop("nboot must be at least 99")
  }
  
  # Convert to matrix
  if (is.data.frame(data)) {
    vnames <- names(data)
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    vnames <- colnames(data)
    if (is.null(vnames)) {
      vnames <- paste0("V", seq_len(ncol(data)))
    }
  } else {
    stop("data must be a data.frame or matrix")
  }
  
  if (ncol(data) < 2) {
    stop("At least 2 variables are required")
  }
  
  # Remove rows with missing values
  complete_rows <- complete.cases(data)
  data <- data[complete_rows, , drop = FALSE]
  n_obs <- nrow(data)
  
  if (n_obs < 2 * pmax + 4) {
    stop(sprintf("Insufficient observations (T=%d) for pmax=%d", n_obs, pmax))
  }
  

  # Auto dmax: 0 for Granger in differences, 1 for TY in levels
  if (is.null(dmax)) {
    dmax <- if (test %in% c(2, 4)) 0 else 1
  }
  
  # kmax only relevant for tests 2-5, 7
  if (!test %in% c(2, 3, 4, 5, 7)) {
    kmax <- 0
  }
  
  # Test names
  test_names <- c(
    "Toda-Yamamoto (1995)",
    "Single Fourier Granger (Enders & Jones, 2016)",
    "Single Fourier Toda-Yamamoto (Nazlioglu et al., 2016)",
    "Cumulative Fourier Granger (Enders & Jones, 2019)",
    "Cumulative Fourier Toda-Yamamoto (Nazlioglu et al., 2019)",
    "Quantile Toda-Yamamoto (Cai et al., 2023)",
    "Bootstrap Fourier Granger Causality in Quantiles (Cheng et al., 2021)"
  )
  
  if (verbose) {
    message(sprintf("\nRunning: %s", test_names[test]))
    message(sprintf("Observations: %d | Max lag: %d | Bootstrap: %d", 
                    n_obs, pmax, nboot))
  }
  
  # Run appropriate test
  if (test %in% 1:5) {
    results <- .run_ols_tests(data, vnames, pmax, ic, test, nboot, kmax, dmax, verbose)
    quantile_results <- NULL
  } else {
    results <- .run_quantile_tests(data, vnames, pmax, ic, test, nboot, kmax, 
                                    dmax, quantiles, verbose)
    quantile_results <- attr(results, "quantile_details")
    attr(results, "quantile_details") <- NULL
  }
  
  # Create output object
  out <- list(
    results = results,
    test = test,
    test_name = test_names[test],
    pmax = pmax,
    ic = ic,
    nboot = nboot,
    kmax = kmax,
    dmax = dmax,
    quantiles = if (test %in% 6:7) quantiles else NULL,
    quantile_results = quantile_results,
    n_obs = n_obs,
    variables = vnames
  )
  
  class(out) <- "caustests"
  return(out)
}


#' @export
print.caustests <- function(x, ...) {
  cat("\n")
  cat(rep("=", 70), "\n", sep = "")
  cat(x$test_name, "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  cat(sprintf("Observations: %d | Variables: %s\n", 
              x$n_obs, paste(x$variables, collapse = ", ")))
  cat(sprintf("Max lag: %d | IC: %s | Bootstrap: %d\n",
              x$pmax, if (x$ic == 1) "AIC" else "SBC", x$nboot))
  
  if (x$test %in% c(2, 3, 4, 5, 7)) {
    cat(sprintf("Max Fourier frequency: %d\n", x$kmax))
  }
  
  cat("\n")
  
  if (x$test %in% 1:5) {
    .print_ols_results(x$results, x$test)
  } else {
    .print_quantile_results(x$results, x$quantile_results, x$quantiles, x$test)
  }
  
  invisible(x)
}


#' @export
summary.caustests <- function(object, ...) {
  print(object, ...)
  
  if (object$test %in% 6:7 && !is.null(object$quantile_results)) {
    cat("\nDetailed Quantile Results:\n")
    cat(rep("-", 50), "\n", sep = "")
    
    for (dir_name in names(object$quantile_results)) {
      cat("\n", dir_name, ":\n", sep = "")
      qres <- object$quantile_results[[dir_name]]
      print(qres, row.names = FALSE)
    }
  }
  
  invisible(object)
}


# Internal: Run OLS-based tests (1-5)
.run_ols_tests <- function(data, vnames, pmax, ic, test, nboot, kmax, dmax, verbose) {
  nvar <- ncol(data)
  n_obs <- nrow(data)
  results <- list()
  
  for (dep_idx in 1:nvar) {
    depvar <- data[, dep_idx]
    dep_name <- vnames[dep_idx]
    
    for (cause_idx in 1:nvar) {
      if (cause_idx == dep_idx) next
      
      cause_name <- vnames[cause_idx]
      direction <- paste0(cause_name, " => ", dep_name)
      
      if (verbose) {
        message(sprintf("  Testing: %s", direction))
      }
      
      # Combine dependent, cause, and controls
      ctrl_idx <- setdiff(1:nvar, c(dep_idx, cause_idx))
      y_combined <- cbind(depvar, data[, cause_idx, drop = FALSE])
      if (length(ctrl_idx) > 0) {
        y_combined <- cbind(y_combined, data[, ctrl_idx, drop = FALSE])
      }
      
      # Select optimal lag (and frequency for Fourier tests)
      sel <- .select_lag_freq(y_combined, pmax, kmax, ic, test, dmax)
      p_opt <- sel$p_opt
      k_opt <- sel$k_opt
      
      # Compute Wald statistic
      wald_result <- .compute_wald_ols(y_combined, p_opt, k_opt, test, dmax, n_obs)
      W <- wald_result$wald
      
      # Asymptotic p-value
      pval_asymp <- if (!is.na(W) && W > 0) {
        stats::pchisq(W, df = p_opt, lower.tail = FALSE)
      } else {
        NA
      }
      
      # Bootstrap p-value
      pval_boot <- .bootstrap_wald_ols(y_combined, p_opt, k_opt, test, dmax, 
                                        n_obs, nboot, W)
      
      results[[direction]] <- data.frame(
        direction = direction,
        wald = W,
        pval_asymp = pval_asymp,
        pval_boot = pval_boot,
        lag = p_opt,
        freq = k_opt,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, results)
}


# Internal: Run quantile tests (6-7)
.run_quantile_tests <- function(data, vnames, pmax, ic, test, nboot, kmax, 
                                 dmax, quantiles, verbose) {
  nvar <- ncol(data)
  results <- list()
  quantile_details <- list()
  
  for (dep_idx in 1:nvar) {
    depvar <- data[, dep_idx]
    dep_name <- vnames[dep_idx]
    
    for (cause_idx in 1:nvar) {
      if (cause_idx == dep_idx) next
      
      cause_name <- vnames[cause_idx]
      direction <- paste0(cause_name, " => ", dep_name)
      
      if (verbose) {
        message(sprintf("  Testing: %s", direction))
      }
      
      # Combine variables
      ctrl_idx <- setdiff(1:nvar, c(dep_idx, cause_idx))
      y_combined <- cbind(depvar, data[, cause_idx, drop = FALSE])
      if (length(ctrl_idx) > 0) {
        y_combined <- cbind(y_combined, data[, ctrl_idx, drop = FALSE])
      }
      
      # Select optimal lag (and frequency)
      sel <- .select_lag_freq_quantile(y_combined, pmax, kmax, ic, test, dmax)
      p_opt <- sel$p_opt
      k_opt <- sel$k_opt
      
      # Test at each quantile
      q_results <- data.frame(
        quantile = quantiles,
        wald = NA_real_,
        pval_boot = NA_real_,
        sig = "",
        stringsAsFactors = FALSE
      )
      
      for (q_idx in seq_along(quantiles)) {
        tau <- quantiles[q_idx]
        qtest <- .quantile_wald_boot(y_combined, p_opt, k_opt, dmax, tau, 
                                      nboot, test)
        q_results$wald[q_idx] <- qtest$wald
        q_results$pval_boot[q_idx] <- qtest$pval_boot
        q_results$sig[q_idx] <- .sig_stars(qtest$pval_boot)
      }
      
      quantile_details[[direction]] <- q_results
      
      # Summary: count significant quantiles
      n_sig_01 <- sum(q_results$pval_boot < 0.01, na.rm = TRUE)
      n_sig_05 <- sum(q_results$pval_boot < 0.05, na.rm = TRUE)
      n_sig_10 <- sum(q_results$pval_boot < 0.10, na.rm = TRUE)
      
      results[[direction]] <- data.frame(
        direction = direction,
        lag = p_opt,
        freq = k_opt,
        n_quantiles = length(quantiles),
        sig_01 = n_sig_01,
        sig_05 = n_sig_05,
        sig_10 = n_sig_10,
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- do.call(rbind, results)
  attr(out, "quantile_details") <- quantile_details
  out
}


# Internal: Select optimal lag and frequency
.select_lag_freq <- function(y, pmax, kmax, ic_type, test, dmax) {
  n_obs <- nrow(y)
  best_ic <- Inf
  p_opt <- 1
  k_opt <- 0
  
  if (kmax == 0) {
    # No Fourier: select p only
    for (p in 1:pmax) {
      ic_val <- .compute_var_ic(y, p + dmax, 0, test, n_obs, ic_type)
      if (!is.na(ic_val) && ic_val < best_ic) {
        best_ic <- ic_val
        p_opt <- p
      }
    }
  } else {
    # Fourier: select (p, k) jointly
    for (k in 1:kmax) {
      for (p in 1:pmax) {
        ic_val <- .compute_var_ic(y, p + dmax, k, test, n_obs, ic_type)
        if (!is.na(ic_val) && ic_val < best_ic) {
          best_ic <- ic_val
          p_opt <- p
          k_opt <- k
        }
      }
    }
  }
  
  list(p_opt = p_opt, k_opt = k_opt)
}


# Internal: Select optimal lag and frequency for quantile tests
.select_lag_freq_quantile <- function(y, pmax, kmax, ic_type, test, dmax) {
  # Use OLS-based selection as initialization
  .select_lag_freq(y, pmax, kmax, ic_type, test, dmax)
}


# Internal: Compute VAR information criterion
.compute_var_ic <- function(y, plags, k_freq, test, T_orig, ic_type) {
  lagged <- .make_lags(y, plags)
  if (is.null(lagged)) return(NA)
  
  dep <- lagged$y
  xl <- lagged$lags
  T_eff <- nrow(dep)
  
  # Add deterministic terms
  det_terms <- .make_deterministic(T_eff, test, k_freq, plags, T_orig)
  z <- cbind(xl, det_terms)
  
  # Check for rank deficiency
  if (ncol(z) >= T_eff) return(NA)
  
  # OLS estimation
  tryCatch({
    qr_z <- qr(z)
    if (qr_z$rank < ncol(z)) return(NA)
    
    b <- qr.solve(qr_z, dep)
    e <- dep - z %*% b
    n <- nrow(e)
    nk <- ncol(dep)
    
    VCV <- crossprod(e) / n
    ldv <- log(det(VCV))
    
    if (!is.finite(ldv)) return(NA)
    
    num_params <- nk * nk * plags + nk
    aic <- ldv + (2 / n) * num_params + nk * (1 + log(2 * pi))
    sbc <- ldv + (log(n) / n) * num_params + nk * (1 + log(2 * pi))
    
    if (ic_type == 1) aic else sbc
  }, error = function(e) NA)
}


# Internal: Create lagged variables
.make_lags <- function(y, plags) {
  if (!is.matrix(y)) y <- as.matrix(y)
  n <- nrow(y)
  k <- ncol(y)
  
  if (n <= plags) return(NULL)
  
  # Create lag matrix
  lag_mat <- matrix(NA_real_, nrow = n, ncol = k * plags)
  
  for (i in 1:k) {
    for (j in 1:plags) {
      col_idx <- j + plags * (i - 1)
      lag_mat[(j + 1):n, col_idx] <- y[1:(n - j), i]
    }
  }
  
  # Trim to valid observations
  valid_rows <- (plags + 1):n
  
  list(
    y = y[valid_rows, , drop = FALSE],
    lags = lag_mat[valid_rows, , drop = FALSE]
  )
}


# Internal: Create deterministic terms
.make_deterministic <- function(T_eff, test, k_freq, plags, T_orig) {
  if (test == 1 || k_freq == 0) {
    # Constant only
    return(matrix(1, nrow = T_eff, ncol = 1))
  }
  
  # Time index for Fourier terms
  t_seq <- seq_len(T_eff)
  
  if (test %in% c(2, 3)) {
    # Single frequency
    fourier <- cbind(
      sin(2 * pi * k_freq * t_seq / T_eff),
      cos(2 * pi * k_freq * t_seq / T_eff)
    )
  } else {
    # Cumulative frequencies (tests 4, 5, 7)
    fourier <- NULL
    for (ki in 1:k_freq) {
      fourier <- cbind(
        fourier,
        sin(2 * pi * ki * t_seq / T_eff),
        cos(2 * pi * ki * t_seq / T_eff)
      )
    }
  }
  
  cbind(1, fourier)
}


# Internal: Compute Wald statistic (OLS)
.compute_wald_ols <- function(y, p_opt, k_opt, test, dmax, T_orig) {
  p_full <- p_opt + dmax
  
  lagged <- .make_lags(y, p_full)
  if (is.null(lagged)) return(list(wald = NA))
  
  dep <- lagged$y[, 1]
  xl <- lagged$lags
  T_eff <- length(dep)
  
  det_terms <- .make_deterministic(T_eff, test, k_opt, p_full, T_orig)
  z <- cbind(xl, det_terms)
  
  nc <- ncol(z)
  df <- T_eff - nc
  
  if (df <= 0) return(list(wald = NA))
  
  tryCatch({
    # Unrestricted model
    qr_z <- qr(z)
    b <- qr.solve(qr_z, dep)
    u <- dep - z %*% b
    RSS_ur <- sum(u^2)
    
    # Restricted model: exclude cause lags 1:p_opt
    # Cause variable lags are columns (p_full + 1) to (p_full + p_opt)
    # But we need to exclude columns for the cause variable (2nd in y)
    # In the lag matrix: first p_full columns are dep lags, next p_full are cause lags
    
    cause_cols <- (p_full + 1):(p_full + p_opt)
    
    # But actually the indexing depends on ncol(y)
    nvar <- ncol(y)
    # Column structure: for each variable i, lags 1:p_full are at positions
    # (i-1)*p_full + 1 to (i-1)*p_full + p_full
    # For cause (variable 2), testing lags 1:p_opt
    cause_test_cols <- p_full + 1:p_opt
    
    zr <- z[, -cause_test_cols, drop = FALSE]
    
    qr_zr <- qr(zr)
    br <- qr.solve(qr_zr, dep)
    ur <- dep - zr %*% br
    RSS_r <- sum(ur^2)
    
    # F-test to Wald
    F_stat <- ((RSS_r - RSS_ur) / p_opt) / (RSS_ur / df)
    W <- F_stat * p_opt
    
    list(wald = max(W, 0))
  }, error = function(e) list(wald = NA))
}


# Internal: Bootstrap Wald (OLS)
.bootstrap_wald_ols <- function(y, p_opt, k_opt, test, dmax, T_orig, nboot, W_obs) {
  if (is.na(W_obs)) return(NA)
  
  p_full <- p_opt + dmax
  
  lagged <- .make_lags(y, p_full)
  if (is.null(lagged)) return(NA)
  
  dep <- lagged$y[, 1]
  xl <- lagged$lags
  T_eff <- length(dep)
  
  det_terms <- .make_deterministic(T_eff, test, k_opt, p_full, T_orig)
  z <- cbind(xl, det_terms)
  
  # Cause test columns
  cause_test_cols <- p_full + 1:p_opt
  zr <- z[, -cause_test_cols, drop = FALSE]
  
  # Restricted model
  tryCatch({
    qr_zr <- qr(zr)
    br <- qr.solve(qr_zr, dep)
    yhat <- zr %*% br
    resid <- dep - yhat
    resid <- resid - mean(resid)
    
    W_boot <- numeric(nboot)
    
    for (b in seq_len(nboot)) {
      # Resample residuals
      idx <- sample.int(T_eff, T_eff, replace = TRUE)
      e_star <- resid[idx]
      e_star <- e_star - mean(e_star)
      y_star <- yhat + e_star
      
      # Compute bootstrap Wald
      W_boot[b] <- tryCatch({
        # Unrestricted
        qr_z <- qr(z)
        b_ur <- qr.solve(qr_z, y_star)
        u_ur <- y_star - z %*% b_ur
        RSS_ur <- sum(u_ur^2)
        
        # Restricted
        b_r <- qr.solve(qr_zr, y_star)
        u_r <- y_star - zr %*% b_r
        RSS_r <- sum(u_r^2)
        
        df <- T_eff - ncol(z)
        F_stat <- ((RSS_r - RSS_ur) / p_opt) / (RSS_ur / df)
        max(F_stat * p_opt, 0)
      }, error = function(e) 0)
    }
    
    mean(W_boot >= W_obs)
  }, error = function(e) NA)
}


# Internal: Quantile Wald with bootstrap
.quantile_wald_boot <- function(y, p_opt, k_opt, dmax, tau, nboot, test) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for quantile tests")
  }
  
  p_full <- p_opt + dmax
  
  lagged <- .make_lags(y, p_full)
  if (is.null(lagged)) return(list(wald = NA, pval_boot = NA))
  
  dep <- lagged$y[, 1]
  xl <- lagged$lags
  T_eff <- length(dep)
  
  # Deterministic terms (only for test 7)
  if (test == 7 && k_opt > 0) {
    t_seq <- seq_len(T_eff)
    fourier <- NULL
    for (ki in 1:k_opt) {
      fourier <- cbind(
        fourier,
        sin(2 * pi * ki * t_seq / T_eff),
        cos(2 * pi * ki * t_seq / T_eff)
      )
    }
    z <- cbind(1, xl, fourier)
  } else {
    z <- cbind(1, xl)
  }
  
  # Cause test columns (after intercept)
  cause_test_cols <- 1 + p_full + 1:p_opt
  
  # Create data frame for qreg
  df_full <- data.frame(y = dep, z)
  col_names <- paste0("x", seq_len(ncol(z)))
  names(df_full) <- c("y", col_names)
  
  # Formula for full model
  form_full <- stats::as.formula(paste("y ~", paste(col_names, collapse = " + "), "- 1"))
  
  # Restricted model columns
  restr_cols <- col_names[-cause_test_cols + 1]  # Adjust for naming
  form_restr <- stats::as.formula(paste("y ~", paste(restr_cols, collapse = " + "), "- 1"))
  
  tryCatch({
    # Full model at tau
    fit_full <- quantreg::rq(form_full, data = df_full, tau = tau)
    
    # Extract coefficients for cause variables
    cause_coef_names <- col_names[cause_test_cols - 1]
    b_cause <- stats::coef(fit_full)[cause_coef_names]
    V_full <- tryCatch({
      summary(fit_full, se = "nid", covariance = TRUE)$cov
    }, error = function(e) {
      summary(fit_full, se = "iid", covariance = TRUE)$cov
    })
    
    if (is.null(V_full)) {
      return(list(wald = NA, pval_boot = NA))
    }
    
    # Extract submatrix for cause coefficients
    cause_idx <- match(cause_coef_names, names(stats::coef(fit_full)))
    V_cause <- V_full[cause_idx, cause_idx, drop = FALSE]
    
    # Wald statistic
    W_obs <- tryCatch({
      as.numeric(t(b_cause) %*% solve(V_cause) %*% b_cause)
    }, error = function(e) NA)
    
    if (is.na(W_obs) || W_obs < 0) W_obs <- 0
    
    # Bootstrap
    fit_restr <- quantreg::rq(form_restr, data = df_full, tau = tau)
    yhat <- stats::fitted(fit_restr)
    resid <- dep - yhat
    resid <- resid - mean(resid)
    
    W_boot_ge <- 0
    
    for (b in seq_len(nboot)) {
      idx <- sample.int(T_eff, T_eff, replace = TRUE)
      e_star <- resid[idx]
      e_star <- e_star - mean(e_star)
      df_boot <- df_full
      df_boot$y <- yhat + e_star
      
      W_b <- tryCatch({
        fit_b <- quantreg::rq(form_full, data = df_boot, tau = tau)
        b_b <- stats::coef(fit_b)[cause_coef_names]
        V_b <- tryCatch({
          summary(fit_b, se = "nid", covariance = TRUE)$cov
        }, error = function(e) {
          tryCatch(summary(fit_b, se = "iid", covariance = TRUE)$cov, error = function(e2) NULL)
        })
        if (is.null(V_b)) return(0)
        V_b_cause <- V_b[cause_idx, cause_idx, drop = FALSE]
        as.numeric(t(b_b) %*% solve(V_b_cause) %*% b_b)
      }, error = function(e) 0)
      
      if (!is.na(W_b) && W_b >= W_obs) {
        W_boot_ge <- W_boot_ge + 1
      }
    }
    
    list(wald = W_obs, pval_boot = W_boot_ge / nboot)
  }, error = function(e) {
    list(wald = NA, pval_boot = NA)
  })
}


# Internal: Significance stars
.sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  ""
}


# Internal: Print OLS results
.print_ols_results <- function(results, test) {
  cat(rep("-", 76), "\n", sep = "")
  cat(sprintf("%-32s %8s %9s %9s %4s %4s\n",
              "Null Hypothesis", "Wald", "Asym.p", "Boot.p", "Lag", "Freq"))
  cat(rep("-", 76), "\n", sep = "")
  
  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    sig <- .sig_stars(r$pval_boot)
    dir_trunc <- substr(r$direction, 1, 32)
    cat(sprintf("%-32s %8.3f %9.3f %9.3f%s %4d %4d\n",
                dir_trunc, r$wald, r$pval_asymp, r$pval_boot, sig,
                r$lag, r$freq))
  }
  
  cat(rep("-", 76), "\n", sep = "")
  cat("Note: *p<0.10, **p<0.05, ***p<0.01 based on bootstrap p-values\n")
  if (test %in% c(2, 3, 4, 5)) {
    cat("      Freq = optimal Fourier frequency (k*)\n")
  }
  cat("\n")
}


# Internal: Print quantile results
.print_quantile_results <- function(summary_results, quantile_details, quantiles, test) {
  cat(rep("-", 60), "\n", sep = "")
  cat(sprintf("%-32s %4s %4s %4s %4s %4s\n",
              "Direction", "Lag", "Freq", "Sig1%", "Sig5%", "Sig10%"))
  cat(rep("-", 60), "\n", sep = "")
  
  for (i in seq_len(nrow(summary_results))) {
    r <- summary_results[i, ]
    dir_trunc <- substr(r$direction, 1, 32)
    cat(sprintf("%-32s %4d %4d %5d %5d %6d\n",
                dir_trunc, r$lag, r$freq, r$sig_01, r$sig_05, r$sig_10))
  }
  
  cat(rep("-", 60), "\n", sep = "")
  cat(sprintf("Quantiles tested: %s\n", paste(quantiles, collapse = ", ")))
  cat("Sig columns show count of significant quantiles at each level\n")
  cat("\n")
}


#' Plot Quantile Causality Results
#'
#' Creates diagnostic plots for quantile causality tests (tests 6-7).
#'
#' @param x An object of class \code{"caustests"} from test 6 or 7.
#' @param which Which direction to plot (default: 1, first direction).
#' @param type Plot type: "wald" for Wald statistics, "pval" for p-values,
#'   or "both" (default).
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisibly returns the plotted data.
#'
#' @examples
#' \donttest{
#' data(caustests_data)
#' result <- caustests(caustests_data, test = 6,
#'                     quantiles = c(0.25, 0.50, 0.75), nboot = 199)
#' plot(result)
#' }
#'
#' @export
plot.caustests <- function(x, which = 1, type = "both", ...) {
  if (!x$test %in% 6:7) {
    message("Plotting is only available for quantile tests (6-7)")
    return(invisible(NULL))
  }
  
  if (is.null(x$quantile_results)) {
    message("No quantile results available")
    return(invisible(NULL))
  }
  
  dirs <- names(x$quantile_results)
  if (which > length(dirs)) {
    which <- 1
  }
  
  qdata <- x$quantile_results[[dirs[which]]]
  dir_name <- dirs[which]
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (type == "both") {
    graphics::par(mfrow = c(1, 2))
  }
  
  if (type %in% c("wald", "both")) {
    graphics::barplot(
      qdata$wald,
      names.arg = qdata$quantile,
      main = paste("Wald Statistics:", dir_name),
      xlab = expression(tau),
      ylab = "Wald statistic",
      col = "steelblue",
      border = "white",
      ...
    )
    graphics::abline(h = stats::qchisq(0.95, df = x$results$lag[which]),
                     col = "red", lty = 2)
  }
  
  if (type %in% c("pval", "both")) {
    graphics::plot(
      qdata$quantile, qdata$pval_boot,
      type = "b",
      main = paste("Bootstrap P-values:", dir_name),
      xlab = expression(tau),
      ylab = "P-value",
      ylim = c(0, 1),
      pch = 19,
      col = "steelblue",
      ...
    )
    graphics::abline(h = 0.05, col = "red", lty = 2)
    graphics::abline(h = 0.10, col = "orange", lty = 3)
  }
  
  invisible(qdata)
}
