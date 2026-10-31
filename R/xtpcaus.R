#' Panel Granger Causality Tests
#'
#' Tests whether \code{x} Granger-causes \code{y} in a balanced panel using
#' either the Panel Fourier Toda-Yamamoto (PFTY) test or the Panel Quantile
#' Causality (PQC) test.
#'
#' @details
#' \strong{PFTY specification.} For each unit \code{i} the test regression is
#' \deqn{y_t = a + \sum_{k=1}^{p_i + dmax} g_k y_{t-k} + \sum_{k=1}^{p_i + dmax}
#' b_k x_{t-k} + c_1 \sin(2 \pi f_i t / T) + c_2 \cos(2 \pi f_i t / T) + e_t,}
#' with \code{p_i} in \code{1..pmax} and \code{f_i} in \code{1..kmax} chosen
#' jointly by the information criterion. The null hypothesis restricts the
#' first \code{p_i} lags of \code{x}. The individual Wald statistic uses the
#' residual variance \code{RSS / n} with \code{n = T - p_i - dmax} the
#' effective sample (the normalisation of the Stata command caustests).
#' \code{kmax = 0} omits the Fourier terms; \code{dmax = 0} omits the
#' Toda-Yamamoto augmentation. With \code{kmax = 0} and \code{dmax = 0} the
#' regression is the plain VAR(p) of Dumitrescu and Hurlin (2012, eq. 1) with
#' a unit-specific lag order \code{K_i = p_i}.
#'
#' \strong{Degrees of freedom.} A lag order \code{p} is admissible only when
#' the unrestricted regression keeps at least 5 residual degrees of freedom,
#' \code{T - (p + dmax) - k >= 5} with \code{k = 2 (p + dmax) + 1 + 2 I(kmax > 0)}
#' regressors. Lag orders that violate this bound are not considered (a
#' warning says so); if no lag order is admissible the function stops. This
#' bound implies the condition \code{T_i > 5 + 2 K_i} of Dumitrescu and Hurlin
#' (2012, eq. 10) for the existence of the second moments of the Wald
#' statistic.
#'
#' \strong{Dumitrescu and Hurlin (2012) statistics.} Only for the plain
#' specification (\code{kmax = 0} and \code{dmax = 0}) the function reports the
#' standardised statistics of Dumitrescu and Hurlin (2012). They are built
#' from the individual Wald statistics with the unbiased residual variance
#' \code{RSS / (T_i - 2 K_i - 1)} (the paper's main definition), where
#' \code{T_i = T - p_i} is the effective sample of unit \code{i}:
#' \code{Z-bar} is eq. (9) with the asymptotic moments \code{K_i} and
#' \code{2 K_i}, and \code{Z-tilde} is eq. (26) with the fixed-T moments of
#' eq. (21) and (22),
#' \deqn{E(W_i) = K_i (T_i - 2 K_i - 1) / (T_i - 2 K_i - 3),}
#' \deqn{Var(W_i) = 2 K_i (T_i - 2 K_i - 1)^2 (T_i - K_i - 3) /
#' ((T_i - 2 K_i - 3)^2 (T_i - 2 K_i - 5)).}
#' When the selected lag orders differ across units the moments are averaged
#' over units: Z-bar uses the averaged asymptotic moments of eq. (11) of the
#' paper (mean of K_i and of 2K_i), Z-tilde the averaged fixed-T moments of
#' eq. (33). Both p-values are one-sided (upper
#' tail), the decision rule of the paper. \code{Z-tilde} is the statistic to
#' use at fixed T; \code{Z-bar} requires T to grow before N and is oversized
#' in short panels (the paper's Table 1). For the Fourier or Toda-Yamamoto
#' augmented specification the Wald statistics come from a different
#' regression with a data-selected frequency, the moments of eq. (21) and (22)
#' do not apply, and \code{zbar}, \code{ztilde} and their p-values are
#' \code{NA}; inference then rests on the bootstrap Fisher p-value.
#'
#' \strong{Bootstrap.} The bootstrap is a residual bootstrap under the null
#' built on Dumitrescu and Hurlin (2012, Section 6.2, steps 3 to 7). For
#' each unit the model is re-estimated without the tested lags of \code{x}
#' (restricted model) and its residuals are kept. In each replication a
#' sequence of time indices is drawn with replacement from the periods common
#' to all units and the same sequence is used for every unit (block of size 1
#' in time and size N in the cross-section), so that cross-sectional
#' dependence of the residuals is preserved. A bootstrap series is generated
#' recursively from the restricted model (intercept, autoregressive part,
#' unrestricted lags of \code{x}, Fourier terms, resampled residuals; the
#' observed first \code{p_i + dmax} values are the initial conditions), and
#' the whole procedure, selection of the lag order and frequency included, is
#' repeated on it. The paper's step 5 holds the lagged values of \code{y}
#' fixed; the recursive construction is used here so that the selection step
#' can be repeated, since the selection on the observed data inflates the
#' Wald statistics in short panels. When there is a single candidate
#' (\code{pmax = 1} and \code{kmax <= 1}) the selection is skipped. The
#' individual bootstrap p-values are combined by the Fisher statistic of
#' Emirmahmutoglu and Kose (2011), \code{-2 sum log p_i} on \code{2 N}
#' degrees of freedom. For the plain specification the bootstrap
#' distribution of \code{Z-tilde} is also recorded and \code{ztilde_boot_pv}
#' is the one-sided bootstrap p-value. Bootstrap regressions that are
#' singular or non-finite are dropped with a warning giving their number.
#' In Monte Carlo experiments under the null with N = 10, pmax = 1, dmax = 1,
#' kmax = 3 and nominal size 5 percent, the bootstrap Fisher test rejected
#' about 12 percent at T = 20 and 9 percent at T = 50 (independent normal
#' errors), while the plain specification was at nominal size; the
#' selection of the Fourier frequency on short panels is absorbed only in
#' part. For the Dumitrescu and Hurlin statistics use the plain
#' specification.
#'
#' \strong{Singular regressions.} If the regressor matrix of a unit is
#' singular for every admissible lag order (for example \code{x} constant or
#' identical to \code{y}), its Wald statistic is \code{NA}, a warning names
#' the unit and the unit is excluded from the panel statistics; \code{N_used}
#' gives the number of units that enter them.
#'
#' \strong{Input checks.} Missing values in \code{y}, \code{x},
#' \code{panel_id} or \code{time_id} stop the function with the count and the
#' row positions. Every unit must be observed at the same set of time values
#' exactly once; duplicated or misaligned rows stop the function.
#'
#' \strong{Random numbers.} When \code{seed >= 0} the bootstrap uses that
#' seed in a local RNG scope: the state of the global random number generator
#' is saved before and restored after the call.
#'
#' @param data A data frame in long format.
#' @param y Character. Name of the dependent (caused) variable.
#' @param x Character. Name of the independent (causing) variable.
#' @param panel_id Character. Name of the panel identifier variable.
#' @param time_id Character. Name of the time variable.
#' @param test Character. Test type: \code{"pfty"} for Panel Fourier
#'   Toda-Yamamoto or \code{"pqc"} for Panel Quantile Causality.
#' @param pmax Integer. Maximum lag order for selection. Default is \code{4}.
#'   Lag orders that leave fewer than 5 residual degrees of freedom are not
#'   considered (see Details).
#' @param dmax Integer. Maximum integration order for Toda-Yamamoto
#'   augmentation. Default is \code{1}. \code{0} omits the augmentation.
#' @param nboot Integer. Number of bootstrap replications. Minimum 99.
#'   Default is \code{499}.
#' @param kmax Integer. Maximum Fourier frequency (PFTY only). Default is
#'   \code{3}. \code{0} omits the Fourier terms.
#' @param ic Character. Information criterion: \code{"aic"} or \code{"bic"}.
#'   Default is \code{"aic"}.
#' @param quantiles Numeric vector. Quantile grid for PQC test (values strictly
#'   between 0 and 1). Default is
#'   \code{c(0.1, 0.25, 0.50, 0.75, 0.90)}.
#' @param seed Integer. Random seed for bootstrap. \code{-1} means no seed.
#'   Default is \code{-1}. The global RNG state is restored on exit.
#'
#' @return An object of class \code{"xtpcaus"} containing:
#' \describe{
#'   \item{test}{Character. \code{"pfty"} or \code{"pqc"}.}
#'   \item{N}{Integer. Number of panel units.}
#'   \item{TT}{Integer. Number of time periods.}
#'   \item{nboot}{Integer. Number of bootstrap replications.}
#'   \item{y}{Character. Name of the y variable.}
#'   \item{x}{Character. Name of the x variable.}
#'   For PFTY:
#'   \item{fisher}{Numeric. Fisher panel statistic, \code{-2 sum log p_i}
#'     over the bootstrap p-values of the units that enter the panel
#'     statistics.}
#'   \item{fisher_df}{Integer. Degrees of freedom (2 N_used).}
#'   \item{fisher_pv}{Numeric. Fisher p-value.}
#'   \item{N_used}{Integer. Number of units entering the panel statistics
#'     (units with a singular regression are excluded).}
#'   \item{wbar}{Numeric. Average individual Wald statistic (variance
#'     normalisation \code{RSS / n}).}
#'   \item{wbar_dh}{Numeric. Average individual Wald statistic with the
#'     unbiased variance \code{RSS / (T_i - 2 K_i - 1)} (Dumitrescu and
#'     Hurlin 2012, eq. 4). \code{NA} unless \code{kmax = 0} and
#'     \code{dmax = 0}.}
#'   \item{zbar}{Numeric. Dumitrescu and Hurlin (2012) Z-bar statistic,
#'     eq. (9), with asymptotic moments. \code{NA} unless \code{kmax = 0}
#'     and \code{dmax = 0}.}
#'   \item{zbar_pv}{Numeric. One-sided (upper tail) asymptotic p-value of
#'     \code{zbar}.}
#'   \item{ztilde}{Numeric. Dumitrescu and Hurlin (2012) Z-tilde statistic,
#'     eq. (26), with the fixed-T moments of eq. (21) and (22). \code{NA}
#'     unless \code{kmax = 0} and \code{dmax = 0}.}
#'   \item{ztilde_pv}{Numeric. One-sided (upper tail) asymptotic p-value of
#'     \code{ztilde}.}
#'   \item{ztilde_boot_pv}{Numeric. One-sided bootstrap p-value of
#'     \code{ztilde} from the cross-sectional block bootstrap. \code{NA}
#'     unless \code{kmax = 0} and \code{dmax = 0}.}
#'   \item{ind_wald}{Numeric vector. Individual Wald statistics (length N;
#'     \code{NA} for a unit with a singular regression).}
#'   \item{ind_wald_dh}{Numeric vector. Individual Wald statistics with the
#'     unbiased variance (length N). \code{NA} unless \code{kmax = 0} and
#'     \code{dmax = 0}.}
#'   \item{ind_pval_a}{Numeric vector. Asymptotic chi-squared p-values of the
#'     individual Wald statistics (\code{ind_wald}) with \code{ind_lags}
#'     degrees of freedom (length N). Not used for the panel statistics; the
#'     bootstrap p-values are the basis of inference.}
#'   \item{ind_freq}{Integer vector. Selected Fourier frequencies (length N;
#'     0 when \code{kmax = 0}).}
#'   \item{ind_pval_b}{Numeric vector. Bootstrap p-values (length N).}
#'   \item{ind_lags}{Integer vector. Selected lag orders (length N).}
#'   \item{ind_df}{Integer vector. Residual degrees of freedom of the
#'     unrestricted regressions (length N).}
#'   For PQC:
#'   \item{quantiles}{Numeric vector. Quantiles tested.}
#'   \item{wald_xy}{Numeric vector. Wald statistics per quantile (x to y).}
#'   \item{pval_xy}{Numeric vector. Bootstrap p-values per quantile (x to y).}
#'   \item{wald_yx}{Numeric vector. Wald statistics per quantile (y to x).}
#'   \item{pval_yx}{Numeric vector. Bootstrap p-values per quantile (y to x).}
#'   \item{supwald_xy}{Numeric. Sup-Wald statistic for x to y.}
#'   \item{supwald_yx}{Numeric. Sup-Wald statistic for y to x.}
#'   \item{p_opt}{Integer. Selected optimal lag.}
#' }
#'
#' @references
#' Chuang, C.C., Kuan, C.M. and Lin, H.Y. (2009).
#' Causality in quantiles and dynamic stock return-volume relations.
#' \emph{Journal of Banking and Finance}, 33(7), 1351--1360.
#' \doi{10.1016/j.jbankfin.2009.02.013}
#'
#' Dumitrescu, E.I. and Hurlin, C. (2012).
#' Testing for Granger non-causality in heterogeneous panels.
#' \emph{Economic Modelling}, 29(4), 1450--1460.
#' \doi{10.1016/j.econmod.2012.02.014}
#'
#' Emirmahmutoglu, F. and Kose, N. (2011).
#' Testing for Granger causality in heterogeneous mixed panels.
#' \emph{Economic Modelling}, 28(3), 870--876.
#' \doi{10.1016/j.econmod.2010.10.018}
#'
#' Toda, H.Y. and Yamamoto, T. (1995).
#' Statistical inference in vector autoregressions with possibly integrated
#' processes. \emph{Journal of Econometrics}, 66(1--2), 225--250.
#' \doi{10.1016/0304-4076(94)01616-8}
#'
#' Wang, K.M. and Nguyen, T.B. (2022).
#' A quantile panel-type analysis of income inequality and healthcare
#' expenditure. \emph{Economic Research}, 35(1), 873--893.
#' \doi{10.1080/1331677X.2021.1948436}
#'
#' Yilanci, V. and Gorus, M.S. (2020).
#' Does economic globalization have predictive power for ecological footprint.
#' \emph{Environmental Science and Pollution Research}, 27, 40552--40562.
#' \doi{10.1007/s11356-020-10092-9}
#'
#' @examples
#' dat <- grunfeld_panel()
#' # PFTY test (quick with few bootstrap reps)
#' \donttest{
#' res <- xtpcaus(dat, y = "invest", x = "mvalue",
#'                panel_id = "firm", time_id = "year",
#'                test = "pfty", pmax = 2L, dmax = 1L,
#'                nboot = 99L, kmax = 2L, seed = 42L)
#' print(res)
#'
#' # Plain VAR(p) specification with the Dumitrescu and Hurlin (2012)
#' # Z-bar and Z-tilde statistics
#' res_dh <- xtpcaus(dat, y = "invest", x = "mvalue",
#'                   panel_id = "firm", time_id = "year",
#'                   test = "pfty", pmax = 2L, dmax = 0L,
#'                   nboot = 99L, kmax = 0L, seed = 42L)
#' print(res_dh)
#' }
#'
#' # PQC test
#' \donttest{
#' res2 <- xtpcaus(dat, y = "invest", x = "mvalue",
#'                 panel_id = "firm", time_id = "year",
#'                 test = "pqc", pmax = 2L, nboot = 99L,
#'                 quantiles = c(0.25, 0.50, 0.75), seed = 42L)
#' print(res2)
#' }
#'
#' @export
xtpcaus <- function(data, y, x, panel_id, time_id,
                    test = c("pfty", "pqc"),
                    pmax = 4L, dmax = 1L, nboot = 499L,
                    kmax = 3L, ic = c("aic", "bic"),
                    quantiles = c(0.1, 0.25, 0.50, 0.75, 0.90),
                    seed = -1L) {

  test <- match.arg(test)
  ic   <- match.arg(ic)

  # Validate inputs
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (v in c(y, x, panel_id, time_id)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found in data.", v))
  }
  if (nboot < 99L) stop("'nboot' must be at least 99.")
  if (any(quantiles <= 0) || any(quantiles >= 1)) {
    stop("All 'quantiles' values must be strictly between 0 and 1.")
  }

  pmax  <- as.integer(pmax)
  dmax  <- as.integer(dmax)
  nboot <- as.integer(nboot)
  kmax  <- as.integer(kmax)
  if (is.na(pmax) || pmax < 1L) stop("'pmax' must be a positive integer.")
  if (is.na(dmax) || dmax < 0L) stop("'dmax' must be a non-negative integer.")
  if (is.na(kmax) || kmax < 0L) {
    stop("'kmax' must be a non-negative integer (0 omits the Fourier terms).")
  }

  # Missing values: refuse, reporting count and positions (rows of 'data')
  for (v in c(y, x, panel_id, time_id)) {
    na_pos <- which(is.na(data[[v]]))
    if (length(na_pos)) {
      stop(sprintf(
        "Variable '%s' has %d missing value(s) at row(s) %s of 'data'. Complete cases are required.",
        v, length(na_pos), .xtpcaus_fmt_pos(na_pos)))
    }
  }

  # Sort and check the (panel_id, time_id) structure
  data <- data[order(data[[panel_id]], data[[time_id]]), , drop = FALSE]
  panels <- unique(data[[panel_id]])
  N <- length(panels)
  times <- sort(unique(data[[time_id]]))
  TT <- length(times)

  if (N < 2L) stop("At least 2 panel units are required.")
  if (TT < 10L) stop("At least 10 time periods are required.")

  time_by_unit <- split(data[[time_id]], data[[panel_id]])
  for (nm in names(time_by_unit)) {
    tu <- time_by_unit[[nm]]
    dup <- unique(tu[duplicated(tu)])
    if (length(dup)) {
      stop(sprintf(
        "Panel unit '%s' has duplicated time value(s): %s. Each (panel_id, time_id) pair must occur once.",
        nm, .xtpcaus_fmt_pos(dup)))
    }
    if (length(tu) != TT || !all(tu == times)) {
      stop(sprintf(
        "Panel must be strongly balanced: unit '%s' has %d observation(s) and the panel has %d time period(s). Every unit must be observed at the same %d time values.",
        nm, length(tu), TT, TT))
    }
  }

  # Local RNG scope: the caller's generator state is restored on exit
  if (!is.na(seed) && seed >= 0L) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }

  Y_mat <- matrix(data[[y]], nrow = TT, ncol = N)
  X_mat <- matrix(data[[x]], nrow = TT, ncol = N)

  if (test == "pfty") {
    res <- .xtpcaus_pfty(Y_mat, X_mat, N, TT, pmax, dmax, kmax, ic, nboot,
                         unit_labels = as.character(panels))
  } else {
    res <- .xtpcaus_pqc(Y_mat, X_mat, N, TT, pmax, dmax, quantiles, nboot)
  }

  res$test      <- test
  res$N         <- N
  res$TT        <- TT
  res$nboot     <- nboot
  res$y         <- y
  res$x         <- x
  res$panel_id  <- panel_id
  res$time_id   <- time_id

  structure(res, class = "xtpcaus")
}


# ============================================================
# PFTY: Panel Fourier Toda-Yamamoto
# ============================================================

# Minimum residual degrees of freedom of the unrestricted regression for a
# lag order to be admissible. It implies T_i > 5 + 2 K_i (Dumitrescu and
# Hurlin 2012, eq. 10) for the plain specification.
.xtpcaus_min_df <- 5L

#' @keywords internal
.xtpcaus_fmt_pos <- function(pos, max_show = 10L) {
  pos <- as.character(pos)
  if (length(pos) > max_show) {
    paste0(paste(pos[seq_len(max_show)], collapse = ", "), ", ... (",
           length(pos) - max_show, " more)")
  } else {
    paste(pos, collapse = ", ")
  }
}

#' @keywords internal
.xtpcaus_n_cols <- function(p, dmax, fourier) {
  1L + 2L * (p + dmax) + if (fourier) 2L else 0L
}

#' @keywords internal
.xtpcaus_admissible_p <- function(TT, pmax, dmax, kmax) {
  p <- seq_len(pmax)
  n_reg <- TT - (p + dmax)
  df <- n_reg - .xtpcaus_n_cols(p, dmax, kmax > 0L)
  list(p = p[df >= .xtpcaus_min_df], df = df)
}

#' @keywords internal
.xtpcaus_design <- function(y_vec, x_vec, TT, p, f, dmax) {
  k_aug <- p + dmax
  t_idx <- seq(k_aug + 1L, TT)
  n_reg <- length(t_idx)
  n_cols <- .xtpcaus_n_cols(p, dmax, f > 0L)
  # Regressors: const, lags 1..k_aug of y, lags 1..k_aug of x, [sin, cos]
  X_reg <- matrix(1, nrow = n_reg, ncol = n_cols)
  for (lag in seq_len(k_aug)) {
    X_reg[, 1L + lag]         <- y_vec[t_idx - lag]
    X_reg[, 1L + k_aug + lag] <- x_vec[t_idx - lag]
  }
  if (f > 0L) {
    X_reg[, 2L * k_aug + 2L] <- sin(2 * pi * f * t_idx / TT)
    X_reg[, 2L * k_aug + 3L] <- cos(2 * pi * f * t_idx / TT)
  }
  # H0: first p lags of x (cols 2 + k_aug to 1 + k_aug + p) = 0
  x_lag_cols <- seq(2L + k_aug, 1L + k_aug + p)
  list(X = X_reg, y = y_vec[t_idx], t_idx = t_idx, n_reg = n_reg,
       n_cols = n_cols, k_aug = k_aug, p = p, f = f, x_lag_cols = x_lag_cols,
       df = n_reg - n_cols)
}

#' @keywords internal
.xtpcaus_select_lag_freq <- function(y_vec, x_vec, TT, p_grid, kmax,
                                     dmax, ic_fn) {
  f_grid  <- if (kmax > 0L) seq_len(kmax) else 0L
  best_ic <- Inf
  best_p  <- NA_integer_
  best_f  <- NA_integer_

  for (p in p_grid) {
    for (f in f_grid) {
      d <- .xtpcaus_design(y_vec, x_vec, TT, p, f, dmax)
      beta <- tryCatch(
        solve(crossprod(d$X), crossprod(d$X, d$y)),
        error = function(e) NULL
      )
      if (is.null(beta)) next
      resid <- d$y - d$X %*% beta
      sigma2 <- sum(resid^2) / d$n_reg
      ic_val <- ic_fn(sigma2, d$n_cols, d$n_reg)
      if (is.finite(ic_val) && ic_val < best_ic) {
        best_ic <- ic_val
        best_p  <- p
        best_f  <- f
      }
    }
  }
  list(p = best_p, f = best_f)
}

#' @keywords internal
.xtpcaus_ic_fn <- function(ic) {
  if (ic == "aic") {
    function(sigma2, k, n) log(sigma2) + 2 * k / n
  } else {
    function(sigma2, k, n) log(sigma2) + log(n) * k / n
  }
}

# Wald statistic for H0: beta[x_lag_cols] = 0 on a fixed design.
# wald    : variance RSS / n_reg (package normalisation)
# wald_dh : variance RSS / (n_reg - n_cols) (Dumitrescu and Hurlin 2012,
#           eq. (5) with sigma2 = e'e / (T - 2K - 1); footnote 3 gives the
#           e'e / T variant, which is the wald above)
#' @keywords internal
.xtpcaus_wald_core <- function(y_dep, d) {
  XtX_inv <- tryCatch(solve(crossprod(d$X)), error = function(e) NULL)
  if (is.null(XtX_inv)) {
    return(list(wald = NA_real_, wald_dh = NA_real_, resid = NULL,
                beta = NULL, singular = TRUE))
  }
  beta  <- XtX_inv %*% crossprod(d$X, y_dep)
  resid <- y_dep - d$X %*% beta
  rss   <- sum(resid^2)
  p     <- d$p
  R_mat <- matrix(0, nrow = p, ncol = d$n_cols)
  for (jj in seq_len(p)) R_mat[jj, d$x_lag_cols[jj]] <- 1
  Rbeta <- R_mat %*% beta
  RVR   <- R_mat %*% XtX_inv %*% t(R_mat)
  core  <- as.numeric(t(Rbeta) %*% solve(RVR) %*% Rbeta)
  list(wald = core / (rss / d$n_reg),
       wald_dh = core / (rss / d$df),
       resid = as.numeric(resid), beta = as.numeric(beta), singular = FALSE)
}

#' @keywords internal
.xtpcaus_wald_stat <- function(y_vec, x_vec, TT, p, f, dmax) {
  d <- .xtpcaus_design(y_vec, x_vec, TT, p, f, dmax)
  if (d$df < .xtpcaus_min_df) {
    stop(sprintf(
      "Lag order p = %d with dmax = %d leaves %d residual degrees of freedom at T = %d; at least %d are required.",
      p, dmax, d$df, TT, .xtpcaus_min_df))
  }
  w <- .xtpcaus_wald_core(d$y, d)
  # Restricted model (tested x lags removed) for the bootstrap under H0
  if (!w$singular) {
    X0 <- d$X[, -d$x_lag_cols, drop = FALSE]
    beta0 <- tryCatch(solve(crossprod(X0), crossprod(X0, d$y)),
                      error = function(e) NULL)
    if (is.null(beta0)) {
      w <- list(wald = NA_real_, wald_dh = NA_real_, resid = NULL,
                beta = NULL, singular = TRUE)
    } else {
      beta0    <- as.numeric(beta0)
      w$fit0   <- as.numeric(X0 %*% beta0)
      w$resid0 <- as.numeric(d$y - w$fit0)
      # Restricted model split for the recursive bootstrap: autoregressive
      # coefficients on y lags 1..k_aug (columns 2..1 + k_aug of X0) and the
      # deterministic part (intercept, unrestricted x lags, Fourier terms)
      y_cols   <- seq(2L, 1L + d$k_aug)
      w$g0     <- beta0[y_cols]
      w$det0   <- as.numeric(X0[, -y_cols, drop = FALSE] %*% beta0[-y_cols])
    }
  }
  c(w, list(design = d, n_reg = d$n_reg, df = d$df, t_idx = d$t_idx,
            p = p, f = f, k_aug = d$k_aug, x_lag_cols = d$x_lag_cols))
}

# Bootstrap series under H0 for one unit: y~_t = det0_t + sum_k g0_k y~_{t-k}
# + e~_t for t = k_aug + 1, ..., T, with the observed first k_aug values as
# initial conditions (restricted model, tested x lags removed).
#' @keywords internal
.xtpcaus_boot_series <- function(wi, e_b, y_vec, TT) {
  k   <- wi$k_aug
  y_b <- y_vec
  for (t in seq(k + 1L, TT)) {
    y_b[t] <- wi$det0[t - k] + sum(wi$g0 * y_b[t - seq_len(k)]) + e_b[t - k]
  }
  y_b
}

# Dumitrescu and Hurlin (2012) moments of the individual Wald statistic
# (unbiased variance) for effective sample T_i and K_i restrictions:
# eq. (21) and (22), valid for T_i > 5 + 2 K_i (eq. 10).
#' @keywords internal
.xtpcaus_dh_moments <- function(T_i, K_i) {
  ok <- T_i > 5 + 2 * K_i
  EW <- VW <- rep(NA_real_, length(T_i))
  EW[ok] <- K_i[ok] * (T_i[ok] - 2 * K_i[ok] - 1) / (T_i[ok] - 2 * K_i[ok] - 3)
  VW[ok] <- 2 * K_i[ok] * (T_i[ok] - 2 * K_i[ok] - 1)^2 *
    (T_i[ok] - K_i[ok] - 3) /
    ((T_i[ok] - 2 * K_i[ok] - 3)^2 * (T_i[ok] - 2 * K_i[ok] - 5))
  list(EW = EW, VW = VW)
}

#' @keywords internal
.xtpcaus_pfty <- function(Y_mat, X_mat, N, TT, pmax, dmax, kmax, ic, nboot,
                          unit_labels = as.character(seq_len(N))) {
  ic_fn <- .xtpcaus_ic_fn(ic)
  fourier <- kmax > 0L
  plain_dh <- (kmax == 0L) && (dmax == 0L)

  # Admissible lag orders: residual df >= .xtpcaus_min_df
  adm <- .xtpcaus_admissible_p(TT, pmax, dmax, kmax)
  if (!length(adm$p)) {
    stop(sprintf(paste0(
      "No admissible lag order: with T = %d, dmax = %d and %s the unrestricted ",
      "regression with p = 1 has %d residual degrees of freedom ",
      "(T - (p + dmax) observations, 2 (p + dmax) + %d regressors); at least %d are required. ",
      "Increase T or reduce dmax."),
      TT, dmax, if (fourier) "Fourier terms" else "no Fourier terms",
      adm$df[1L], if (fourier) 3L else 1L, .xtpcaus_min_df))
  }
  if (length(adm$p) < pmax) {
    warning(sprintf(paste0(
      "Lag orders above %d are not considered: at T = %d with dmax = %d they leave ",
      "fewer than %d residual degrees of freedom in the unrestricted regression."),
      max(adm$p), TT, dmax, .xtpcaus_min_df), call. = FALSE)
  }

  ind_wald    <- rep(NA_real_, N)
  ind_wald_dh <- rep(NA_real_, N)
  ind_freq    <- rep(NA_integer_, N)
  ind_lags    <- rep(NA_integer_, N)
  ind_df      <- rep(NA_integer_, N)
  ind_pval_a  <- rep(NA_real_, N)
  ind_pval_b  <- rep(NA_real_, N)
  wald_info   <- vector("list", N)

  for (i in seq_len(N)) {
    sel <- .xtpcaus_select_lag_freq(Y_mat[, i], X_mat[, i],
                                    TT, adm$p, kmax, dmax, ic_fn)
    if (is.na(sel$p)) next           # singular for every candidate
    wi <- .xtpcaus_wald_stat(Y_mat[, i], X_mat[, i], TT, sel$p, sel$f, dmax)
    if (wi$singular) next
    ind_wald[i]    <- wi$wald
    ind_wald_dh[i] <- wi$wald_dh
    ind_freq[i]    <- sel$f
    ind_lags[i]    <- sel$p
    ind_df[i]      <- wi$df
    ind_pval_a[i]  <- stats::pchisq(wi$wald, df = sel$p, lower.tail = FALSE)
    wald_info[[i]] <- wi
  }

  ok <- !is.na(ind_wald)
  N_used <- sum(ok)
  if (N_used < N) {
    warning(sprintf(paste0(
      "Unit(s) %s: singular regressor matrix (X'X not invertible, for example x ",
      "constant or collinear with the lags of y); the individual Wald statistic is NA ",
      "and the unit is excluded from the panel statistics (N used = %d)."),
      paste(unit_labels[!ok], collapse = ", "), N_used), call. = FALSE)
  }

  # Dumitrescu and Hurlin (2012) moments (plain specification only)
  dh_EW <- dh_VW <- rep(NA_real_, N)
  if (plain_dh && N_used > 0L) {
    T_i <- TT - ind_lags            # effective sample of unit i
    mom <- .xtpcaus_dh_moments(T_i[ok], ind_lags[ok])
    dh_EW[ok] <- mom$EW
    dh_VW[ok] <- mom$VW
  }
  dh_ztilde <- function(w_dh, K_vec) {
    mom <- .xtpcaus_dh_moments(TT - K_vec, K_vec)
    sqrt(length(w_dh)) * (mean(w_dh) - mean(mom$EW)) / sqrt(mean(mom$VW))
  }

  # Bootstrap under H0 (Dumitrescu and Hurlin 2012, Section 6.2, steps 3-7):
  # residuals of the restricted model; the same time indices are drawn for
  # all units (block of size 1 in time and size N in the cross-section); the
  # bootstrap series is generated recursively from the restricted model; the
  # lag and frequency selection and the Wald statistic are then repeated on
  # the bootstrap series, so that the selection step is part of the bootstrap
  # distribution. With a single candidate (pmax = 1, kmax <= 1) the
  # selection is skipped.
  reselect     <- length(adm$p) > 1L || kmax > 1L
  boot_wald    <- matrix(NA_real_, nrow = nboot, ncol = N)
  boot_wald_dh <- matrix(NA_real_, nrow = nboot, ncol = N)
  boot_lags    <- matrix(NA_integer_, nrow = nboot, ncol = N)
  boot_ztilde  <- rep(NA_real_, nboot)
  if (N_used > 0L) {
    k_aug_i  <- vapply(which(ok), function(i) wald_info[[i]]$k_aug, integer(1))
    n_reg_i  <- vapply(which(ok), function(i) wald_info[[i]]$n_reg, integer(1))
    t_common <- seq(max(k_aug_i) + 1L, TT)   # periods with residuals in every unit
    n_common <- length(t_common)
    n_draw   <- max(n_reg_i)
    for (b in seq_len(nboot)) {
      t_b <- t_common[sample.int(n_common, n_draw, replace = TRUE)]
      for (i in which(ok)) {
        wi  <- wald_info[[i]]
        # residual at time t sits at position t - k_aug in resid0
        e_b <- wi$resid0[t_b[seq_len(wi$n_reg)] - wi$k_aug]
        y_b <- .xtpcaus_boot_series(wi, e_b, Y_mat[, i], TT)
        if (!all(is.finite(y_b))) next
        if (reselect) {
          sel_b <- .xtpcaus_select_lag_freq(y_b, X_mat[, i], TT, adm$p,
                                            kmax, dmax, ic_fn)
          if (is.na(sel_b$p)) next
        } else {
          sel_b <- list(p = wi$p, f = wi$f)
        }
        d_b <- .xtpcaus_design(y_b, X_mat[, i], TT, sel_b$p, sel_b$f, dmax)
        wb  <- .xtpcaus_wald_core(d_b$y, d_b)
        if (!is.finite(wb$wald)) next
        boot_wald[b, i]    <- wb$wald
        boot_wald_dh[b, i] <- wb$wald_dh
        boot_lags[b, i]    <- sel_b$p
      }
      if (plain_dh && all(!is.na(boot_wald_dh[b, ok]))) {
        boot_ztilde[b] <- dh_ztilde(boot_wald_dh[b, ok], boot_lags[b, ok])
      }
    }
    n_boot_na <- sum(is.na(boot_wald[, ok]))
    if (n_boot_na > 0L) {
      warning(sprintf(paste0(
        "%d of %d bootstrap Wald statistics could not be computed (singular or ",
        "non-finite bootstrap regression) and were dropped from the bootstrap p-values."),
        n_boot_na, nboot * N_used), call. = FALSE)
    }
  }

  for (i in which(ok)) {
    ind_pval_b[i] <- mean(boot_wald[, i] >= ind_wald[i], na.rm = TRUE)
    # Avoid exact 0
    ind_pval_b[i] <- max(ind_pval_b[i], 1 / (nboot + 1))
  }

  # Fisher panel statistic (Emirmahmutoglu and Kose 2011)
  if (N_used > 0L) {
    fisher_stat <- -2 * sum(log(ind_pval_b[ok]))
    fisher_df   <- 2L * N_used
    fisher_pv   <- stats::pchisq(fisher_stat, df = fisher_df, lower.tail = FALSE)
    wbar        <- mean(ind_wald[ok])
  } else {
    fisher_stat <- fisher_pv <- wbar <- NA_real_
    fisher_df   <- 0L
  }

  # Dumitrescu and Hurlin (2012) standardised statistics, plain
  # specification only: Z-bar (eq. 9) and Z-tilde (eq. 26), with the moments
  # averaged over units when K_i differs (Z-bar: eq. 11; Z-tilde: eq. 33);
  # one-sided p-values.
  wbar_dh <- zbar <- zbar_pv <- ztilde <- ztilde_pv <- ztilde_boot_pv <- NA_real_
  if (plain_dh && N_used > 0L) {
    K_i     <- ind_lags[ok]
    wbar_dh <- mean(ind_wald_dh[ok])
    zbar    <- sqrt(N_used) * (wbar_dh - mean(K_i)) / sqrt(mean(2 * K_i))
    zbar_pv <- stats::pnorm(zbar, lower.tail = FALSE)
    if (all(!is.na(dh_VW[ok]))) {
      ztilde    <- dh_ztilde(ind_wald_dh[ok], ind_lags[ok])
      ztilde_pv <- stats::pnorm(ztilde, lower.tail = FALSE)
      zt_b <- boot_ztilde[!is.na(boot_ztilde)]
      if (length(zt_b)) {
        ztilde_boot_pv <- max(mean(zt_b >= ztilde), 1 / (nboot + 1))
      }
    } else {
      warning(sprintf(
        "Z-tilde not computed: T_i > 5 + 2 K_i (Dumitrescu and Hurlin 2012, eq. 10) fails for unit(s) %s.",
        paste(unit_labels[ok][is.na(dh_VW[ok])], collapse = ", ")), call. = FALSE)
    }
  } else if (!plain_dh) {
    ind_wald_dh[] <- NA_real_
  }

  list(
    fisher = fisher_stat, fisher_df = fisher_df, fisher_pv = fisher_pv,
    N_used = N_used,
    wbar = wbar, wbar_dh = wbar_dh,
    zbar = zbar, zbar_pv = zbar_pv,
    ztilde = ztilde, ztilde_pv = ztilde_pv, ztilde_boot_pv = ztilde_boot_pv,
    ind_wald = ind_wald, ind_wald_dh = ind_wald_dh, ind_freq = ind_freq,
    ind_pval_a = ind_pval_a, ind_pval_b = ind_pval_b,
    ind_lags = ind_lags, ind_df = ind_df,
    pmax = pmax, dmax = dmax, kmax = kmax, ic = ic,
    plain_dh = plain_dh, min_df = .xtpcaus_min_df
  )
}


# ============================================================
# PQC: Panel Quantile Causality
# ============================================================

#' @keywords internal
.xtpcaus_qreg_coef <- function(y, X, tau) {
  # Barrodale-Roberts quantile regression (L1)
  n  <- length(y)
  k  <- ncol(X)
  # Use stats::optim with Koenker-Bassett loss for small problems
  # For efficiency, implement direct iteratively-reweighted approach
  # Start from OLS
  beta0 <- tryCatch(
    as.numeric(solve(crossprod(X), crossprod(X, y))),
    error = function(e) rep(0, k)
  )
  # IRLS for quantile regression
  for (iter in seq_len(50L)) {
    resid <- as.numeric(y - X %*% beta0)
    w <- pmax(abs(resid), 1e-6)
    W <- diag(1 / w)
    # Tilted weight
    wt <- ifelse(resid >= 0, tau / w, (1 - tau) / w)
    beta1 <- tryCatch(
      as.numeric(solve(t(X) %*% W %*% X, t(X) %*% (W %*% y + wt * w - wt * resid))),
      error = function(e) beta0
    )
    if (max(abs(beta1 - beta0)) < 1e-8) { beta0 <- beta1; break }
    beta0 <- beta1
  }
  beta0
}

#' @keywords internal
.xtpcaus_pqc_wald <- function(y_vec, x_vec, N, TT, p, tau) {
  # Panel quantile VAR with individual fixed effects (within transformation)
  # Stack: outcome = y, regressors = lags of y + lags of x + individual dummies
  n_reg <- TT - p
  if (n_reg < p + 2L) return(list(wald = 0, coef_sum = 0))

  t_idx <- seq(p + 1L, TT)

  # Build stacked (N * n_reg) x (2p + N) design matrix
  n_obs <- N * n_reg
  n_cols <- 2L * p + N  # p lags y + p lags x + N fixed effects
  X_all <- matrix(0, nrow = n_obs, ncol = n_cols)
  y_all <- numeric(n_obs)

  for (i in seq_len(N)) {
    rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
    y_all[rows_i] <- y_vec[t_idx + (i - 1L) * TT]
    for (lag in seq_len(p)) {
      X_all[rows_i, lag]     <- y_vec[t_idx - lag + (i - 1L) * TT]
      X_all[rows_i, p + lag] <- x_vec[t_idx - lag + (i - 1L) * TT]
    }
    X_all[rows_i, 2L * p + i] <- 1  # panel fixed effect
  }

  beta_q <- tryCatch(
    .xtpcaus_qreg_coef(y_all, X_all, tau),
    error = function(e) NULL
  )
  if (is.null(beta_q)) return(list(wald = 0, coef_sum = 0))

  # Extract x-lag coefficients (cols p+1 to 2p)
  x_lag_cols <- seq(p + 1L, 2L * p)
  beta_x <- beta_q[x_lag_cols]

  # Wald test: H0: beta_x = 0
  resid <- y_all - X_all %*% beta_q
  sigma2 <- sum(resid^2) / (n_obs - ncol(X_all))
  XtX_inv <- tryCatch(solve(crossprod(X_all)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(wald = 0, coef_sum = sum(beta_x)))

  R_mat <- matrix(0, nrow = p, ncol = n_cols)
  for (jj in seq_len(p)) R_mat[jj, x_lag_cols[jj]] <- 1
  RVR   <- R_mat %*% XtX_inv %*% t(R_mat)
  wald  <- tryCatch(
    as.numeric(t(beta_x) %*% solve(sigma2 * RVR) %*% beta_x),
    error = function(e) 0
  )
  list(wald = max(wald, 0), coef_sum = sum(beta_x))
}

#' @keywords internal
.xtpcaus_pqc <- function(Y_mat, X_mat, N, TT, pmax, dmax, quantiles, nboot) {
  # Select optimal lag via AIC (OLS)
  y_vec <- as.vector(Y_mat)
  x_vec <- as.vector(X_mat)

  best_aic <- Inf
  p_opt <- 1L
  for (p in seq_len(pmax)) {
    n_reg <- TT - p
    if (n_reg < p + 2L) next
    t_idx <- seq(p + 1L, TT)
    n_obs_p <- N * n_reg
    X_ls <- matrix(0, nrow = n_obs_p, ncol = 2L * p + N)
    y_ls <- numeric(n_obs_p)
    for (i in seq_len(N)) {
      rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
      y_ls[rows_i] <- Y_mat[t_idx, i]
      for (lag in seq_len(p)) {
        X_ls[rows_i, lag]     <- Y_mat[t_idx - lag, i]
        X_ls[rows_i, p + lag] <- X_mat[t_idx - lag, i]
      }
      X_ls[rows_i, 2L * p + i] <- 1
    }
    beta_ols <- tryCatch(
      solve(crossprod(X_ls), crossprod(X_ls, y_ls)),
      error = function(e) NULL
    )
    if (is.null(beta_ols)) next
    rss <- sum((y_ls - X_ls %*% beta_ols)^2)
    aic_val <- log(rss / n_obs_p) + 2 * ncol(X_ls) / n_obs_p
    if (aic_val < best_aic) { best_aic <- aic_val; p_opt <- p }
  }

  nq <- length(quantiles)
  wald_xy  <- numeric(nq)
  coef_xy  <- numeric(nq)
  wald_yx  <- numeric(nq)
  coef_yx  <- numeric(nq)
  pval_xy  <- numeric(nq)
  pval_yx  <- numeric(nq)

  # Vectorized Y and X for stacking
  y_stacked <- as.vector(Y_mat)
  x_stacked <- as.vector(X_mat)

  for (q_idx in seq_len(nq)) {
    tau <- quantiles[q_idx]
    res_xy <- .xtpcaus_pqc_wald(y_stacked, x_stacked, N, TT, p_opt, tau)
    res_yx <- .xtpcaus_pqc_wald(x_stacked, y_stacked, N, TT, p_opt, tau)
    wald_xy[q_idx] <- res_xy$wald
    coef_xy[q_idx] <- res_xy$coef_sum
    wald_yx[q_idx] <- res_yx$wald
    coef_yx[q_idx] <- res_yx$coef_sum
  }

  # Bootstrap p-values
  n_reg <- TT - p_opt
  t_idx_b <- seq(p_opt + 1L, TT)
  n_obs_b <- N * n_reg
  X_bs <- matrix(0, nrow = n_obs_b, ncol = 2L * p_opt + N)
  y_bs <- numeric(n_obs_b)
  for (i in seq_len(N)) {
    rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
    y_bs[rows_i] <- Y_mat[t_idx_b, i]
    for (lag in seq_len(p_opt)) {
      X_bs[rows_i, lag]           <- Y_mat[t_idx_b - lag, i]
      X_bs[rows_i, p_opt + lag]   <- X_mat[t_idx_b - lag, i]
    }
    X_bs[rows_i, 2L * p_opt + i] <- 1
  }

  # OLS under H0 (restrict x-lag coefs to 0)
  beta_ols_full <- tryCatch(
    as.numeric(solve(crossprod(X_bs), crossprod(X_bs, y_bs))),
    error = function(e) rep(0, ncol(X_bs))
  )
  beta_null <- beta_ols_full
  x_cols_null <- seq(p_opt + 1L, 2L * p_opt)
  beta_null[x_cols_null] <- 0
  resid_null <- y_bs - X_bs %*% beta_null

  # Bootstrap
  boot_wald_xy <- matrix(0, nrow = nboot, ncol = nq)
  boot_wald_yx <- matrix(0, nrow = nboot, ncol = nq)

  for (b in seq_len(nboot)) {
    idx_b <- sample.int(n_obs_b, n_obs_b, replace = TRUE)
    e_b   <- resid_null[idx_b]
    y_boot <- X_bs %*% beta_null + e_b
    y_boot_vec <- as.numeric(y_boot)

    for (q_idx in seq_len(nq)) {
      tau <- quantiles[q_idx]
      wb <- tryCatch(
        .xtpcaus_pqc_wald(y_boot_vec, x_stacked, N, TT, p_opt, tau)$wald,
        error = function(e) 0
      )
      boot_wald_xy[b, q_idx] <- if (is.na(wb)) 0 else wb
      # For the y to x direction, bootstrap uses same y (x causes y direction uses x as dependent)
      wb2 <- tryCatch(
        .xtpcaus_pqc_wald(x_stacked, y_boot_vec, N, TT, p_opt, tau)$wald,
        error = function(e) 0
      )
      boot_wald_yx[b, q_idx] <- if (is.na(wb2)) 0 else wb2
    }
  }

  for (q_idx in seq_len(nq)) {
    pval_xy[q_idx] <- max(mean(boot_wald_xy[, q_idx] >= wald_xy[q_idx]),
                          1 / (nboot + 1))
    pval_yx[q_idx] <- max(mean(boot_wald_yx[, q_idx] >= wald_yx[q_idx]),
                          1 / (nboot + 1))
  }

  supwald_xy <- max(wald_xy)
  supwald_yx <- max(wald_yx)

  list(
    quantiles = quantiles, p_opt = p_opt,
    wald_xy = wald_xy, coef_xy = coef_xy, pval_xy = pval_xy,
    wald_yx = wald_yx, coef_yx = coef_yx, pval_yx = pval_yx,
    supwald_xy = supwald_xy, supwald_yx = supwald_yx
  )
}


# ============================================================
# S3 Methods
# ============================================================

#' Print Method for xtpcaus Objects
#'
#' @param x An object of class \code{"xtpcaus"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtpcaus <- function(x, ...) {
  message(strrep("-", 65))
  message(sprintf("Panel Causality Test: %s", toupper(x$test)))
  message(sprintf("H0: %s does NOT Granger-cause %s", x$x, x$y))
  message(sprintf("N = %d panels, T = %d periods, nboot = %d",
    x$N, x$TT, x$nboot))
  message(strrep("-", 65))

  stars_of <- function(pv) {
    if (is.na(pv)) "" else if (pv < 0.01) "***" else if (pv < 0.05) "**" else
      if (pv < 0.10) "*" else ""
  }
  fmt <- function(v) if (is.na(v)) "NA" else sprintf("%.4f", v)

  if (x$test == "pfty") {
    if (!is.null(x$N_used) && x$N_used < x$N) {
      message(sprintf("Units entering the panel statistics: %d of %d (singular regressions excluded)",
        x$N_used, x$N))
    }
    message(sprintf("Fisher stat: %s  df: %d  p-value: %s%s",
      fmt(x$fisher), x$fisher_df, fmt(x$fisher_pv), stars_of(x$fisher_pv)))
    message(sprintf("W-bar:       %s", fmt(x$wbar)))
    if (isTRUE(x$plain_dh)) {
      message("Dumitrescu and Hurlin (2012) statistics (plain VAR(p) specification,")
      message("  unbiased residual variance, one-sided p-values):")
      message(sprintf("W-bar (DH):  %s", fmt(x$wbar_dh)))
      message(sprintf("Z-bar:       %s  p-value: %s%s  (asymptotic in T, eq. 9)",
        fmt(x$zbar), fmt(x$zbar_pv), stars_of(x$zbar_pv)))
      message(sprintf("Z-tilde:     %s  p-value: %s%s  (fixed T, eq. 26)",
        fmt(x$ztilde), fmt(x$ztilde_pv), stars_of(x$ztilde_pv)))
      message(sprintf("             bootstrap p-value: %s%s",
        fmt(x$ztilde_boot_pv), stars_of(x$ztilde_boot_pv)))
    } else {
      message("Z-bar, Z-tilde: not reported for the Fourier or Toda-Yamamoto augmented")
      message("  specification (kmax > 0 or dmax > 0); the Dumitrescu and Hurlin (2012)")
      message("  moments do not apply. Inference rests on the bootstrap Fisher p-value.")
    }
    message("")
    message("Individual results:")
    message(sprintf("  %-6s %-6s %-8s %-6s %-12s %-12s",
      "Panel", "Lag", "Freq", "df", "Wald", "Boot p-val"))
    for (i in seq_len(x$N)) {
      if (is.na(x$ind_wald[i])) {
        message(sprintf("  %-6d %-6s %-8s %-6s %12s %8s  (singular regression)",
          i, "NA", "NA", "NA", "NA", "NA"))
      } else {
        message(sprintf("  %-6d %-6d %-8d %-6d %12.4f %8.4f%s",
          i, x$ind_lags[i], x$ind_freq[i], x$ind_df[i],
          x$ind_wald[i], x$ind_pval_b[i], stars_of(x$ind_pval_b[i])))
      }
    }
  } else {
    message(sprintf("Sup-Wald (%s to %s): %.4f", x$x, x$y, x$supwald_xy))
    message(sprintf("Sup-Wald (%s to %s): %.4f", x$y, x$x, x$supwald_yx))
    message(sprintf("Optimal lag: %d", x$p_opt))
    message("")
    message(sprintf("  %-6s %-12s %-10s %-12s %-10s",
      "Tau", paste0("Wald(", x$x, " to ", x$y, ")"),
      "p-val", paste0("Wald(", x$y, " to ", x$x, ")"), "p-val"))
    for (q_idx in seq_along(x$quantiles)) {
      message(sprintf("  %-6.2f %12.4f %10.4f%s %12.4f %10.4f%s",
        x$quantiles[q_idx],
        x$wald_xy[q_idx], x$pval_xy[q_idx],
        if (x$pval_xy[q_idx] < 0.05) "*" else " ",
        x$wald_yx[q_idx], x$pval_yx[q_idx],
        if (x$pval_yx[q_idx] < 0.05) "*" else " "))
    }
  }
  message(strrep("-", 65))
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}

#' Summary Method for xtpcaus Objects
#'
#' @param object An object of class \code{"xtpcaus"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtpcaus <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Example Panel Data for xtpcaus
#'
#' Returns a small balanced panel dataset (subset of Grunfeld 1958) for
#' use in examples and testing.
#'
#' @return A data frame with columns \code{firm}, \code{year},
#'   \code{invest}, and \code{mvalue}.
#'
#' @examples
#' dat <- grunfeld_panel()
#' head(dat)
#'
#' @export
grunfeld_panel <- function() {
  data.frame(
    firm = rep(1:5, each = 20L),
    year = rep(1935L:1954L, times = 5L),
    invest = c(
      317.6, 391.8, 410.6, 257.7, 330.8, 461.2, 512.0, 448.0, 499.6, 547.5,
      561.2, 688.1, 568.9, 529.2, 555.1, 642.9, 755.9, 891.2, 1304.4, 1486.7,
      40.29, 72.76, 66.26, 65.67, 76.04, 87.99, 100.0, 95.1, 104.4, 118.2,
      114.0, 135.7, 130.7, 169.6, 162.7, 162.0, 190.2, 181.9, 232.8, 256.7,
      209.9, 355.3, 318.9, 267.4, 339.5, 400.3, 419.5, 417.3, 347.2, 364.2,
      361.1, 312.2, 273.1, 264.0, 187.0, 145.8, 208.4, 162.8, 184.8, 152.1,
      33.1, 45.0, 77.2, 44.6, 48.1, 74.4, 113.0, 91.9, 61.3, 56.8,
      93.6, 159.9, 147.2, 146.3, 98.3, 93.5, 135.2, 157.3, 179.5, 189.6,
      40.29, 72.76, 66.26, 65.67, 76.04, 87.99, 100.0, 95.1, 104.4, 118.2,
      114.0, 135.7, 130.7, 169.6, 162.7, 162.0, 190.2, 181.9, 232.8, 256.7
    ),
    mvalue = c(
      2792.7, 2759.9, 2132.0, 1834.1, 1588.0, 1749.4, 1687.2, 2007.7, 2177.3, 2011.6,
      2533.2, 2459.4, 2157.7, 2082.0, 1808.1, 1786.0, 1975.7, 2021.7, 2390.4, 2613.6,
      182.8, 213.3, 206.6, 209.0, 225.5, 226.7, 226.7, 222.3, 227.9, 237.9,
      226.7, 232.6, 262.1, 302.1, 324.3, 335.9, 351.9, 344.9, 365.5, 382.1,
      1362.4, 1807.1, 1952.1, 2080.0, 1968.1, 1795.5, 1666.6, 1634.8, 1640.0, 1671.7,
      1547.5, 1560.0, 1401.4, 1364.0, 1216.8, 1099.8, 1131.3, 1169.7, 1289.5, 1467.1,
      95.3, 107.5, 136.3, 113.4, 141.1, 168.6, 214.8, 218.0, 186.0, 230.3,
      228.7, 293.2, 289.0, 271.7, 246.0, 239.6, 281.0, 265.3, 270.8, 290.5,
      182.8, 213.3, 206.6, 209.0, 225.5, 226.7, 226.7, 222.3, 227.9, 237.9,
      226.7, 232.6, 262.1, 302.1, 324.3, 335.9, 351.9, 344.9, 365.5, 382.1
    ),
    stringsAsFactors = FALSE
  )
}
