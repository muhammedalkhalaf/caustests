# Hand recomputations of xtpcaus(test = "pfty") against lm() and the
# formulas of Dumitrescu and Hurlin (2012); input validation; RNG scope.
# Tolerance 1e-10: the quantities are closed-form OLS functions of the same
# data and the two routes differ only by floating point rounding (observed
# differences of order 1e-14 in the audit).

sim_panel <- function(N, TT, seed, burn = 100L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    g <- runif(1, -0.9, 0.9); a <- rnorm(1); s <- sqrt(runif(1, 0.5, 1.5))
    y <- numeric(TT + burn); x <- numeric(TT + burn)
    for (t in 2:(TT + burn)) {
      y[t] <- a + g * y[t - 1] + s * rnorm(1)
      x[t] <- 0.5 * x[t - 1] + rnorm(1)
    }
    data.frame(id = i, time = seq_len(TT),
               y = y[(burn + 1):(TT + burn)], x = x[(burn + 1):(TT + burn)])
  }))
}

# lm() version of the PFTY regression: selection by AIC over p and f,
# Wald on the first p lags of x with variance RSS / n (package) and
# RSS / (n - k) (Dumitrescu and Hurlin, unbiased).
hand_unit <- function(y, x, TT, pmax, dmax, kmax, min_df = 5L) {
  best <- list(ic = Inf)
  f_grid <- if (kmax > 0) seq_len(kmax) else 0L
  for (p in seq_len(pmax)) for (f in f_grid) {
    k <- p + dmax; t_idx <- (k + 1):TT; n <- length(t_idx)
    ncols <- 1 + 2 * k + if (f > 0) 2 else 0
    if (n - ncols < min_df) next
    df <- data.frame(yt = y[t_idx])
    for (l in seq_len(k)) {
      df[[paste0("yl", l)]] <- y[t_idx - l]
      df[[paste0("xl", l)]] <- x[t_idx - l]
    }
    if (f > 0) {
      df$s <- sin(2 * pi * f * t_idx / TT)
      df$c <- cos(2 * pi * f * t_idx / TT)
    }
    fit <- lm(yt ~ ., data = df)
    if (any(is.na(coef(fit)))) next
    ic <- log(sum(resid(fit)^2) / n) + 2 * length(coef(fit)) / n
    if (ic < best$ic) best <- list(ic = ic, p = p, f = f, fit = fit, n = n)
  }
  fit <- best$fit; p <- best$p; n <- best$n
  X <- model.matrix(fit); XtXi <- solve(crossprod(X))
  idx <- match(paste0("xl", seq_len(p)), colnames(X)); b <- coef(fit)[idx]
  core <- as.numeric(t(b) %*% solve(XtXi[idx, idx, drop = FALSE]) %*% b)
  rss <- sum(resid(fit)^2)
  list(W = core / (rss / n), W_dh = core / (rss / (n - ncol(X))),
       p = p, f = best$f, df = n - ncol(X))
}

test_that("PFTY W_i, W_bar and Fisher match lm() (Fourier, dmax = 1)", {
  N <- 6L; TT <- 30L
  dat <- sim_panel(N, TT, 20261030)
  res <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty",
                 pmax = 2L, dmax = 1L, kmax = 3L, nboot = 99L, seed = 1L)
  hw <- lapply(seq_len(N), function(i) {
    d <- dat[dat$id == i, ]; hand_unit(d$y, d$x, TT, 2L, 1L, 3L)
  })
  expect_equal(res$ind_wald, sapply(hw, `[[`, "W"), tolerance = 1e-10)
  expect_equal(res$ind_lags, sapply(hw, `[[`, "p"))
  expect_equal(res$ind_freq, sapply(hw, `[[`, "f"))
  expect_equal(res$ind_df, sapply(hw, `[[`, "df"))
  expect_equal(res$wbar, mean(sapply(hw, `[[`, "W")), tolerance = 1e-10)
  expect_equal(res$fisher, -2 * sum(log(res$ind_pval_b)), tolerance = 1e-10)
  expect_equal(res$fisher_pv, pchisq(res$fisher, 2 * N, lower.tail = FALSE),
               tolerance = 1e-10)
  # DH statistics are not reported for this specification
  expect_true(is.na(res$zbar) && is.na(res$zbar_pv) && is.na(res$ztilde))
  expect_true(all(is.na(res$ind_wald_dh)))
  expect_equal(res$N_used, N)
})

test_that("plain specification: Z_bar (eq. 9) and Z_tilde (eq. 21, 22, 26) by hand", {
  N <- 8L; TT <- 25L
  dat <- sim_panel(N, TT, 777)
  res <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty",
                 pmax = 2L, dmax = 0L, kmax = 0L, nboot = 99L, seed = 3L)
  hw <- lapply(seq_len(N), function(i) {
    d <- dat[dat$id == i, ]; hand_unit(d$y, d$x, TT, 2L, 0L, 0L)
  })
  W_dh <- sapply(hw, `[[`, "W_dh"); K_i <- sapply(hw, `[[`, "p")
  expect_equal(res$ind_wald_dh, W_dh, tolerance = 1e-10)
  expect_true(all(res$ind_freq == 0L))
  T_i <- TT - K_i                                   # effective sample
  wbar_dh <- mean(W_dh)
  # eq. (9), eq. (33) form with asymptotic moments K_i, 2 K_i
  zbar <- sqrt(N) * (wbar_dh - mean(K_i)) / sqrt(mean(2 * K_i))
  # eq. (21), (22)
  EW <- K_i * (T_i - 2 * K_i - 1) / (T_i - 2 * K_i - 3)
  VW <- 2 * K_i * (T_i - 2 * K_i - 1)^2 * (T_i - K_i - 3) /
    ((T_i - 2 * K_i - 3)^2 * (T_i - 2 * K_i - 5))
  ztilde <- sqrt(N) * (wbar_dh - mean(EW)) / sqrt(mean(VW))   # eq. (25)/(33)
  expect_equal(res$wbar_dh, wbar_dh, tolerance = 1e-10)
  expect_equal(res$zbar, zbar, tolerance = 1e-10)
  expect_equal(res$ztilde, ztilde, tolerance = 1e-10)
  # one-sided (upper tail) p-values, the decision rule of the paper
  expect_equal(res$zbar_pv, pnorm(zbar, lower.tail = FALSE), tolerance = 1e-10)
  expect_equal(res$ztilde_pv, pnorm(ztilde, lower.tail = FALSE), tolerance = 1e-10)
  expect_true(res$ztilde_boot_pv >= 1 / 100 && res$ztilde_boot_pv <= 1)
  # when all K_i are equal, eq. (26) closed form equals eq. (25)
  if (length(unique(K_i)) == 1L) {
    K <- K_i[1]; Tt <- T_i[1]
    z26 <- sqrt(N / (2 * K) * (Tt - 2 * K - 5) / (Tt - K - 3)) *
      ((Tt - 2 * K - 3) / (Tt - 2 * K - 1) * wbar_dh - K)
    expect_equal(res$ztilde, z26, tolerance = 1e-10)
  }
  # paper's condition T_i > 5 + 2 K_i holds for every unit
  expect_true(all(T_i > 5 + 2 * K_i))
})

test_that("plain specification with pmax = 1: eq. (26) closed form", {
  N <- 5L; TT <- 20L
  dat <- sim_panel(N, TT, 4242)
  res <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty",
                 pmax = 1L, dmax = 0L, kmax = 0L, nboot = 99L, seed = 3L)
  K <- 1; Tt <- TT - 1
  z26 <- sqrt(N / (2 * K) * (Tt - 2 * K - 5) / (Tt - K - 3)) *
    ((Tt - 2 * K - 3) / (Tt - 2 * K - 1) * res$wbar_dh - K)
  expect_equal(res$ztilde, z26, tolerance = 1e-10)
  expect_equal(res$zbar, sqrt(N / (2 * K)) * (res$wbar_dh - K), tolerance = 1e-10)
})

test_that("degrees-of-freedom guard: T = 12 with pmax = 2 stops; lags above the bound warn", {
  dat <- sim_panel(4L, 12L, 11)
  expect_error(
    xtpcaus(dat, "y", "x", "id", "time", test = "pfty",
            pmax = 2L, dmax = 1L, kmax = 3L, nboot = 99L),
    "No admissible lag order")
  # T = 20, pmax = 4, dmax = 1: p = 4 leaves 15 - 13 = 2 df, p = 3 leaves 5
  dat20 <- sim_panel(4L, 20L, 12)
  expect_warning(
    r <- xtpcaus(dat20, "y", "x", "id", "time", test = "pfty",
                 pmax = 4L, dmax = 1L, kmax = 3L, nboot = 99L, seed = 1L),
    "Lag orders above 3 are not considered")
  expect_true(all(r$ind_lags <= 3L))
  expect_true(all(r$ind_df >= 5L))
})

test_that("missing values stop with count and positions", {
  dat <- sim_panel(3L, 15L, 5)
  dat$y[c(4, 20)] <- NA
  expect_error(
    xtpcaus(dat, "y", "x", "id", "time", test = "pfty", nboot = 99L),
    "Variable 'y' has 2 missing value\\(s\\) at row\\(s\\) 4, 20")
  dat <- sim_panel(3L, 15L, 5)
  dat$x[7] <- NA
  expect_error(
    xtpcaus(dat, "y", "x", "id", "time", test = "pfty", nboot = 99L),
    "Variable 'x' has 1 missing value\\(s\\) at row\\(s\\) 7")
})

test_that("duplicated or misaligned units are refused", {
  dat <- sim_panel(5L, 15L, 9)
  # unit 5 loses its last 3 rows, unit 1 gets 3 duplicated rows: N*T rows
  bad <- rbind(dat[dat$id != 5 | dat$time <= 12, ],
               dat[dat$id == 1 & dat$time <= 3, ])
  expect_error(
    xtpcaus(bad, "y", "x", "id", "time", test = "pfty", nboot = 99L),
    "duplicated time value")
  short <- dat[!(dat$id == 5 & dat$time > 12), ]
  expect_error(
    xtpcaus(short, "y", "x", "id", "time", test = "pfty", nboot = 99L),
    "strongly balanced")
})

test_that("singular unit: warning, NA Wald, excluded from panel statistics", {
  N <- 4L; TT <- 25L
  dat <- sim_panel(N, TT, 31)
  dat$x[dat$id == 2] <- 1          # constant x: collinear with the intercept
  expect_warning(
    res <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty",
                   pmax = 1L, dmax = 0L, kmax = 0L, nboot = 99L, seed = 2L),
    "Unit\\(s\\) 2: singular regressor matrix")
  expect_true(is.na(res$ind_wald[2]) && is.na(res$ind_pval_b[2]))
  expect_equal(res$N_used, 3L)
  expect_equal(res$fisher_df, 6L)
  expect_equal(res$fisher, -2 * sum(log(res$ind_pval_b[-2])), tolerance = 1e-10)
  expect_equal(res$wbar, mean(res$ind_wald[-2]), tolerance = 1e-10)
  expect_invisible(print(res))
})

test_that("seed argument leaves the global RNG state untouched", {
  dat <- sim_panel(3L, 15L, 8)
  set.seed(123); a <- runif(1)
  set.seed(123)
  r1 <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty", pmax = 1L,
                nboot = 99L, seed = 7L)
  b <- runif(1)
  expect_identical(a, b)
  # reproducible for the same seed
  r2 <- xtpcaus(dat, "y", "x", "id", "time", test = "pfty", pmax = 1L,
                nboot = 99L, seed = 7L)
  expect_identical(r1$ind_pval_b, r2$ind_pval_b)
})

test_that("no arrows or ampersands between author names in printed output", {
  dat <- grunfeld_panel()
  res <- xtpcaus(dat, y = "invest", x = "mvalue", panel_id = "firm",
                 time_id = "year", test = "pfty", pmax = 1L, dmax = 1L,
                 nboot = 99L, kmax = 1L, seed = 1L)
  out <- capture.output(print(res), type = "message")
  expect_false(any(grepl("=>|->", out)))
  expect_true(any(grepl("not reported", out)))
})
