## Size check of caustests::xtpcaus(test = "pfty") under homogeneous non-causality
## (HNC), audit data-generating process (Dumitrescu and Hurlin 2012 Section 5 and
## steps A1-A2, as in the 1.1.6 audit): y_it = a_i + g_i y_{i,t-1} + s_i e_it,
## a_i ~ N(0, 1), g_i ~ U(-0.99, 0.99), s_i^2 ~ U(0.5, 1.5), x_it ~ N(0, 1)
## independent of y, 100 burn-in periods.
##
## Two specifications are run on each replication:
##   plain:   kmax = 0, dmax = 0 (eq. 1 regression, K_i = p_i): analytic
##            Z-tilde (eq. 26), Z-bar (eq. 11 form with K_i = 1), and the
##            bootstrap Z-tilde and Fisher p-values;
##   default: dmax = 1, kmax = 3 (Fourier Toda-Yamamoto): bootstrap Fisher.
##
## Usage (from any directory, with caustests installed):
##   Rscript size_rerun_xtpcaus.R [R] [output.csv]
## Defaults: R = 300 replications, output file in tempdir(). Runtime is about
## 10 to 15 minutes for the four cells on one core.

suppressPackageStartupMessages(library(caustests))

args    <- commandArgs(trailingOnly = TRUE)
R       <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
out_csv <- if (length(args) >= 2L) args[2L] else
  file.path(tempdir(), "size_rerun_xtpcaus.csv")

N     <- 10L
nboot <- 99L

gen_panel <- function(N, TT, err, burn = 100L) {
  Tn  <- TT + burn
  a   <- rnorm(N)
  g   <- runif(N, -0.99, 0.99)
  s   <- sqrt(runif(N, 0.5, 1.5))
  fac <- if (err == "factor") rnorm(Tn) else NULL
  do.call(rbind, lapply(seq_len(N), function(i) {
    e <- if (err == "factor") sqrt(0.5) * fac + sqrt(0.5) * rnorm(Tn) else rnorm(Tn)
    y <- numeric(Tn)
    y[1] <- a[i]
    for (t in 2:Tn) y[t] <- a[i] + g[i] * y[t - 1] + s[i] * e[t]
    data.frame(id = i, time = seq_len(TT),
               y = y[(burn + 1):Tn], x = rnorm(TT))
  }))
}

set.seed(20261030)
out <- list()
for (TT in c(20L, 50L)) for (err in c("normal", "factor")) {
  rej <- matrix(NA_real_, R, 6, dimnames = list(NULL,
    c("plain_ztilde", "plain_zbar", "plain_ztilde_boot", "plain_fisher",
      "default_fisher", "plain_wbar_dh")))
  for (r in seq_len(R)) {
    d  <- gen_panel(N, TT, err)
    rp <- suppressWarnings(xtpcaus(d, "y", "x", "id", "time", test = "pfty",
                                   pmax = 1L, dmax = 0L, kmax = 0L, nboot = nboot))
    rd <- suppressWarnings(xtpcaus(d, "y", "x", "id", "time", test = "pfty",
                                   pmax = 1L, dmax = 1L, kmax = 3L, nboot = nboot))
    rej[r, ] <- c(rp$ztilde_pv < 0.05, rp$zbar_pv < 0.05,
                  rp$ztilde_boot_pv < 0.05, rp$fisher_pv < 0.05,
                  rd$fisher_pv < 0.05, rp$wbar_dh)
  }
  out[[length(out) + 1L]] <- data.frame(
    N = N, TT = TT, K = 1L, err = err, R = R, t(colMeans(rej)),
    se = sqrt(0.05 * 0.95 / R))
  message(sprintf("T = %d, errors = %s done", TT, err))
}
res <- do.call(rbind, out)
write.csv(res, out_csv, row.names = FALSE)
print(res)
