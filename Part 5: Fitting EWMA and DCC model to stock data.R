## ============================================
## EWMA vs DCC without 'ccgarch' (robust script)
## - Uses rmgarch if available; otherwise manual DCC(1,1)
## ============================================

suppressPackageStartupMessages({
  library(tseries)
  library(zoo)
  # fGarch is only needed for fallback; rmgarch path uses rugarch/rmgarch
  # install.packages("fGarch")   # <- uncomment if you need fallback
})

## -----------------------------
## 1) Data: MSFT & IBM (AdjClose)
## -----------------------------
p1 <- get.hist.quote("msft", start="2000-01-01", end="2009-12-31",
                     quote="AdjClose", quiet=TRUE)
p2 <- get.hist.quote("ibm",  start="2000-01-01", end="2009-12-31",
                     quote="AdjClose", quiet=TRUE)

p <- na.omit(merge(p1, p2))               # zoo 2-column
colnames(p) <- c("MSFT","IBM")
y <- diff(log(p)) * 100                   # returns in %
y <- na.omit(y)

## -----------------------------
## 2) EWMA vols/corr (lambda=0.94)
## -----------------------------
lambda <- 0.94
Tn <- nrow(y)

# EWMA covariance recursion
EWMA <- matrix(NA_real_, nrow=Tn, ncol=3)   # [var1, var2, cov12]
S <- cov(coredata(y))                        # initialize at sample cov
EWMA[1,] <- c(S[1,1], S[2,2], S[1,2])

for (t in 2:Tn) {
  yt <- matrix(as.numeric(y[t,]), ncol=1)
  S  <- lambda * S + (1 - lambda) * (yt %*% t(yt))
  EWMA[t,] <- c(S[1,1], S[2,2], S[1,2])
}

EWMArho <- EWMA[,3] / sqrt(EWMA[,1] * EWMA[,2])

ewma_var <- zoo(cbind(MSFT = EWMA[,1], IBM = EWMA[,2]), order.by = index(y))
rho_ewma <- zoo(EWMArho, order.by = index(y))

## -----------------------------
## 3) Try rmgarch; else fallback
## -----------------------------
have_rugarch <- requireNamespace("rugarch", quietly = TRUE)
have_rmgarch <- requireNamespace("rmgarch", quietly = TRUE)

if (have_rugarch && have_rmgarch) {
  # ---------- CORRECT rmgarch DCC PATH ----------
  library(rugarch); library(rmgarch)
  
  # Univariate sGARCH(1,1) spec (zero mean)
  ug <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
    distribution.model = "norm"
  )
  
  # NOTE: multispec is from rugarch, NOT rmgarch
  uspec <- rugarch::multispec( replicate(2, ug) )
  
  # DCC(1,1) spec + fit
  dccspec <- rmgarch::dccspec(uspec = uspec,
                              dccOrder = c(1,1),
                              distribution = "mvnorm")
  dccfit  <- rmgarch::dccfit(dccspec, data = as.matrix(zoo::coredata(y)))
  
  # Extract dynamic correlation ρ12,t
  Rt    <- rmgarch::rcor(dccfit)         # array 2 x 2 x T
  rho12 <- Rt[1, 2, ]
  rho_dcc <- zoo(as.numeric(rho12), order.by = index(y))
  
  dcc_label <- "DCC (rmgarch)"
  
} else {
  # ---------- MANUAL DCC FALLBACK (no rmgarch) ----------
  stopifnot(requireNamespace("fGarch", quietly = TRUE))
  library(fGarch)
  
  y1 <- as.numeric(y[,"MSFT"]); y2 <- as.numeric(y[,"IBM"])
  
  fit_g11 <- function(x) {
    g <- fGarch::garchFit(~ garch(1,1), data = x, include.mean = FALSE,
                          cond.dist = "norm", trace = FALSE)
    list(sigma = volatility(g), z = residuals(g, standardize = TRUE))
  }
  
  g1 <- fit_g11(y1)
  g2 <- fit_g11(y2)
  
  Z <- cbind(g1$z, g2$z)
  Z <- Z[complete.cases(Z), , drop = FALSE]
  idx <- tail(index(y), nrow(Z))
  
  Sbar <- cov(Z)
  
  pack   <- function(a,b) c(qlogis(a/0.999), qlogis(b/(0.999 - a)))
  unpack <- function(th) {
    a <- 0.999 * plogis(th[1])
    b <- (0.999 - a) * plogis(th[2])
    c(a=a, b=b)
  }
  
  dcc_nll <- function(th) {
    ab <- unpack(th); a <- ab["a"]; b <- ab["b"]
    if (a <= 0 || b <= 0 || a + b >= 0.999) return(1e12)
    Q <- Sbar; ll <- 0; Tz <- nrow(Z)
    for (t in 2:Tz) {
      e <- matrix(Z[t-1,], ncol=1)
      Q <- (1 - a - b) * Sbar + a * (e %*% t(e)) + b * Q
      d <- sqrt(diag(Q)); R <- diag(1/d) %*% Q %*% diag(1/d)
      x <- matrix(Z[t,], ncol=1)
      detR <- determinant(R, logarithm = TRUE)$modulus[1]
      ll <- ll + 0.5 * (detR + t(x) %*% solve(R) %*% x)
    }
    as.numeric(ll)
  }
  
  th0   <- pack(0.05, 0.94)
  opt   <- optim(th0, dcc_nll, method = "BFGS", control = list(maxit = 2000))
  abhat <- unpack(opt$par); a_hat <- abhat["a"]; b_hat <- abhat["b"]
  
  Q <- Sbar
  rho <- rep(NA_real_, nrow(Z))
  rho[1] <- Sbar[1,2] / sqrt(Sbar[1,1] * Sbar[2,2])
  for (t in 2:nrow(Z)) {
    e <- matrix(Z[t-1,], ncol=1)
    Q <- (1 - a_hat - b_hat) * Sbar + a_hat * (e %*% t(e)) + b_hat * Q
    d <- sqrt(diag(Q)); R <- diag(1/d) %*% Q %*% diag(1/d)
    rho[t] <- R[1,2]
  }
  rho_dcc <- zoo(as.numeric(rho), order.by = idx)
  dcc_label <- sprintf("DCC (manual) α=%.3f, β=%.3f", a_hat, b_hat)
}


## -----------------------------
## 4) Plots
## -----------------------------
op <- par(mfrow = c(3,1), mar = c(4,4,3,1))

# (A) EWMA conditional variances
plot.zoo(ewma_var, plot.type = "single",
         col = c("steelblue","firebrick"), lwd = 1.6,
         xlab = "", ylab = "Variance (%^2)",
         main = sprintf("EWMA Conditional Variances (λ = %.2f)", lambda))
grid(col="grey80", lty=1, lwd=0.8)
legend("topleft", bty="n", lwd=1.6, col=c("steelblue","firebrick"),
       legend=c("MSFT var","IBM var"))

# (B) EWMA correlation only
plot.zoo(rho_ewma, col = "darkgreen", lwd = 1.6,
         xlab = "", ylab = "Correlation", ylim = c(-1,1),
         main = "EWMA Conditional Correlation — MSFT vs IBM")
abline(h=0, col="grey50"); grid(col="grey80", lty=1, lwd=0.8)

# (C) EWMA vs DCC overlay (aligned indices)
rho_e <- rho_ewma[index(rho_ewma) %in% index(rho_dcc)]
rho_d <- rho_dcc[index(rho_dcc) %in% index(rho_ewma)]
# re-align in case of small length mismatch
common_idx <- intersect(index(rho_e), index(rho_d))
rho_e <- rho_e[common_idx]; rho_d <- rho_d[common_idx]

plot(index(rho_e), as.numeric(rho_e), type = "l",
     col = "darkgreen", lwd = 1.6, xlab = "", ylab = "Correlation",
     ylim = c(-1,1),
     main = "Conditional Correlation: EWMA vs DCC")
lines(index(rho_d), as.numeric(rho_d), col = "royalblue", lwd = 1.4)
abline(h=0, col="grey50"); grid(col="grey80", lty=1, lwd=0.8)
legend("topleft", bty="n", lwd=c(1.6,1.4),
       col=c("darkgreen","royalblue"),
       legend=c("EWMA","DCC"))

par(op)

## ============================================================
## 5) EWMA vs DCC — hard comparison (MSE, R2, ACF/PACF/QQ, whiteness, LL)
## ============================================================

stopifnot(exists("y"), exists("rho_ewma"), exists("ewma_var"), exists("rho_dcc"))

library(zoo)
library(stats)

## --- helper: rolling realized correlation (proxy target) ---
roll_realized_corr <- function(y, window = 60L) {
  stopifnot(ncol(y) == 2)
  Y <- zoo::coredata(y)
  n <- nrow(Y)
  out <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    if (t >= window) {
      out[t] <- cor(Y[(t - window + 1):t, 1],
                    Y[(t - window + 1):t, 2])
    }
  }
  zoo(out, order.by = index(y))
}

## --- helper: build H_t (2x2) from sigmas and rho ---
Ht_from_sigmarho <- function(sig1, sig2, rho) {
  D <- diag(c(sig1, sig2))
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  D %*% R %*% D
}

## --- helper: whiten residual y_t with D_t and rho_t ---
whiten_one <- function(y_t, sig1, sig2, rho) {
  s1 <- max(sig1, 1e-8); s2 <- max(sig2, 1e-8)
  D  <- diag(c(s1, s2))
  u  <- solve(D, y_t)                          # remove vol
  R  <- matrix(c(1, rho, rho, 1), 2, 2)
  U  <- chol(R)                                 # R = U'U (upper tri)
  z  <- solve(U, u)                             # decorrelate: Var(z)=I
  as.numeric(z)
}

## --- helper: per-time bivariate normal log-likelihood ---
ll_one <- function(y_t, H) {
  # constant term cancels in comparisons, but we include it
  d <- 2L
  log2pi <- log(2*pi)
  detH   <- determinant(H, logarithm = TRUE)$modulus[1]
  quad   <- drop(t(y_t) %*% solve(H) %*% y_t)
  -0.5 * (d*log2pi + detH + quad)
}

## ============= Prepare inputs =============
# Align model series on a common index
idx_common <- intersect(index(rho_ewma), index(rho_dcc))
rho_e <- rho_ewma[idx_common]
rho_d <- rho_dcc[idx_common]
yy    <- y[idx_common, ]

# EWMA sigmas (sqrt of EWMA variances)
sig_e1 <- sqrt(coredata(ewma_var[idx_common, 1]))
sig_e2 <- sqrt(coredata(ewma_var[idx_common, 2]))

# DCC sigmas:
if (exists("dccfit")) {
  # rmgarch path
  library(rmgarch)
  Rc  <- rmgarch::rcov(dccfit)                 # 2x2xT array
  # extract sigmas from variances
  sig_d1 <- sqrt(Rc[1,1,])
  sig_d2 <- sqrt(Rc[2,2,])
  # ensure same length as rho_d
  len <- min(length(sig_d1), length(rho_d))
  sig_d1 <- sig_d1[(length(sig_d1)-len+1):length(sig_d1)]
  sig_d2 <- sig_d2[(length(sig_d2)-len+1):length(sig_d2)]
  rho_d  <- rho_d[(length(rho_d)-len+1):length(rho_d)]
  yy     <- yy[(nrow(yy)-len+1):nrow(yy), ]
  rho_e  <- rho_e[index(yy)]
  sig_e1 <- sig_e1[index(yy)]
  sig_e2 <- sig_e2[index(yy)]
  dcc_path_label <- "DCC (rmgarch)"
} else if (exists("g1") && exists("g2")) {
  # manual path (univariate GARCH -> DCC recursion)
  sig_d1 <- as.numeric(g1$sigma)
  sig_d2 <- as.numeric(g2$sigma)
  # align lengths
  L <- min(length(sig_d1), length(sig_d2), length(rho_d), nrow(yy))
  sig_d1 <- tail(sig_d1, L); sig_d2 <- tail(sig_d2, L)
  rho_d  <- tail(as.numeric(rho_d), L)
  yy     <- tail(yy, L)
  rho_e  <- tail(as.numeric(rho_e), L)
  sig_e1 <- tail(sig_e1, L); sig_e2 <- tail(sig_e2, L)
  dcc_path_label <- "DCC (manual)"
} else {
  stop("Cannot locate DCC sigmas — run DCC step first (rmgarch or manual).")
}

## ============= (1) Accuracy vs realized correlation =============
W <- 60L                            # rolling window length (days)
rho_real <- roll_realized_corr(yy, window = W)
# align all to common valid (non-NA) window
common <- Reduce(intersect, list(index(rho_real), index(yy)))
rho_real <- rho_real[common]
yy       <- yy[common, ]
rho_e    <- zoo(as.numeric(rho_e[common]), order.by = common)
rho_d    <- zoo(as.numeric(rho_d[common]), order.by = common)
sig_e1   <- as.numeric(sig_e1[common]); sig_e2 <- as.numeric(sig_e2[common])
# sig_d1/2 might be numeric vectors; trim to match
L <- length(common)
sig_d1 <- tail(sig_d1, L); sig_d2 <- tail(sig_d2, L)

mse <- function(pred, target) mean((pred - target)^2, na.rm = TRUE)
r2  <- function(pred, target) {
  t <- as.numeric(target); p <- as.numeric(pred)
  ok <- is.finite(t) & is.finite(p)
  t <- t[ok]; p <- p[ok]
  1 - sum((t - p)^2) / sum((t - mean(t))^2)
}

mse_EWMA <- mse(rho_e, rho_real)
mse_DCC  <- mse(rho_d, rho_real)
r2_EWMA  <- r2(rho_e, rho_real)
r2_DCC   <- r2(rho_d, rho_real)

cat(sprintf("\n[Correlation vs realized (W=%d)]  MSE  — EWMA: %.6f | DCC: %.6f\n", W, mse_EWMA, mse_DCC))
cat(sprintf("[Correlation vs realized (W=%d)]  R^2  — EWMA: %.4f  | DCC: %.4f\n", W, r2_EWMA, r2_DCC))

## ---------- BULLETPROOF ALIGNMENT (EWMA vs DCC inputs) ----------
library(zoo)

# Ensure model pieces are zoo with dates
rho_e_z <- zoo(as.numeric(rho_ewma), order.by = index(rho_ewma))   # EWMA rho
rho_d_z <- zoo(as.numeric(rho_dcc),  order.by = index(rho_dcc))    # DCC  rho
sig_e1_z <- sqrt(ewma_var[, 1])   # EWMA sigma1 (already zoo via ewma_var)
sig_e2_z <- sqrt(ewma_var[, 2])   # EWMA sigma2

if (exists("dccfit")) {
  # rmgarch path: rcov aligned to the input data (y)
  library(rmgarch)
  Rc <- rmgarch::rcov(dccfit)                 # 2 x 2 x T array
  sig_d1_z <- zoo(sqrt(Rc[1,1,]), order.by = index(y))
  sig_d2_z <- zoo(sqrt(Rc[2,2,]), order.by = index(y))
  dcc_path_label <- "DCC (rmgarch)"
} else if (exists("g1") && exists("g2")) {
  # manual path: g1/g2 from fallback
  sig_d1_z <- zoo(as.numeric(g1$sigma), order.by = index(y))
  sig_d2_z <- zoo(as.numeric(g2$sigma), order.by = index(y))
  dcc_path_label <- "DCC (manual)"
} else {
  stop("DCC components not found. Run the DCC step first.")
}

# Merge everything on common dates; drop any NA row
M <- na.omit(merge(
  y[,1], y[,2],       # returns
  rho_e_z, rho_d_z,   # correlations
  sig_e1_z, sig_e2_z, # EWMA sigmas
  sig_d1_z, sig_d2_z, # DCC  sigmas
  all = FALSE
))
colnames(M) <- c("y1","y2","rho_e","rho_d","sig_e1","sig_e2","sig_d1","sig_d2")

if (nrow(M) <= 10L) stop("Too few aligned observations after merge. Check data span/availability.")

# Extract aligned objects for downstream code
yy     <- M[, c("y1","y2")]
rho_e  <- M[, "rho_e"]
rho_d  <- M[, "rho_d"]
sig_e1 <- as.numeric(M[, "sig_e1"])
sig_e2 <- as.numeric(M[, "sig_e2"])
sig_d1 <- as.numeric(M[, "sig_d1"])
sig_d2 <- as.numeric(M[, "sig_d2"])


## --- PATCH 1: safer whitener (clamp rho) ---
whiten_one <- function(y_t, sig1, sig2, rho) {
  s1 <- max(sig1, 1e-8); s2 <- max(sig2, 1e-8)
  rho <- max(min(rho, 0.999), -0.999)          # clamp to avoid singular R
  D  <- diag(c(s1, s2))
  u  <- solve(D, y_t)                           # de-vol
  R  <- matrix(c(1, rho, rho, 1), 2, 2)
  U  <- chol(R)                                 # R = U'U
  z  <- solve(U, u)                             # decorrelate: Var(z)=I
  as.numeric(z)
}

## --- PATCH 2: build a single NA-free mask for BOTH models ---
# Convert zoo to numeric vectors for masking
rho_e_num <- as.numeric(rho_e)
rho_d_num <- as.numeric(rho_d)

ok_all <- complete.cases(yy) &
  is.finite(sig_e1) & is.finite(sig_e2) & is.finite(rho_e_num) &
  is.finite(sig_d1) & is.finite(sig_d2) & is.finite(rho_d_num)

# Apply mask consistently
yy       <- yy[ok_all, ]
rho_e    <- zoo(rho_e_num[ok_all], order.by = index(yy))
rho_d    <- zoo(rho_d_num[ok_all], order.by = index(yy))
sig_e1   <- sig_e1[ok_all]; sig_e2 <- sig_e2[ok_all]
sig_d1   <- sig_d1[ok_all]; sig_d2 <- sig_d2[ok_all]

stopifnot(nrow(yy) > 10)  # sanity

## ============= (2) Whitened residuals + (3) tests =============
# Build whitened residuals z_t for both models
Z_E <- matrix(NA_real_, nrow = nrow(yy), ncol = 2)
Z_D <- matrix(NA_real_, nrow = nrow(yy), ncol = 2)
Ymat <- as.matrix(yy)

for (t in seq_len(nrow(yy))) {
  # EWMA
  Z_E[t, ] <- whiten_one(Ymat[t, ], sig_e1[t], sig_e2[t], as.numeric(rho_e[t]))
  # DCC
  Z_D[t, ] <- whiten_one(Ymat[t, ], sig_d1[t], sig_d2[t], as.numeric(rho_d[t]))
}
colnames(Z_E) <- colnames(Z_D) <- c("MSFT","IBM")

# Plots: ACF, PACF, QQ for MSFT (left) and IBM (right), comparing EWMA vs DCC
op <- par(mfrow = c(3, 2), mar = c(4,4,3,1))
# (A) ACF of z_t (EWMA vs DCC) — MSFT
acf(Z_E[,1], 40, main = "EWMA: ACF(z) MSFT"); grid(col="grey80", lty=1, lwd=0.8)
acf(Z_D[,1], 40, main = paste(dcc_path_label, ": ACF(z) MSFT")); grid(col="grey80", lty=1, lwd=0.8)
# (B) PACF of z_t
pacf(Z_E[,1], 40, main = "EWMA: PACF(z) MSFT"); grid(col="grey80", lty=1, lwd=0.8)
pacf(Z_D[,1], 40, main = paste(dcc_path_label, ": PACF(z) MSFT")); grid(col="grey80", lty=1, lwd=0.8)
# (C) QQ plot
qqnorm(Z_E[,1], main = "EWMA: QQ(z) MSFT"); qqline(Z_E[,1], col="red", lwd=1.2)
qqnorm(Z_D[,1], main = paste(dcc_path_label, ": QQ(z) MSFT")); qqline(Z_D[,1], col="red", lwd=1.2)
par(op)

# Same for IBM if you want:
op <- par(mfrow = c(3, 2), mar = c(4,4,3,1))
acf(Z_E[,2], 40, main = "EWMA: ACF(z) IBM");  grid(col="grey80", lty=1, lwd=0.8)
acf(Z_D[,2], 40, main = paste(dcc_path_label, ": ACF(z) IBM")); grid(col="grey80", lty=1, lwd=0.8)
pacf(Z_E[,2], 40, main = "EWMA: PACF(z) IBM"); grid(col="grey80", lty=1, lwd=0.8)
pacf(Z_D[,2], 40, main = paste(dcc_path_label, ": PACF(z) IBM")); grid(col="grey80", lty=1, lwd=0.8)
qqnorm(Z_E[,2], main = "EWMA: QQ(z) IBM"); qqline(Z_E[,2], col="red", lwd=1.2)
qqnorm(Z_D[,2], main = paste(dcc_path_label, ": QQ(z) IBM")); qqline(Z_D[,2], col="red", lwd=1.2)
par(op)

## Whiteness tests (Ljung–Box) on z and z^2 (lags 10 and 20)
lb_report <- function(z, name) {
  res <- list(
    LB10  = Box.test(z,    lag = 10, type = "Ljung-Box"),
    LB20  = Box.test(z,    lag = 20, type = "Ljung-Box"),
    LB10s = Box.test(z^2,  lag = 10, type = "Ljung-Box"),
    LB20s = Box.test(z^2,  lag = 20, type = "Ljung-Box")
  )
  cat(sprintf("\n[%s] Ljung–Box p-values:\n", name))
  cat(sprintf("  LB(10):  z  p=%.3f | z^2 p=%.3f\n",
              res$LB10$p.value, res$LB10s$p.value))
  cat(sprintf("  LB(20):  z  p=%.3f | z^2 p=%.3f\n",
              res$LB20$p.value, res$LB20s$p.value))
  invisible(res)
}

cat("\n--- Whiteness (MSFT) ---\n")
lb_E_MSFT <- lb_report(Z_E[,1], "EWMA MSFT")
lb_D_MSFT <- lb_report(Z_D[,1], paste(dcc_path_label, "MSFT"))

cat("\n--- Whiteness (IBM) ---\n")
lb_E_IBM  <- lb_report(Z_E[,2], "EWMA IBM")
lb_D_IBM  <- lb_report(Z_D[,2], paste(dcc_path_label, "IBM"))

## ============= (4) Log-likelihood comparison =============
LL_model <- function(sig1, sig2, rho, yy) {
  stopifnot(length(sig1)==nrow(yy), length(sig2)==nrow(yy), length(rho)==nrow(yy))
  ll <- 0
  for (t in seq_len(nrow(yy))) {
    Ht <- Ht_from_sigmarho(sig1[t], sig2[t], rho[t])
    ll <- ll + ll_one(as.numeric(yy[t, ]), Ht)
  }
  ll
}

LL_EWMA <- LL_model(sig_e1, sig_e2, as.numeric(rho_e), as.matrix(yy))
LL_DCC  <- LL_model(sig_d1, sig_d2, as.numeric(rho_d), as.matrix(yy))

cat(sprintf("\n[Gaussian log-likelihood over %d days]\n  EWMA: %.2f\n  %s: %.2f\n",
            nrow(yy), LL_EWMA, dcc_path_label, LL_DCC))

## ===========================
## 5) Information Criteria: AIC, BIC, KLIC
## ===========================
# Uses LL_EWMA, LL_DCC, yy, and (if available) dccfit created above.

n_obs <- nrow(yy)

# Parameter counts:
# - EWMA: set k_EWMA = 1 if you regard λ as a fitted parameter;
#         set k_EWMA = 0 if λ was fixed a priori.
k_EWMA <- 1

# - DCC: if rmgarch path, count directly; else manual fallback = 2×GARCH(1,1) + DCC(α,β) = 8
if (exists("dccfit")) {
  k_DCC <- length(coef(dccfit))
  dcc_ic_label <- "DCC (rmgarch)"
} else if (exists("g1") && exists("g2")) {
  k_DCC <- 8L
  dcc_ic_label <- dcc_label
} else {
  stop("Cannot determine parameter count for DCC — run the DCC step first.")
}

# AIC / BIC
AIC_EWMA <- -2 * as.numeric(LL_EWMA) + 2 * k_EWMA
BIC_EWMA <- -2 * as.numeric(LL_EWMA) + k_EWMA * log(n_obs)

AIC_DCC  <- -2 * as.numeric(LL_DCC)  + 2 * k_DCC
BIC_DCC  <- -2 * as.numeric(LL_DCC)  + k_DCC  * log(n_obs)

# Kullback–Leibler IC (per-observation deviance)
KLIC_EWMA <- (-2 * as.numeric(LL_EWMA)) / n_obs
KLIC_DCC  <- (-2 * as.numeric(LL_DCC))  / n_obs

cat("\n--- Information Criteria ---\n")
cat(sprintf("EWMA  : AIC = %.2f | BIC = %.2f | KLIC (per obs) = %.6f\n",
            AIC_EWMA, BIC_EWMA, KLIC_EWMA))
cat(sprintf("%s: AIC = %.2f | BIC = %.2f | KLIC (per obs) = %.6f\n",
            dcc_ic_label, AIC_DCC, BIC_DCC, KLIC_DCC))


## ============= (Nice-to-have) Overlay realized vs model =============
op <- par(mfrow = c(1,1), mar = c(4,4,3,1))
plot(rho_real, col="black", lwd=1.2, ylim=c(-1,1),
     main = sprintf("Realized corr (W=%d) vs EWMA vs %s", W, dcc_path_label),
     xlab = "", ylab = "Correlation")
lines(rho_e, col="darkgreen", lwd=1.4)
lines(rho_d, col="royalblue", lwd=1.4)
abline(h=0, col="grey50"); grid(col="grey80", lty=1, lwd=0.8)
legend("topleft", bty="n", lwd=c(1.2,1.4,1.4),
       col=c("black","darkgreen","royalblue"),
       legend=c("Realized (roll W)","EWMA","DCC"))
par(op)
