# =========================================================
# Portfolio VaR (3 stocks: MS, IBM, AAPL) with EWMA-t, MA-t,
# HS, GARCH-t (on portfolio) + Copula-FHS (per-asset GARCH-t + t-copula)
# =========================================================

library(tseries)
library(zoo)
library(fGarch)
library(copula)   # for the t-copula

# ---------- Data ----------
tickers <- c("MS","IBM","AAPL")
get_px <- function(sym) get.hist.quote(instrument=sym,
                                       start="1994-02-11", end="2009-12-31",
                                       quote="AdjClose", quiet=TRUE)
px_list <- lapply(tickers, get_px)

# align to common dates
px <- do.call(merge, lapply(px_list, function(z) as.zoo(z)))
colnames(px) <- tickers
px <- na.locf(px, na.rm = FALSE)
px <- na.omit(px)

# log-returns per asset
R <- diff(log(px))
R <- na.omit(R)                      # zoo matrix (T x 3)
Tn <- nrow(R)

# portfolio weights (daily rebalanced)
w <- c(1/3, 1/3, 1/3); w <- w / sum(w)

# portfolio return series
y <- as.numeric(R %*% w)             # vector length Tn

# ---------- Prelim ----------
WE <- 1000
p  <- 0.01
l1 <- ceiling((WE + 1) * p)          # HS index (type-1)
value <- 1

USE_COPULA <- TRUE                   # toggle Copula-FHS

# Allocate VaR matrix (add CopulaFHS column if enabled)
k_methods <- if (USE_COPULA) 5 else 4
VaR <- matrix(NA_real_, nrow = Tn, ncol = k_methods)
colnames(VaR) <- c("EWMA","MA","HS","GARCH", if (USE_COPULA) "CopulaFHS" else NULL)

# EWMA warm-up on portfolio
lambda <- 0.94
s11 <- var(y[1:30])
for (t in 2:WE) s11 <- lambda*s11 + (1 - lambda)*y[t-1]^2

# helpers: excess kurtosis -> df for t
excess_kurt <- function(x){
  m <- mean(x); s <- sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean(((x - m)/s)^4) - 3
}
df_from_kurt <- function(k){
  if (!is.finite(k) || k <= 0) return(1000)  # ~Normal if no/neg excess kurtosis
  min(1000, 4 + 6 / k)                       # t-kurtosis: 6/(nu-4)
}

# ---------- Rolling VaR on portfolio ----------
for (t in (WE + 1):Tn) {
  t1 <- t - WE; t2 <- t - 1
  window_p <- y[t1:t2]                       # portfolio window
  window_R <- coredata(R[t1:t2, , drop=FALSE])  # WE x 3 matrix (assets)
  
  # --- EWMA sigma (on portfolio) ---
  s11 <- lambda*s11 + (1 - lambda)*y[t - 1]^2
  sigma_ewma <- sqrt(s11)
  
  # heavier-tailed quantile from portfolio window
  k_win <- excess_kurt(window_p)
  nu_hat <- df_from_kurt(k_win)
  q_t <- qt(p, df = nu_hat)      # negative
  
  # EWMA-t & MA-t VaR (on portfolio)
  VaR[t, "EWMA"] <- -sigma_ewma * q_t * value
  VaR[t, "MA"]   <- -sd(window_p) * q_t * value
  
  # HS VaR on portfolio
  ys <- sort(window_p)
  VaR[t, "HS"] <- -ys[l1] * value
  
  # GARCH(1,1)-t on portfolio; use model's 1-step sigma & t-quantile
  g_p <- garchFit(~ garch(1,1), data = window_p, include.mean = FALSE,
                  cond.dist = "std", trace = FALSE)
  pred_p <- predict(g_p, n.ahead = 1)
  if ("standardDeviation" %in% names(pred_p)) {
    sig1_p <- as.numeric(pred_p$standardDeviation[1])
  } else if ("sigma.t" %in% names(pred_p)) {
    sig1_p <- as.numeric(pred_p$sigma.t[1])
  } else {
    pars <- coef(g_p); h_t <- tail(g_p@h.t, 1); eps2 <- tail(window_p, 1)^2
    h1 <- pars["omega"] + pars["alpha1"] * eps2 + pars["beta1"] * h_t
    sig1_p <- sqrt(h1)
  }
  nu_g_p <- as.numeric(coef(g_p)["shape"]); if (!is.finite(nu_g_p)) nu_g_p <- 8
  VaR[t, "GARCH"] <- -qstd(p, mean=0, sd=sig1_p, nu=nu_g_p) * value
  
  # --- Copula-FHS (per-asset GARCH-t + t-copula simulation) ---
  if (USE_COPULA) {
    # Fit GARCH(1,1)-t per asset on window
    gfits <- vector("list", 3)
    for (j in 1:3) {
      gfits[[j]] <- garchFit(~ garch(1,1), data = window_R[, j],
                             include.mean = FALSE, cond.dist = "std", trace = FALSE)
    }
    
    # Standardized residuals Z_j = res_j / sigma_j
    Z <- sapply(1:3, function(j) {
      gj  <- gfits[[j]]
      res <- residuals(gj)
      sig <- sqrt(gj@h.t)
      as.numeric(res / sig)
    })
    
    # Uniforms via marginal t CDFs (df = shape_j)
    nus <- sapply(gfits, function(gj) as.numeric(coef(gj)["shape"]))
    U   <- sapply(1:3, function(j) pt(Z[, j], df = nus[j]))
    U   <- pmin(pmax(U, 1e-6), 1 - 1e-6)   # clamp
    
    # Fit t-copula on U
    tc_init <- tCopula(dim = 3, dispstr = "un", df = 8, df.fixed = FALSE)
    fit <- fitCopula(tc_init, data = U, method = "ml")
    tc  <- fit@copula
    
    # Next-day conditional sigmas per asset
    sig1 <- sapply(1:3, function(j) {
      gj <- gfits[[j]]
      pr <- predict(gj, n.ahead = 1)
      if ("standardDeviation" %in% names(pr)) as.numeric(pr$standardDeviation[1])
      else if ("sigma.t" %in% names(pr))      as.numeric(pr$sigma.t[1])
      else {
        pars <- coef(gj); h_t <- tail(gj@h.t, 1); eps2 <- tail(window_R[, j], 1)^2
        sqrt(pars["omega"] + pars["alpha1"]*eps2 + pars["beta1"]*h_t)
      }
    })
    
    # Simulate joint shocks and map to portfolio
    Nsim <- 20000
    Usim <- rCopula(Nsim, tc)                 # N x 3 uniforms
    Qsim <- sapply(1:3, function(j) qt(Usim[, j], df = nus[j]))  # t-shocks
    Rsim <- sweep(Qsim, 2, sig1, `*`)         # scale by next-day sigmas
    Psim <- as.numeric(Rsim %*% w)            # simulated portfolio shocks
    
    VaR[t, "CopulaFHS"] <- -quantile(Psim, probs = p, type = 1)
  }
}

W1 <- WE + 1

# ---------- Backtest prints ----------
for (i in 1:ncol(VaR)) {
  VR <- sum(y[W1:Tn] < -VaR[W1:Tn, i], na.rm = TRUE) / (p * (Tn - WE))
  s  <- sd(VaR[W1:Tn, i], na.rm = TRUE)
  cat(i, colnames(VaR)[i], "VR", round(VR,3), "VaR vol", round(s,6), "\n")
}

# ---------- Plot (auto-legend incl. Copula if present) ----------
cols <- c("black","turquoise","green","blue","orange","red3")
matplot(cbind(y[W1:Tn], -VaR[W1:Tn, ]), type='l', lty=1,
        col=cols[seq_len(ncol(VaR)+1)],
        lwd=c(1.6, rep(1.2, ncol(VaR))), xlab="Zeit", ylab="Portf.-Return / -VaR (1%)")
legend("bottomleft",
       legend=c("Portfolio Returns", colnames(VaR)),
       col=cols[seq_len(ncol(VaR)+1)], lty=1,
       lwd=c(1.6, rep(1.2, ncol(VaR))), bty="n")
abline(h=0, col="grey80")

# ---------- Kupiec & Christoffersen on portfolio ----------
bern_test <- function (p, v) {
  a <- p^(sum(v)) * (1-p)^(length(v)-sum(v))
  b <- (sum(v)/length(v))^(sum(v)) * (1-(sum(v)/length(v)))^(length(v)-sum(v))
  -2*log(a/b)
}
ind_test <- function (V) {
  J <- matrix(ncol=4, nrow=length(V))
  for (i in 2:length(V)) {
    J[i,1] <- (V[i-1]==0 & V[i]==0)
    J[i,2] <- (V[i-1]==0 & V[i]==1)
    J[i,3] <- (V[i-1]==1 & V[i]==0)
    J[i,4] <- (V[i-1]==1 & V[i]==1)
  }
  V_00 <- sum(J[,1],na.rm=TRUE); V_01 <- sum(J[,2],na.rm=TRUE)
  V_10 <- sum(J[,3],na.rm=TRUE); V_11 <- sum(J[,4],na.rm=TRUE)
  p_00 <- V_00/(V_00+V_01); p_01 <- V_01/(V_00+V_01)
  p_10 <- V_10/(V_10+V_11); p_11 <- V_11/(V_10+V_11)
  hat_p <- (V_01+V_11)/(V_00+V_01+V_10+V_11)
  a <- (1-hat_p)^(V_00+V_10)*(hat_p)^(V_01+V_11)
  b <- (p_00)^(V_00)*(p_01)^(V_01)*(p_10)^(V_10)*p_11^(V_11)
  -2*log(a/b)
}

ya  <- y[W1:Tn]; VaRa <- VaR[W1:Tn, , drop=FALSE]
q   <- sweep(VaRa, 1, ya, function(th, ret) ret < -th)

for (i in 1:ncol(VaR)) {
  v   <- as.integer(q[, i])
  ber <- bern_test(p, v); ind <- ind_test(v)
  cat(colnames(VaR)[i], "Bernoulli", ber, 1-pchisq(ber,1),
      "independence", ind, 1-pchisq(ind,1), "\n")
}

# ---------- Exceedance counts (robust) ----------
labs <- colnames(VaR)
yy <- ya
X  <- VaRa
ok <- is.finite(yy) & rowSums(is.finite(X)) == ncol(X)
yy <- yy[ok]; X <- X[ok, , drop = FALSE]
E  <- yy < -X
exceed_counts <- colSums(E, na.rm = TRUE)
N  <- length(yy); expected <- p * N; hit_rate <- exceed_counts / N
cat("\nExceedances by model (", N, " days; expected ≈ ", round(expected, 2), "):\n", sep="")
for (j in seq_along(labs)) {
  cat(sprintf("  %-10s: %4d   (hit rate = %5.2f%%)\n",
              labs[j], exceed_counts[j], 100*hit_rate[j]))
}
