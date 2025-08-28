# =========================================================
# Portfolio VaR (3 stocks: MS, IBM, AAPL) with EWMA-t, MA-t,
# HS, and GARCH-t.
# =========================================================

library(tseries)
library(zoo)
library(fGarch)

# ---------- Data ----------
tickers <- c("MS","IBM","AAPL")
get_px <- function(sym) get.hist.quote(instrument=sym,
                                       start="1994-02-11", end="2009-12-31",
                                       quote="AdjClose", quiet=TRUE)
px_list <- lapply(tickers, get_px)

# align to common dates
px <- do.call(merge, lapply(px_list, function(z) as.zoo(z)))
colnames(px) <- tickers
px <- na.locf(px, na.rm=FALSE)
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
VaR <- matrix(NA_real_, nrow=Tn, ncol=4)
colnames(VaR) <- c("EWMA","MA","HS","GARCH")

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
  window <- y[t1:t2]
  
  # EWMA sigma (on portfolio)
  s11 <- lambda*s11 + (1 - lambda)*y[t - 1]^2
  sigma_ewma <- sqrt(s11)
  
  # heavier-tailed quantile from window
  k_win <- excess_kurt(window)
  nu_hat <- df_from_kurt(k_win)
  q_t <- qt(p, df = nu_hat)      # negative
  
  # EWMA-t & MA-t VaR
  VaR[t, "EWMA"] <- -sigma_ewma * q_t * value
  VaR[t, "MA"]   <- -sd(window)  * q_t * value
  
  # HS VaR on portfolio
  ys <- sort(window)
  VaR[t, "HS"] <- -ys[l1] * value
  
  # GARCH(1,1)-t on portfolio; use model's 1-step sigma & t-quantile
  g <- garchFit(~ garch(1,1), data = window, include.mean = FALSE,
                cond.dist = "std", trace = FALSE)
  pred <- predict(g, n.ahead = 1)
  if ("standardDeviation" %in% names(pred)) {
    sig1 <- as.numeric(pred$standardDeviation[1])
  } else if ("sigma.t" %in% names(pred)) {
    sig1 <- as.numeric(pred$sigma.t[1])
  } else {
    pars <- coef(g); h_t <- tail(g@h.t, 1); eps2 <- tail(window, 1)^2
    h1 <- pars["omega"] + pars["alpha1"] * eps2 + pars["beta1"] * h_t
    sig1 <- sqrt(h1)
  }
  nu_g <- as.numeric(coef(g)["shape"]); if (!is.finite(nu_g)) nu_g <- 8
  VaR[t, "GARCH"] <- -qstd(p, mean=0, sd=sig1, nu=nu_g) * value
}

W1 <- WE + 1

# ---------- Backtest prints ----------
for (i in 1:4) {
  VR <- sum(y[W1:Tn] < -VaR[W1:Tn, i]) / (p * (Tn - WE))
  s  <- sd(VaR[W1:Tn, i])
  cat(i, colnames(VaR)[i], "VR", round(VR,3), "VaR vol", round(s,6), "\n")
}

# ---------- Plot ----------
matplot(cbind(y[W1:Tn], -VaR[W1:Tn, ]), type='l', lty=1,
        col=c("black","turquoise","green","blue","orange"),
        lwd=c(1.6,1.2,1.2,1.2,1.2), xlab="Zeit", ylab="Portf.-Return / -VaR (1%)")
legend("bottomleft",
       legend=c("Portfolio Returns","VaR EWMA (t)","VaR MA (t)","VaR HS","VaR GARCH-t"),
       col=c("black","turquoise","green","blue","orange"), lty=1,
       lwd=c(1.6,1.2,1.2,1.2,1.2), bty="n")
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

ya  <- y[W1:Tn]; VaRa <- VaR[W1:Tn,]
q   <- ya < -VaRa
for (i in 1:4) {
  v   <- as.integer(q[,i])
  ber <- bern_test(p, v); ind <- ind_test(v)
  cat(colnames(VaR)[i], "Bernoulli", ber, 1-pchisq(ber,1),
      "independence", ind, 1-pchisq(ind,1), "\n")
}

kupiec_p <- 1 - pchisq(
  -2 * ( exceed_counts * log(p / hit_rate) + (N - exceed_counts) * log((1 - p) / (1 - hit_rate)) ),
  df = 1
)
cat("\nKupiec (UC) p-values:\n")
for (j in seq_along(labs)) {
  cat(sprintf("  %-10s: %.4g\n", labs[j], kupiec_p[j]))
}

# ---------- Exceedance counts (robust) ----------
W1 <- WE + 1
labs <- colnames(VaR)
if (is.null(labs)) labs <- paste0("M", seq_len(ncol(VaR)))

yy <- y[W1:Tn]
X  <- VaR[W1:Tn, , drop = FALSE]

# keep only rows where portfolio return and all VaR columns are finite
ok <- is.finite(yy) & rowSums(is.finite(X)) == ncol(X)
yy <- yy[ok]
X  <- X[ok, , drop = FALSE]

E  <- yy < -X                          # logical matrix (breach if return < -VaR)
exceed_counts <- colSums(E, na.rm = TRUE)
N  <- length(yy)
expected <- p * N
hit_rate  <- exceed_counts / N

cat("\nExceedances by model (", N, " days; expected ≈ ", round(expected, 2), "):\n", sep = "")
for (j in seq_along(labs)) {
  cat(sprintf("  %-10s: %4d   (hit rate = %5.2f%%)\n", labs[j], exceed_counts[j], 100*hit_rate[j]))
}
