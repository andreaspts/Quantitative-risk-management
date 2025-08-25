#improvement of first file so that exceedances are more acceptable:
#based on stock index data VaR are calculated in 4 ways
#all in sample
#plots of log returns vs. respective VaRs computed
#backtesting to check if VaR exceedances are Bernoulli distributed and independent
#exceedance count to get a feeling if evaluation seems reasoable

library("tseries")
library("zoo")  
p=get.hist.quote(instrument="^gspc",
                 start="1994-02-11",end="2009-12-31", quote="AdjClose", quiet=TRUE)

# log-returns
y <- diff(log(p))
y <- coredata(y)

# prelim
T <- length(y)
WE <- 1000
p <- 0.01                       # tail prob (keep name 'p')
l1 <- ceiling((WE + 1) * p)     # HS index (type-1)
value <- 1
VaR <- matrix(nrow=T, ncol=4); colnames(VaR) <- c("EWMA","MA","HS","GARCH")

# EWMA setup
lambda <- 0.94
s11 <- var(y[1:30])
for (t in 2:WE) s11 <- lambda*s11 + (1-lambda)*y[t-1]^2

# helper: excess kurtosis -> df for t
excess_kurt <- function(x){
  m <- mean(x); s <- sd(x)
  if (!is.finite(s) || s == 0) return(0)
  mean(((x - m)/s)^4) - 3
}
df_from_kurt <- function(k){
  if (!is.finite(k) || k <= 0) return(1000)  # ~normal if no/heavy neg kurtosis
  min(1000, 4 + 6 / k)                       # t-kurtosis: 6/(nu-4)
}

library(fGarch)

# rolling backtest
for (t in (WE+1):T){
  t1 <- t-WE; t2 <- t-1
  window <- y[t1:t2]
  
  # update EWMA sigma
  s11 <- lambda*s11 + (1-lambda)*y[t-1]^2
  sigma_ewma <- sqrt(s11)
  
  # ---- Heavier tails for EWMA & MA (Student-t quantile) ----
  k_win <- excess_kurt(window)
  nu_hat <- df_from_kurt(k_win)
  q_t <- qt(p, df = nu_hat)                  # negative number
  
  VaR[t, "EWMA"] <- -sigma_ewma * q_t * value
  VaR[t, "MA"]   <- -sd(window) * q_t * value
  
  # HS (unchanged)
  ys <- sort(window)
  VaR[t, "HS"] <- -ys[l1] * value
  
  # ---- GARCH(1,1) with t-innovations; use model's 1-step sigma & t-quantile ----
  g <- garchFit(~ garch(1,1), data = window, include.mean = FALSE,
                cond.dist = "std", trace = FALSE)
  
  # 1) preferred: use predict() as a data.frame
  pred <- predict(g, n.ahead = 1)
  if ("standardDeviation" %in% names(pred)) {
    sig1 <- as.numeric(pred$standardDeviation[1])
  } else if ("sigma.t" %in% names(pred)) {         # some versions use 'sigma.t'
    sig1 <- as.numeric(pred$sigma.t[1])
  } else {
    # 2) robust fallback: one-step variance via recursion
    pars <- coef(g)
    h_t  <- tail(g@h.t, 1)
    eps2 <- tail(window, 1)^2
    h1   <- pars["omega"] + pars["alpha1"] * eps2 + pars["beta1"] * h_t
    sig1 <- sqrt(h1)
  }
  
  nu_g <- as.numeric(coef(g)["shape"])
  if (!is.finite(nu_g)) nu_g <- 8
  VaR[t, "GARCH"] <- -qstd(p, mean = 0, sd = sig1, nu = nu_g) * value
  
}

W1 <- WE+1
for (i in 1:4){
  VR <- sum(y[W1:T] < -VaR[W1:T, i]) / (p * (T - WE))
  s  <- sd(VaR[W1:T, i])
  cat(i, colnames(VaR)[i], "VR", round(VR,3), "VaR vol", round(s,6), "\n")
}

# plots
matplot(cbind(y[W1:T], -VaR[W1:T, ]), type='l', lty=1,
        col=c("black","turquoise","green","blue","orange"),
        lwd=c(1.6,1.2,1.2,1.2,1.2), xlab="Zeit", ylab="log-Returns / -VaR (1%)")
legend("bottomleft",
       legend=c("Log-Returns","VaR EWMA (t)","VaR MA (t)","VaR HS","VaR GARCH-t"),
       col=c("black","turquoise","green","blue","orange"), lty=1,
       lwd=c(1.6,1.2,1.2,1.2,1.2), bty="n")
abline(h=0, col="grey80")

# Kupiec & Christoffersen (unchanged)
bern_test=function (p, v) {
  a=p^(sum(v))*(1-p)^(length(v)-sum(v))
  b=(sum(v)/length(v))^(sum(v))*(1-(sum(v)/length(v)))^(length(v)-sum(v))
  -2*log(a/b)
}
ind_test=function (V) {
  J=matrix(ncol=4, nrow=length(V))
  for (i in 2:length(V)){
    J[i,1]=(V[i-1]==0 & V[i]==0)
    J[i,2]=(V[i-1]==0 & V[i]==1)
    J[i,3]=(V[i-1]==1 & V[i]==0)
    J[i,4]=(V[i-1]==1 & V[i]==1)
  }
  V_00=sum(J[,1],na.rm=TRUE); V_01=sum(J[,2],na.rm=TRUE)
  V_10=sum(J[,3],na.rm=TRUE); V_11=sum(J[,4],na.rm=TRUE)
  p_00=V_00/(V_00+V_01); p_01=V_01/(V_00+V_01)
  p_10=V_10/(V_10+V_11); p_11=V_11/(V_10+V_11)
  hat_p=(V_01+V_11)/(V_00+V_01+V_10+V_11)
  a=(1-hat_p)^(V_00+V_10)*(hat_p)^(V_01+V_11)
  b=(p_00)^(V_00)*(p_01)^(V_01)*(p_10)^(V_10)*p_11^(V_11)
  -2*log(a/b)
}

ya <- y[W1:T]; VaRa <- VaR[W1:T,]
m <- colnames(VaR)
for (i in 1:4){
  q <- ya < -VaRa[,i]
  v <- 0*ya; v[q] <- 1
  ber <- bern_test(p, v); ind <- ind_test(v)
  cat(m[i], "Bernoulli", ber, 1-pchisq(ber,1),
      "independence", ind, 1-pchisq(ind,1), "\n")
}

# Exceedance counts
E <- (ya < -VaRa); exceed_counts <- colSums(E, na.rm=TRUE)
cat("Exceedances by model (out of", (T-WE), "; expected ≈", p*(T-WE), "):\n")
for (j in 1:4) cat(m[j], ":", exceed_counts[j], "\n")
