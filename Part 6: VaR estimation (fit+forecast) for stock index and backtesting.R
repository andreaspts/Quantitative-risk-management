#based on stock index data VaR are calculated in 4 ways
#all in sample
#plots of log returns vs. respective VaRs computed
#backtesting to check if VaR exceedances are Bernoulli distributed and independent
#exceedance count to get a feeling if evaluation seems reasoable

library("tseries")
library("zoo")  
p=get.hist.quote(instrument="^gspc",
                 start="1994-02-11",end="2009-12-31", quote="AdjClose", quiet=TRUE)

# download prices
y=diff(log(p))
y=coredata(y)

#VaR estimation and backtesting
#preliminary work
T=length (y)
#estimation window
WE=1000
p=0.01
l1=ceiling((WE + 1) * p)#WE*p 
value=1;
VaR=matrix(nrow=T, ncol=4)

# EWMA setup
lambda= 0.94;
s11=var(y[1:30]);

for(t in 2:WE) {
  s11=lambda*s11+(1-lambda)*y[t-1]^2
}

#backtesting
library (fGarch)

for(t in (WE+1):T){
  t1=t-WE;
  t2=t-1;
  window=y[t1:t2]
  
  # EWMA
  s11=lambda*s11+ (1-lambda) *y [t-1] ^2
  VaR[t, 1] =-qnorm(p)*sqrt(s11) *value
  
  # MA
  VaR[t, 2] =-sd(window)*qnorm(p)*value
  
  # HS (historical simulation)
  ys=sort(window) # sort returns
  VaR[t, 3]=-ys[l1] *value # Calculate the level of VaR
  
  # GARCH (1, 1)
  g=garchFit (formula=~garch (1,1), window,trace=FALSE,
              include.mean=FALSE)
  par=g@fit$matcoef
  s4=par[1]+par[2]*window[WE]^2+par[3]*g@h.t[WE]
  VaR[t,4]=-sqrt(s4)*qnorm(p)*value
}

W1=WE+1
for (i in 1:4){
  VR = sum(y[W1:T] < -VaR[W1:T, i]) / (p * (T - WE))
  s  = sd(VaR[W1:T, i])
  cat (i, "VR", VR, "VaR vol", s, "\n")
}

#Simple plot
matplot(cbind (y[W1:T],-VaR [W1:T,]),type='l') # but still have to backtest and check for exceedances
## Plot with colors and legend
matplot(
  cbind(y[W1:T], -VaR[W1:T, ]),
  type = "l", lty = 1,
  col  = c("black", "turquoise", "green", "blue", "orange"),
  lwd  = c(1.6, 1.2, 1.2, 1.2, 1.2),
  xlab = "Zeit",
  ylab = "log-Returns / -VaR (1%)"
)
legend(
  "bottomleft",
  legend = c("Log-Returns",
             "VaR EWMA (1%)",
             "VaR MA (1%)",
             "VaR HS (1%)",
             "VaR GARCH(1,1) (1%)"),
  col = c("black", "turquoise", "green", "blue", "orange"),
  lty = 1, lwd = c(1.6, 1.2, 1.2, 1.2, 1.2),
  bty = "n"
)
abline(h = 0, col = "grey80")

#Bernoulli backtesting (want to check if exceedances are Bernoulli distributed)
bern_test=function (p, v) {
  a=p^(sum(v))*(1-p)^(length(v)-sum(v))
  b=(sum(v)/length(v))^(sum(v)) * (1-(sum(v)/length(v)))^(length(v)-sum(v))
  return(-2*log(a/b))
}

#Independence coverage
ind_test=function (V) {
  J=matrix(ncol=4, nrow=length(V))
  for (i in 2: length (V)){
    J[i, 1] =V[i-1]==0 & V[i]==0
    J[i, 2] =V[i-1] ==0 & V[i]==1
    J[i,3]=V[i-1]==1 & V[i]==0
    J[i,4]=V[i-1]==1 & V[i]==1
  }
  
  V_00=sum(J[,1], na.rm=TRUE)
  V_01=sum(J[,2], na.rm=TRUE)
  V_10=sum(J[,3], na.rm=TRUE)
  V_11=sum(J[,4], na.rm=TRUE)
  p_00=V_00 / (V_00+V_01)
  p_01=V_01/ (V_00+V_01) 
  p_10=V_10/ (V_10+V_11)
  p_11=V_11/ (V_10+V_11)
  hat_p=(V_01+V_11)/(V_00+V_01+V_10+V_11)
  a=(1-hat_p)^(V_00+V_10)*(hat_p)^(V_01+V_11)
  b= (p_00)^(V_00)*(p_01)^(V_01)*(p_10)^(V_10)*p_11^(V_11)
  return(-2*log(a/b))
}

W1=WE+1
ya=y[W1: T]
VaRa=VaR[W1: T,]
m=c("EWMA", "MA", "HS", "GARCH" )
for (i in 1:4){
  q = y[W1:T] < -VaR[W1:T, i]
  v = VaRa*0
  v[q, i] = 1
  ber=bern_test(p, v[,i])
  ind=ind_test(v[,i])
  cat(i, m[i], 'Bernoulli', ber, 1-pchisq(ber, 1),
      "independence", ind, 1-pchisq(ind, 1), "\n")
}

# --- Exceedance counts ---
E <- (y[W1:T] < -VaR[W1:T, ])          # logical matrix: breach if return < -VaR
exceed_counts <- colSums(E, na.rm=TRUE)
expected <- p * (T - WE)

colnames(VaR) <- c("EWMA","MA","HS","GARCH")  # in case not set
cat("Exceedances by model (out of", (T-WE), "days; expected ≈", expected, "):\n")
for (j in 1:4) {
  cat(colnames(VaR)[j], ":", exceed_counts[j], "\n")
}
