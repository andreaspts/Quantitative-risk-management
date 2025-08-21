# Download the index prices (S&P 500)
library(tseries)
library(zoo)  # for plot.zoo

p <- get.hist.quote(instrument="^gspc",
                    start="2005-01-01", end="2009-12-31",
                    quote="AdjClose", quiet=TRUE)

# Determination of log returns (in %), drop initial NA from diff
y <- na.omit(diff(log(p)) * 100)

# Clean the data from its mean value (demeaned returns)
Y <- y - mean(y)

# --- Plots: price, returns y, demeaned returns Y --------------------------
op <- par(mfrow = c(3,1), mar = c(4,4,3,1))

# (1) Closing price
plot.zoo(p, col = "black", lwd = 1.6, xlab = "", ylab = "Index level",
         main = "S&P 500 — Adjusted Close")
grid(col = "grey80", lty = 1, lwd = 0.8)

# (2) Daily log returns (in %)
plot.zoo(y, col = "steelblue", lwd = 1.4, xlab = "", ylab = "Return (%)",
         main = "S&P 500 — Daily Log Returns (in %)")
abline(h = 0, col = "grey50")
grid(col = "grey80", lty = 1, lwd = 0.8)

# (3) Demeaned daily log returns
plot.zoo(Y, col = "firebrick", lwd = 1.4, xlab = "", ylab = "Return (%)",
         main = "S&P 500 — Demeaned Daily Log Returns (Y = y − mean(y))")
abline(h = 0, col = "grey50")
grid(col = "grey80", lty = 1, lwd = 0.8)

par(op)

#Fitting GARCH models to this:
#install.packages("fGarch")
library(fGarch)

garchFit(~garch(1, 0),data=y,include.mean=FALSE)
garchFit(~garch(4,0),data=y,include.mean=FALSE)
garchFit(~garch(4,1),data=y,include.mean=FALSE)
garchFit(~garch(1,1),data=y,include.mean=FALSE)
garchFit(~garch(1,1),data=y,include.mean=FALSE,cond.dist="std", trace=F) #student t distribution

# Save output in res

res=garchFit(~garch (1,1),data=y, include.mean=FALSE,cond.dist="sstd", trace=F) #skewed student t distribution

plot (res)
