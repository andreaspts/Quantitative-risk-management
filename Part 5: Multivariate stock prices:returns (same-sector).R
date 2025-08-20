# install.packages(c("quantmod","vars"))  # if needed
library(quantmod); library(vars); library(zoo)

tickers <- c("GM","BMW.DE","MBG.DE")      # Mercedes: MBG.DE on XETRA
getSymbols(tickers, src="yahoo", from="2018-01-01", auto.assign=TRUE)

# Adjusted close prices merged
P <- na.omit(merge(Ad(GM), Ad(BMW.DE), Ad(MBG.DE)))
colnames(P) <- c("GM","BMW","MBG")

# Daily log returns
R <- na.omit(diff(log(P)))                # stationary

# ---------- NEW: plot the 3 price series in one plot ----------
cols <- c("firebrick","dodgerblue3","darkgreen")
op <- par(mfrow = c(1,1), mar = c(4,4,3,1))
plot.zoo(P, plot.type = "single",
         col = cols, lwd = 1.6,
         xlab = "", ylab = "Adjusted close",
         main = "GM, BMW, MBG — Adjusted Closing Prices")
grid(col = "grey80", lty = 1, lwd = 0.8)
legend("topleft", bty = "n", lwd = 1.6, col = cols, legend = colnames(P))
par(op)

# ---------- NEW: ACF & PACF of returns for each stock ----------
op <- par(mfrow = c(3, 2), mar = c(4,4,3,1))

acf(as.numeric(R$GM),  lag.max = 40, main = "ACF: GM returns");  grid(col="grey80", lty=1, lwd=0.8)
pacf(as.numeric(R$GM), lag.max = 40, main = "PACF: GM returns"); grid(col="grey80", lty=1, lwd=0.8)

acf(as.numeric(R$BMW),  lag.max = 40, main = "ACF: BMW returns");  grid(col="grey80", lty=1, lwd=0.8)
pacf(as.numeric(R$BMW), lag.max = 40, main = "PACF: BMW returns"); grid(col="grey80", lty=1, lwd=0.8)

acf(as.numeric(R$MBG),  lag.max = 40, main = "ACF: MBG returns");  grid(col="grey80", lty=1, lwd=0.8)
pacf(as.numeric(R$MBG), lag.max = 40, main = "PACF: MBG returns"); grid(col="grey80", lty=1, lwd=0.8)

par(op)

# ---------- Existing quick look ----------
pairs(coredata(R))                        # scatter matrix
print(cor(R))                             # correlation matrix

# ---------- Existing VAR on returns ----------
sel <- VARselect(R, lag.max = 10, type = "const")
fit <- VAR(R, p = sel$selection["AIC(n)"], type = "const")
summary(fit)
