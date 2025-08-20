# install.packages(c("quantmod","vars"))  # if needed
library(quantmod)
library(vars)

getSymbols("AAPL", src = "yahoo", from = "2008-01-01")

px  <- Ad(AAPL)   # adjusted close
vol <- Vo(AAPL)   # volume

# -------- NEW: Dual-axis plot of raw Price (left) and Volume (right) ------
RAW <- na.omit(merge(px, vol))
idx_raw <- index(RAW)

op <- par(mfrow = c(1,1), mar = c(4, 4, 3, 4) + 0.1)
plot(idx_raw, as.numeric(RAW[, 1]),
     type = "l", col = "black", lwd = 1.4,
     xlab = "", ylab = "Adjusted close",
     main = "AAPL: Adjusted Close (left) vs Volume (right)")
grid(col = "grey80", lty = 1, lwd = 0.8)

par(new = TRUE)
plot(idx_raw, as.numeric(RAW[, 2]),
     type = "l", col = "darkorange3", lwd = 1.2,
     axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = "darkorange3")
mtext("Volume", side = 4, line = 3, col = "darkorange3")
legend("topleft", bty = "n", lwd = c(1.4, 1.2),
       col = c("black", "darkorange3"),
       legend = c("Adjusted close", "Volume"))
par(op)

# -------- Logged series ----------------------------------------------------
ret     <- na.omit(diff(log(px)))     # log returns
dlogVol <- na.omit(diff(log(vol)))    # change in log volume

# Align in one object (drop any remaining NA mismatch)
X <- na.omit(merge(ret, dlogVol))
colnames(X) <- c("ret", "dlogVol")

# -------- Dual-axis plot: returns (left) vs Δlog volume (right) -----------
idx <- index(X)

op <- par(mfrow = c(1,1), mar = c(4, 4, 3, 4) + 0.1)
plot(idx, as.numeric(X$ret),
     type = "l", col = "steelblue", lwd = 1.4,
     xlab = "", ylab = "Log return",
     main = "AAPL: Log returns (left) vs Δlog(volume) (right)")
grid(col = "grey80", lty = 1, lwd = 0.8)

par(new = TRUE)
plot(idx, as.numeric(X$dlogVol),
     type = "l", col = "firebrick", lwd = 1.4,
     axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = "firebrick")
mtext("Δ log(volume)", side = 4, line = 3, col = "firebrick")
legend("topleft", bty = "n", lwd = 1.4,
       col = c("steelblue", "firebrick"),
       legend = c("Log return", "Δ log(volume)"))
par(op)

# -------- ACF & PACF for both series (2x2 grid) ---------------------------
op <- par(mfrow = c(2, 2), mar = c(4,4,3,1))

acf(as.numeric(X$ret), lag.max = 40, main = "ACF: log returns")
grid(col = "grey80", lty = 1, lwd = 0.8)
pacf(as.numeric(X$ret), lag.max = 40, main = "PACF: log returns")
grid(col = "grey80", lty = 1, lwd = 0.8)

acf(as.numeric(X$dlogVol), lag.max = 40, main = "ACF: Δ log(volume)")
grid(col = "grey80", lty = 1, lwd = 0.8)
pacf(as.numeric(X$dlogVol), lag.max = 40, main = "PACF: Δ log(volume)")
grid(col = "grey80", lty = 1, lwd = 0.8)

par(op)

# -------- Small VAR on [ret, dlogVol] -------------------------------------
sel <- VARselect(X, lag.max = 10, type = "const")
fit <- VAR(X, p = sel$selection["AIC(n)"], type = "const")
summary(fit)

