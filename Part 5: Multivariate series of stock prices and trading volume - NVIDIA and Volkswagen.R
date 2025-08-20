# install.packages(c("quantmod","vars"))  # if needed
library(quantmod)
library(vars)

analyze_symbol <- function(symbol, start = "2008-01-01", label = symbol) {
  cat("\n==========", label, "(", symbol, ") ==========\n")
  
  # Fetch from Yahoo as a single xts object (no auto-assign side effects)
  XTS <- getSymbols(symbol, src = "yahoo", from = start, auto.assign = FALSE)
  px  <- Ad(XTS)   # adjusted close
  vol <- Vo(XTS)   # volume
  
  ## 1) Raw Price (left) vs Volume (right)
  RAW <- na.omit(merge(px, vol))
  idx_raw <- index(RAW)
  op <- par(mfrow = c(1,1), mar = c(4,4,3,4) + 0.1)
  plot(idx_raw, as.numeric(RAW[,1]),
       type = "l", col = "black", lwd = 1.4,
       xlab = "", ylab = "Adjusted close",
       main = sprintf("%s: Adjusted Close (left) vs Volume (right)", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  par(new = TRUE)
  plot(idx_raw, as.numeric(RAW[,2]),
       type = "l", col = "darkorange3", lwd = 1.2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4, col.axis = "darkorange3")
  mtext("Volume", side = 4, line = 3, col = "darkorange3")
  legend("topleft", bty = "n", lwd = c(1.4, 1.2),
         col = c("black", "darkorange3"),
         legend = c("Adjusted close", "Volume"))
  par(op)
  
  ## 2) Logged series
  ret     <- na.omit(diff(log(px)))     # log returns
  dlogVol <- na.omit(diff(log(vol)))    # Δ log(volume)
  X <- na.omit(merge(ret, dlogVol))
  colnames(X) <- c("ret", "dlogVol")
  idx <- index(X)
  
  # Dual-axis: returns (left) vs Δlog(volume) (right)
  op <- par(mfrow = c(1,1), mar = c(4,4,3,4) + 0.1)
  plot(idx, as.numeric(X$ret),
       type = "l", col = "steelblue", lwd = 1.4,
       xlab = "", ylab = "Log return",
       main = sprintf("%s: Log returns (left) vs Δlog(volume) (right)", label))
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
  
  # 3) ACF & PACF (2×2)
  op <- par(mfrow = c(2,2), mar = c(4,4,3,1))
  acf(as.numeric(X$ret), lag.max = 40, main = sprintf("ACF: log returns — %s", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(as.numeric(X$ret), lag.max = 40, main = sprintf("PACF: log returns — %s", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(as.numeric(X$dlogVol), lag.max = 40, main = sprintf("ACF: Δ log(volume) — %s", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(as.numeric(X$dlogVol), lag.max = 40, main = sprintf("PACF: Δ log(volume) — %s", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  par(op)
  
  # 4) Small VAR on [ret, dlogVol]
  sel <- VARselect(X, lag.max = 10, type = "const")
  p_sel <- sel$selection["AIC(n)"]
  fit <- VAR(X, p = p_sel, type = "const")
  cat(sprintf("%s — VAR lag chosen by AIC: p = %s\n", label, p_sel))
  print(summary(fit))
}

# ---- Run for NVIDIA and Volkswagen ----
analyze_symbol("NVDA",    start = "2008-01-01", label = "NVIDIA")
analyze_symbol("VOW3.DE", start = "2008-01-01", label = "Volkswagen (VOW3.DE)")
