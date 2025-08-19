# install.packages("quantmod")  # uncomment if needed
library(quantmod)
library(zoo)

# --- Fetch S&P 500 closing prices ---
start_date <- "2025-01-01"
getSymbols("^GSPC", src = "yahoo", from = start_date, auto.assign = TRUE)
px_close <- Cl(GSPC)

# --- Log returns from closing prices ---
ret <- na.omit(diff(log(px_close)))

# --- Moving averages (window 2 and 10) ---
MA2_px  <- TTR::SMA(px_close, n = 2)
MA10_px <- TTR::SMA(px_close, n = 10)
MA2_r   <- TTR::SMA(ret,      n = 2)
MA10_r  <- TTR::SMA(ret,      n = 10)

# Align series
PLOT_PX  <- na.omit(merge(px_close, MA2_px, MA10_px, all = FALSE))
PLOT_RET <- na.omit(merge(ret,      MA2_r,  MA10_r,  all = FALSE))

# --- Two-panel plot ----------------------------------------------------------
op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

## (1) Closing prices + MA(2) + MA(10)
plot.zoo(PLOT_PX, plot.type = "single",
         col = c("black", "dodgerblue3", "firebrick"),
         lwd = c(1.4, 1.6, 1.6),
         xlab = "", ylab = "Index level",
         main = "S&P 500 Close with MA(2) & MA(10)")
grid(col = "grey80", lty = 1, lwd = 0.8)   # solid gray lines
legend("topleft", bty = "n", lwd = c(1.4, 1.6, 1.6),
       col = c("black", "dodgerblue3", "firebrick"),
       legend = c("Close", "MA(2)", "MA(10)"))

## (2) Log returns + MA(2) + MA(10)
plot.zoo(PLOT_RET, plot.type = "single",
         col = c("grey30", "dodgerblue3", "firebrick"),
         lwd = c(1.4, 1.6, 1.6),
         xlab = "", ylab = "Log return",
         main = "S&P 500 Log Returns with MA(2) & MA(10)")
grid(col = "grey80", lty = 1, lwd = 0.8)   # solid gray lines
legend("topleft", bty = "n", lwd = c(1.4, 1.6, 1.6),
       col = c("grey30", "dodgerblue3", "firebrick"),
       legend = c("Log return", "MA(2)", "MA(10)"))

par(op)
