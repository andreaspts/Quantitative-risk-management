# install.packages("quantmod")   # uncomment if needed
library(quantmod)

# --- fetch S&P 500 levels ---
start_date <- "2010-01-01"
getSymbols("^GSPC", src = "yahoo", from = start_date, auto.assign = TRUE)

# --- daily log returns ---
ret_xts <- na.omit(dailyReturn(Ad(GSPC), type = "log"))
r <- as.numeric(ret_xts)  # numeric vector for acf/pacf/qq

# --- time series plot of returns ---
plot(ret_xts, main = "S&P 500: Daily Log Returns",
     ylab = "log return", xlab = "", col = "steelblue")

# --- Diagnostics: ACF, PACF, QQ (returns) + ACF, PACF (squared returns) ---
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Returns: ACF, PACF, QQ
acf(r,  lag.max = 40, main = "ACF: returns")
pacf(r, lag.max = 40, main = "PACF: returns")
qqnorm(r, main = "Normal QQ: returns"); qqline(r, col = "red")

# Squared returns: ACF, PACF (volatility clustering)
acf(r^2,  lag.max = 40, main = "ACF: squared returns")
pacf(r^2, lag.max = 40, main = "PACF: squared returns")

# Bonus panel: histogram with kernel density (heavy tails check)
hist(r, breaks = 60, probability = TRUE,
     main = "Histogram: returns", xlab = "log return")
lines(density(r, na.rm = TRUE))

par(op)
