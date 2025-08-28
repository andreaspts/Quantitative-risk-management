# === NASDAQ (^IXIC): Close+monthly max (plot 1), -Log returns+monthly max (plot 2),
#                     then GEV fit to -monthly_max_ret ============================

# install.packages(c("quantmod","evd"))  # <- run once if needed
library(quantmod)
library(evd)

## --- Data ---
ixic <- getSymbols("^IXIC", src = "yahoo",
                   from = "1994-01-01", to = "2005-01-01",#Sys.Date(),
                   auto.assign = FALSE)

cl  <- Cl(ixic)
ret <- na.omit(diff(log(cl)))  # daily log returns

## Monthly maxima of closing data
monthly_max_cl  <- apply.monthly(cl,  function(x) max(x, na.rm = TRUE))
## Monthly maxima of negative log returns (maxima of losses)
monthly_max_ret <- apply.monthly(-ret, function(x) max(x, na.rm = TRUE))

## --- Two stacked plots (only 2 figures total) ---
op <- par(mfrow = c(2,1), mar = c(4,4,3,1))

# (1) Close + monthly max(Close)
plot.zoo(merge(cl, monthly_max_cl),
         plot.type = "single",
         col = c("grey40","tomato"),
         lwd = c(1,2),
         type = c("l","o"),     # line for close, line+points for monthly max
         pch  = c(NA,16),
         cex  = c(NA,0.6),
         xlab = "Date", ylab = "Index Level",
         main = "NASDAQ (^IXIC): Daily Close with Monthly Max(Close)")
legend("topleft", c("Daily Close","Monthly Max (Close)"),
       col = c("grey40","tomato"), lty = 1, pch = c(NA,16), lwd = c(1,2), bty = "n")

# (2) Log returns + monthly max(Log return)
plot.zoo(merge(ret, monthly_max_ret),
         plot.type = "single",
         col = c("grey60","royalblue"),
         lwd = c(1,2),
         type = c( "l","o"),
         pch  = c(NA,16),
         cex  = c(NA,0.6),
         xlab = "Date", ylab = "-Log Return (daily)",
         main = "NASDAQ (^IXIC): Daily -Log Returns with Monthly Max(Return)")
abline(h = 0, col = "grey80")
legend("topleft", c("Daily -Log Return","Monthly Max (Log Return)"),
       col = c("grey60","royalblue"), lty = 1, pch = c(NA,16), lwd = c(1,2), bty = "n")

par(op)

## --- GEV fit to monthly_max_ret (block maxima of daily -log returns) ---
x <- as.numeric(coredata(monthly_max_ret))
x <- x[is.finite(x)]

fit_gev <- fgev(x)  # MLE

# Parameters μ, σ, ξ (+ SE and 95% CI if available)
pars  <- fit_gev$estimate
mu    <- as.numeric(pars["loc"])
sigma <- as.numeric(pars["scale"])
xi    <- as.numeric(pars["shape"])

vc  <- tryCatch(vcov(fit_gev), error = function(e) NULL)
se  <- if (!is.null(vc)) sqrt(diag(vc)) else rep(NA_real_, 3)
se_mu    <- se[match("loc",   names(pars))]
se_sigma <- se[match("scale", names(pars))]
se_xi    <- se[match("shape", names(pars))]

param_tbl <- data.frame(
  parameter = c("mu (loc)", "sigma (scale)", "xi (shape)"),
  estimate  = c(mu, sigma, xi),
  se        = c(se_mu, se_sigma, se_xi),
  lower95   = c(mu, sigma, xi) - 1.96 * c(se_mu, se_sigma, se_xi),
  upper95   = c(mu, sigma, xi) + 1.96 * c(se_mu, se_sigma, se_xi),
  stringsAsFactors = FALSE
)
num_cols <- sapply(param_tbl, is.numeric)
param_tbl_print <- param_tbl
param_tbl_print[num_cols] <- lapply(param_tbl_print[num_cols], function(v) round(v, 6))
cat("\nGEV parameters (monthly maxima of daily -log returns):\n")
print(param_tbl_print, row.names = FALSE)

# Return levels (per block = per month)
period_months <- c(12, 60, 120, 240)   # ~1y, 5y, 10y, 20y
p_block <- 1 - 1/period_months
rl <- qgev(p_block, loc = mu, scale = sigma, shape = xi)
rl_tbl <- data.frame(period_months = period_months,
                     return_level  = round(rl, 6))
cat("\nGEV return levels (per month/block):\n")
print(rl_tbl, row.names = FALSE)

# Diagnostic plot: histogram + fitted GEV density
op2 <- par(mfrow = c(1,1), mar = c(4,4,3,1))
hist(x, breaks = "FD", freq = FALSE, col = "grey90", border = "white",
     main = "GEV fit to Monthly Max of Daily -Log Returns",
     xlab = "Monthly max (daily -log return)",
     xlim = c(-0.1, 0.2))                     # <-- extended x-axis range
xx <- seq(-0.2, 1.0, length.out = 400)        # <-- match the plotting range
lines(xx, dgev(xx, loc = mu, scale = sigma, shape = xi),
      lwd = 2, col = "tomato")
legend("topright", c("Empirical (hist)","Fitted GEV density"),
       lty = c(NA,1), lwd = c(NA,2), pch = c(15, NA),
       pt.cex = 1.2, col = c("grey70","tomato"), bty = "n")
par(op2)

#yields a Frechet distribution