# === NASDAQ (^IXIC): Close+monthly max (plot 1), -Log returns+monthly max (plot 2),
#                     then GEV fit to -monthly_max_ret + POT-GPD on daily losses ===

# install.packages(c("quantmod","evd"))  # <- run once if needed
library(quantmod)
library(evd)

## --- Data ---
ixic <- getSymbols("^IXIC", src = "yahoo",
                   from = "1994-01-01", to = "2025-01-01",
                   auto.assign = FALSE)

cl  <- Cl(ixic)
ret <- na.omit(diff(log(cl)))              # daily log returns
loss <- as.numeric(-ret)                    # daily "losses" (positive = down day)

## Monthly maxima
monthly_max_cl  <- apply.monthly(cl,  function(x) max(x, na.rm = TRUE))
monthly_max_ret <- apply.monthly(-ret, function(x) max(x, na.rm = TRUE))  # max monthly loss

## --- Two stacked plots (only 2 figures total) ---
op <- par(mfrow = c(2,1), mar = c(4,4,3,1))

# (1) Close + monthly max(Close)
plot.zoo(merge(cl, monthly_max_cl),
         plot.type = "single",
         col = c("grey40","tomato"),
         lwd = c(1,2),
         type = c("l","o"),
         pch  = c(NA,16),
         cex  = c(NA,0.6),
         xlab = "Date", ylab = "Index Level",
         main = "NASDAQ (^IXIC): Daily Close with Monthly Max(Close)")
legend("topleft", c("Daily Close","Monthly Max (Close)"),
       col = c("grey40","tomato"), lty = 1, pch = c(NA,16), lwd = c(1,2), bty = "n")

# (2) -Log returns + monthly max(-Log return)
plot.zoo(merge(-ret, monthly_max_ret),
         plot.type = "single",
         col = c("grey60","royalblue"),
         lwd = c(1,2),
         type = c("l","o"),
         pch  = c(NA,16),
         cex  = c(NA,0.6),
         xlab = "Date", ylab = "-Log Return (daily)",
         main = "NASDAQ (^IXIC): Daily -Log Returns with Monthly Max(-Return)")
abline(h = 0, col = "grey80")
legend("topleft", c("Daily -Log Return","Monthly Max (-Log Return)"),
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

# GEV diagnostic plot: histogram + fitted GEV density
op2 <- par(mfrow = c(1,1), mar = c(4,4,3,1))
hist(x, breaks = "FD", freq = FALSE, col = "grey90", border = "white",
     main = "GEV fit to Monthly Max of Daily -Log Returns",
     xlab = "Monthly max (daily -log return)")
xx <- seq(min(x), max(x), length.out = 400)
lines(xx, dgev(xx, loc = mu, scale = sigma, shape = xi),
      lwd = 2, col = "tomato")
legend("topright", c("Empirical (hist)","Fitted GEV density"),
       lty = c(NA,1), lwd = c(NA,2), pch = c(15, NA),
       pt.cex = 1.2, col = c("grey70","tomato"), bty = "n")
par(op2)

## --- POT (Peaks-Over-Threshold) GPD on daily losses ---------------------------
## Pickands–Balkema–de Haan: Excesses above a high threshold ~ GPD

# 1) Choose a high threshold (simple heuristic: 95% quantile of daily losses)
u <- as.numeric(quantile(loss, probs = 0.95, na.rm = TRUE))

# 2) Fit GPD to exceedances over u
fit_gpd   <- fpot(loss, threshold = u)   # evd::fpot
gpd_pars  <- fit_gpd$estimate            # scale & shape at threshold u
gpd_sigma <- as.numeric(gpd_pars["scale"])
gpd_xi    <- as.numeric(gpd_pars["shape"])

# SEs (if available)
gpd_vc <- tryCatch(vcov(fit_gpd), error = function(e) NULL)
gpd_se <- if (!is.null(gpd_vc)) sqrt(diag(gpd_vc)) else rep(NA_real_, 2)
se_sigma_gpd <- gpd_se[match("scale", names(gpd_pars))]
se_xi_gpd    <- gpd_se[match("shape", names(gpd_pars))]

n_tot <- length(loss)
n_exc <- sum(loss > u, na.rm = TRUE)
cat("\nPOT-GPD on daily losses (X = -log return):\n")
cat("  Threshold u =", round(u, 6),
    " | exceedances =", n_exc, "/", n_tot,
    sprintf("(%.2f%%)\n", 100 * n_exc / n_tot))

gpd_tbl <- data.frame(
  parameter = c("sigma_u (scale)", "xi (shape)"),
  estimate  = c(gpd_sigma, gpd_xi),
  se        = c(se_sigma_gpd, se_xi_gpd),
  lower95   = c(gpd_sigma, gpd_xi) - 1.96 * c(se_sigma_gpd, se_xi_gpd),
  upper95   = c(gpd_sigma, gpd_xi) + 1.96 * c(se_sigma_gpd, se_xi_gpd),
  stringsAsFactors = FALSE
)
num_cols <- sapply(gpd_tbl, is.numeric)
gpd_tbl_print <- gpd_tbl
gpd_tbl_print[num_cols] <- lapply(gpd_tbl_print[num_cols], function(v) round(v, 6))
print(gpd_tbl_print, row.names = FALSE)

# 3) Diagnostic: histogram of exceedances (Y = X - u), overlay GPD density
y_exc <- loss[loss > u] - u    # exceedances above u
op3 <- par(mfrow = c(1,1), mar = c(4,4,3,1))
hist(y_exc, breaks = "FD", freq = FALSE, col = "grey90", border = "white",
     main = "POT-GPD fit on Exceedances of Daily Losses",
     xlab = "Exceedance (loss - u)")
yy <- seq(0, max(y_exc, na.rm = TRUE), length.out = 400)

# Fitted GPD exceedance density
lines(yy, dgpd(yy, scale = gpd_sigma, shape = gpd_xi),
      lwd = 2, col = "royalblue")

# --- OVERLAY: GEV-implied conditional exceedance density f(Y=y | X>u)
# NOTE: GEV fitted to monthly maxima, not daily losses → interpret cautiously.
denom <- 1 - pgev(u, loc = mu, scale = sigma, shape = xi)
if (!is.finite(denom) || denom <= 0) denom <- .Machine$double.eps
gev_cond <- dgev(u + yy, loc = mu, scale = sigma, shape = xi) / denom
lines(yy, gev_cond, lwd = 2, col = "tomato", lty = 2)

legend("topright",
       c("Empirical exceedances", "Fitted GPD density", "GEV-implied cond. density"),
       lty = c(NA,1,2), lwd = c(NA,2,2), pch = c(15, NA, NA),
       pt.cex = 1.2, col = c("grey70","royalblue","tomato"), bty = "n")
par(op3)

#Note: While according to the theorem, the ξ values should match, 
#they don’t have to match exactly here.They should be compatible, 
#because both methods estimate the same tail index ξ 
#under the same data-generating process (Fisher–Tippett–Gnedenko & Pickands–Balkema–de Haan).
#The numbers obtained are consistent. 
#The intervals overlap widely, so statistically one can’t reject 
#that both procedures share the same ξ. 
#Practically, they point to at most mildly heavy tails.