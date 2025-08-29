#plot danish fire data

#load QRM library
library(QRM) 

plot(danish,
     xlab = "Year",
     ylab = "Insurance claims in million DKK")


#flag unusually large claims in the Danish fire losses data

#1. Threshold-based detection
#idea: Picks a high quantile or 
#an absolute value cutoff and marks points above it.

data(danish)
y <- as.numeric(danish)
dates <- time(danish)

# Example: threshold at 99th percentile
thr <- quantile(y, 0.99)
high_idx <- y > thr

plot(dates, y, type = "h",
     xlab = "Year", ylab = "Insurance claims (million DKK)",
     main = "Danish Fire Losses")
points(dates[high_idx], y[high_idx], col = "red", pch = 19)
abline(h = thr, col = "blue", lty = 2)

#2. Median Absolute Deviation (MAD) method
#idea: Identifies outliers.

med <- median(y)
mad_val <- mad(y)
high_idx <- y > med + 5 * mad_val  # "5" = sensitivity factor

plot(dates, y, type = "h",
     main = "Outliers via MAD",
     xlab = "Year", ylab = "Insurance claims (million DKK)")
points(dates[high_idx], y[high_idx], col = "red", pch = 19)

#3. EVT-based approach (fits a Generalized Pareto Distribution to tail)
#idea: Gives statistical parameters for the tail 
#and helps to define "unusually high".

library(evir)
gpd_fit <- gpd(y, threshold = thr)  # Fit GPD above 99% quantile
gpd_fit

# ===== Tail analysis additions =====

# Make sure we have y (in million DKK) and dates (fractional years) from above:
# y <- as.numeric(danish); dates <- time(danish)

# --- 1) Count of claims > 1 million DKK between 1980 and 1990 (inclusive) ---
idx_8090 <- (dates >= 1980) & (dates <= 1990)
count_over_1m <- sum(y[idx_8090] > 1, na.rm = TRUE)
cat("\nCount of Danish fire claims > 1 million DKK in [1980, 1990]:",
    count_over_1m, "\n")

# --- 2) POT threshold at 10 million DKK; fit GPD by MLE ---
# We'll use evd::fpot which does MLE and returns 'scale' and 'shape'
suppressWarnings(require(evd))
u <- 10  # threshold in million DKK
fit_gpd <- fpot(y, threshold = u)   # MLE fit on exceedances y > u
gpd_pars <- fit_gpd$estimate
sigma_u  <- as.numeric(gpd_pars["scale"])  # GPD scale at threshold u
xi_hat   <- as.numeric(gpd_pars["shape"])  # GPD shape (tail index)

# also print basic fit info
vc  <- tryCatch(vcov(fit_gpd), error = function(e) NULL)
se  <- if (!is.null(vc)) sqrt(diag(vc)) else rep(NA_real_, 2)
se_sigma <- se[match("scale", names(gpd_pars))]
se_xi    <- se[match("shape", names(gpd_pars))]
gpd_tbl <- data.frame(
  parameter = c("sigma_u (scale)", "xi (shape)"),
  estimate  = c(sigma_u, xi_hat),
  se        = c(se_sigma, se_xi),
  lower95   = c(sigma_u, xi_hat) - 1.96 * c(se_sigma, se_xi),
  upper95   = c(sigma_u, xi_hat) + 1.96 * c(se_sigma, se_xi)
)
num_cols <- sapply(gpd_tbl, is.numeric)
gpd_tbl[num_cols] <- lapply(gpd_tbl[num_cols], function(v) round(v, 6))
cat("\nGPD MLE at u = 10 million DKK:\n")
print(gpd_tbl, row.names = FALSE)

# beta is the scale parameter from fpot()
beta_hat <- sigma_u
cat("\nbeta (scale at threshold u) =", round(beta_hat, 6), "\n")
cat("\nxi (form at threshold u) =", round(xi_hat, 6), "\n")

# --- 3) Plot: empirical distribution (staircase) of excess losses vs fitted GPD CDF ---
excess <- y[y > u] - u   # exceedances X - u  (X in million DKK)
Fhat   <- ecdf(excess)   # empirical CDF (staircase)

# GPD CDF function
pgpd_cdf <- function(z, xi, beta) {
  z <- pmax(z, 0)
  if (abs(xi) < 1e-8) {
    1 - exp(-z / beta)
  } else {
    pmax(0, pmin(1, 1 - (1 + xi * z / beta)^(-1/xi)))
  }
}

# Plot CDFs
plot(Fhat, do.points = FALSE, verticals = TRUE,
     main = "Excesses above 10m DKK: Empirical CDF vs fitted GPD CDF",
     xlab = "Excess (million DKK above u = 10)", ylab = "CDF",
     col = "grey30", lwd = 1.5)
zgrid <- seq(0, max(excess, na.rm = TRUE)*1.05, length.out = 600)
lines(zgrid, pgpd_cdf(zgrid, xi_hat, sigma_u), col = "tomato", lwd = 2)
rug(excess, col = "grey70")
legend("bottomright",
       legend = c("Empirical CDF (exceedances)", "Fitted GPD CDF", "Exceedances (rug)"),
       col    = c("grey30", "tomato", "grey70"),
       lty    = c(1, 1, NA), lwd = c(1.5, 2, NA), pch = c(NA, NA, 124), bty = "n")


# ===== 4) VaR and ES from the fitted GPD (POT at u = 10) =====================

# Exceedance rate above the threshold
n_tot <- length(y)
n_exc <- sum(y > u, na.rm = TRUE)
p_u   <- n_exc / n_tot  # P(Y > u)

cat(sprintf("\nExceedance rate above u=%.2f: %d / %d = %.4f\n",
            u, n_exc, n_tot, p_u))

# Vector of confidence levels to report (feel free to change)
probs <- c(0.99, 0.995, 0.999)

# Helper: VaR from POT-GPD for p in probs (only valid if (1-p) < p_u)
gpd_var <- function(p, u, xi, beta, p_u) {
  # returns NA if quantile is below threshold
  below <- (1 - p) >= p_u
  out   <- rep(NA_real_, length(p))
  ok    <- !below
  
  if (any(ok)) {
    if (abs(xi) < 1e-8) {
      # xi -> 0 (Exponential limit)
      z <- beta * log(p_u / (1 - p[ok]))
    } else {
      z <- (beta/xi) * ((p_u / (1 - p[ok]))^xi - 1)
    }
    out[ok] <- u + z
  }
  out
}

# Helper: ES from POT-GPD (finite only if xi < 1); same applicability as VaR
gpd_es <- function(p, u, xi, beta, p_u) {
  VaR <- gpd_var(p, u, xi, beta, p_u)
  out <- rep(NA_real_, length(p))
  ok  <- is.finite(VaR) & (xi < 1)
  
  if (any(ok)) {
    z <- VaR[ok] - u
    if (abs(xi) < 1e-8) {
      # xi -> 0: ES = u + (beta + z)
      out[ok] <- u + (beta + z)
    } else {
      out[ok] <- u + (beta + z) / (1 - xi)
    }
  }
  # xi >= 1 => ES is infinite; mark as Inf where VaR is above threshold
  inf_idx <- is.finite(VaR) & (xi >= 1)
  out[inf_idx] <- Inf
  
  out
}

VaR_p <- gpd_var(probs, u = u, xi = xi_hat, beta = sigma_u, p_u = p_u)
ES_p  <- gpd_es (probs, u = u, xi = xi_hat, beta = sigma_u, p_u = p_u)

# Nicely formatted output (units: million DKK)
out_tbl <- data.frame(
  prob = probs,
  VaR_million_DKK = round(VaR_p, 4),
  ES_million_DKK  = round(ES_p, 4),
  note = ifelse((1 - probs) >= p_u,
                "below threshold u (POT not used)", "")
)
cat("\nPOT–GPD Risk Measures (units: million DKK)\n")
print(out_tbl, row.names = FALSE)

# Optional: quick message for any probs below threshold
if (any((1 - probs) >= p_u)) {
  cat("\nNote: For entries marked 'below threshold u', the requested quantile lies",
      "in the body (not the POT tail). Consider higher p or a lower threshold.\n")
}
