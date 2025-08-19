#For AR/MA models, dependence shows up in the series ACF/PACF.
#For ARCH/GARCH, linear correlations of returns are (near) zero; 
#the dependence is in the volatility, so we examine ACF/PACF of squared (or absolute) returns.

set.seed(123)
n    <- 60      # length of each series
lags <- 40       # ACF/PACF max lag

# ---------- helpers ----------
triptych <- function(x, label, lags = 40, on_squared = FALSE) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "value",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  y <- if (on_squared) x^2 else x
  acf(y,  lag.max = lags, main = paste("ACF —", if (on_squared) "squared" else "series", "—", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(y, lag.max = lags, main = paste("PACF —", if (on_squared) "squared" else "series", "—", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
}

simulate_arch1 <- function(n, omega, alpha, burn = 500) {
  stopifnot(omega > 0, alpha >= 0, alpha < 1)
  z <- rnorm(n + burn)
  eps <- sigma2 <- numeric(n + burn)
  sigma2[1] <- omega / (1 - alpha)
  for (t in 2:(n + burn)) {
    sigma2[t] <- omega + alpha * eps[t - 1]^2
    eps[t]    <- sqrt(sigma2[t]) * z[t]
  }
  tail(eps, n)
}

simulate_garch_pq <- function(n, omega, alpha, beta, burn = 500) {
  # GARCH(p,q): sigma_t^2 = omega + sum alpha_i * eps^2_{t-i} + sum beta_j * sigma^2_{t-j}
  stopifnot(omega > 0, all(alpha >= 0), all(beta >= 0),
            sum(alpha) + sum(beta) < 1)
  q <- length(alpha); p <- length(beta)
  z <- rnorm(n + burn)
  eps <- sigma2 <- numeric(n + burn)
  sigma2[1] <- omega / (1 - sum(alpha) - sum(beta))
  start <- max(2, p + 1, q + 1)
  for (t in start:(n + burn)) {
    arch_term  <- if (q > 0) sum(alpha * eps[(t - 1):(t - q)]^2) else 0
    garch_term <- if (p > 0) sum(beta  * sigma2[(t - 1):(t - p)]) else 0
    sigma2[t] <- omega + arch_term + garch_term
    eps[t]    <- sqrt(sigma2[t]) * z[t]
  }
  tail(eps, n)
}

##### =========================================
#####  AR(1): three series (series ACF/PACF)
##### =========================================
phi_vals <- c(-0.6, 0.4, 0.8)   # stationary: |phi|<1
ar_series <- lapply(phi_vals, function(ph) arima.sim(model = list(ar = ph), n = n))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(ar_series)) {
  triptych(ar_series[[k]], sprintf("AR(1): φ=%.1f", phi_vals[k]), lags = lags, on_squared = FALSE)
}
par(op)

##### =========================================
#####  MA(1): three series (series ACF/PACF)
##### =========================================
theta_vals <- c(-0.6, 0.3, 0.8)  # invertible region
ma_series <- lapply(theta_vals, function(th) arima.sim(model = list(ma = th), n = n))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(ma_series)) {
  triptych(ma_series[[k]], sprintf("MA(1): θ=%.1f", theta_vals[k]), lags = lags, on_squared = FALSE)
}
par(op)

##### ===================================================
#####  ARCH(1): three series (ACF/PACF of SQUARED series)
##### ===================================================
# choose ω, α with α<1
arch_params <- list(
  c(omega = 0.10, alpha = 0.50),
  c(omega = 0.05, alpha = 0.80),
  c(omega = 0.20, alpha = 0.30)
)
arch_series <- lapply(arch_params, function(p) simulate_arch1(n, p["omega"], p["alpha"]))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(arch_series)) {
  lbl <- sprintf("ARCH(1): ω=%.2f, α=%.2f", arch_params[[k]]["omega"], arch_params[[k]]["alpha"])
  triptych(arch_series[[k]], lbl, lags = lags, on_squared = TRUE)
}
par(op)

##### ===============================================================
#####  GARCH(0,1): three series (≡ ARCH(1); ACF/PACF of SQUARED series)
##### ===============================================================
g01_params <- list(
  c(omega = 0.10, alpha = 0.40),
  c(omega = 0.10, alpha = 0.60),
  c(omega = 0.10, alpha = 0.80)
)
g01_series <- lapply(g01_params, function(p)
  simulate_garch_pq(n, omega = p["omega"], alpha = c(p["alpha"]), beta = numeric(0)))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(g01_series)) {
  lbl <- sprintf("GARCH(0,1): ω=%.2f, α1=%.2f", g01_params[[k]]["omega"], g01_params[[k]]["alpha"])
  triptych(g01_series[[k]], lbl, lags = lags, on_squared = TRUE)
}
par(op)
