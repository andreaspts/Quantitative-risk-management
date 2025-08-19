#for ARCH(1) and GARCH(0,1) I plot the series ACF/PACF on the raw returns (no squaring), 
#just like for AR(1) and MA(1)

set.seed(123)
n    <- 600      # length of each series
lags <- 40       # ACF/PACF max lag

# ---------- helpers ----------
triptych_simple <- function(x, label, lags = 40) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "value",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(x,  lag.max = lags, main = paste("ACF —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(x, lag.max = lags, main = paste("PACF —", label))
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

simulate_garch_pq <- function(n, omega, alpha, beta = numeric(0), burn = 500) {
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

##### =========================
#####  AR(1): three series
##### =========================
phi_vals <- c(-0.6, 0.4, 0.8)
ar_series <- lapply(phi_vals, \(ph) arima.sim(model = list(ar = ph), n = n))

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
for (k in seq_along(ar_series)) {
  triptych_simple(ar_series[[k]], sprintf("AR(1): φ=%.1f", phi_vals[k]), lags)
}
par(op)

##### =========================
#####  MA(1): three series
##### =========================
theta_vals <- c(-0.6, 0.3, 0.8)
ma_series <- lapply(theta_vals, \(th) arima.sim(model = list(ma = th), n = n))

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
for (k in seq_along(ma_series)) {
  triptych_simple(ma_series[[k]], sprintf("MA(1): θ=%.1f", theta_vals[k]), lags)
}
par(op)

##### =========================
#####  ARCH(1): three series (RAW ACF/PACF)
##### =========================
arch_params <- list(
  c(omega = 0.10, alpha = 0.50),
  c(omega = 0.05, alpha = 0.80),
  c(omega = 0.20, alpha = 0.30)
)
arch_series <- lapply(arch_params, \(p) simulate_arch1(n, p["omega"], p["alpha"]))

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
for (k in seq_along(arch_series)) {
  lbl <- sprintf("ARCH(1): ω=%.2f, α=%.2f", arch_params[[k]]["omega"], arch_params[[k]]["alpha"])
  triptych_simple(arch_series[[k]], lbl, lags)
}
par(op)

##### =========================
#####  GARCH(0,1): three series (RAW ACF/PACF)
##### =========================
g01_params <- list(
  c(omega = 0.10, alpha = 0.40),
  c(omega = 0.10, alpha = 0.60),
  c(omega = 0.10, alpha = 0.80)
)
g01_series <- lapply(g01_params, \(p) simulate_garch_pq(n, omega = p["omega"], alpha = c(p["alpha"]), beta = numeric(0)))

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
for (k in seq_along(g01_series)) {
  lbl <- sprintf("GARCH(0,1): ω=%.2f, α1=%.2f", g01_params[[k]]["omega"], g01_params[[k]]["alpha"])
  triptych_simple(g01_series[[k]], lbl, lags)
}
par(op)
