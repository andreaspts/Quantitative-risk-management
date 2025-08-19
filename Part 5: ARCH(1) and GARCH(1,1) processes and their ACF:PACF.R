set.seed(123)
n    <- 800          # length of each series
lags <- 40           # ACF/PACF max lag

# ---------- simulators (no extra packages) ----------
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

simulate_garch11 <- function(n, omega, alpha, beta, burn = 500) {
  stopifnot(omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1)
  z <- rnorm(n + burn)
  eps <- sigma2 <- numeric(n + burn)
  sigma2[1] <- omega / (1 - alpha - beta)
  for (t in 2:(n + burn)) {
    sigma2[t] <- omega + alpha * eps[t - 1]^2 + beta * sigma2[t - 1]
    eps[t]    <- sqrt(sigma2[t]) * z[t]
  }
  tail(eps, n)
}

# ---------- plotting helper (series + ACF|PACF of squared) ----------
triptych_vol <- function(x, label, lags = 40) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "return",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(x^2, lag.max = lags, main = paste("ACF of squared —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(x^2, lag.max = lags, main = paste("PACF of squared —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
}

##### =========================
#####  ARCH(1): three settings
##### =========================
arch_params <- list(
  c(omega = 0.10, alpha = 0.50),
  c(omega = 0.05, alpha = 0.80),
  c(omega = 0.20, alpha = 0.30)
)

arch_series <- lapply(arch_params, function(p)
  simulate_arch1(n, omega = p["omega"], alpha = p["alpha"]))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(arch_series)) {
  lbl <- sprintf("ARCH(1): ω=%.2f, α=%.2f", arch_params[[k]]["omega"], arch_params[[k]]["alpha"])
  triptych_vol(arch_series[[k]], lbl, lags = lags)
}
par(op)

##### =========================
#####  GARCH(1,1): three settings
##### =========================
garch_params <- list(
  c(omega = 0.05, alpha = 0.05, beta = 0.90),
  c(omega = 0.05, alpha = 0.10, beta = 0.85),
  c(omega = 0.05, alpha = 0.20, beta = 0.75)
)
# note: alpha + beta < 1 for covariance-stationarity

garch_series <- lapply(garch_params, function(p)
  simulate_garch11(n, omega = p["omega"], alpha = p["alpha"], beta = p["beta"]))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in seq_along(garch_series)) {
  lbl <- sprintf("GARCH(1,1): ω=%.2f, α=%.2f, β=%.2f",
                 garch_params[[k]]["omega"], garch_params[[k]]["alpha"], garch_params[[k]]["beta"])
  triptych_vol(garch_series[[k]], lbl, lags = lags)
}
par(op)
