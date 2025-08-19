set.seed(123)
n    <- 60        # length of kept series
burn <- 800        # burn-in to reach stationarity
lags <- 40         # max lag for ACF/PACF

# -------- simulator: ARMA(1,1)–GARCH(1,1) -------------------------------
# Model:
#   (x_t - mu) = phi * (x_{t-1} - mu) + eps_t + theta * eps_{t-1}
#   eps_t = sigma_t * z_t,  z_t ~ N(0,1)
#   sigma_t^2 = omega + alpha * eps_{t-1}^2 + beta * sigma_{t-1}^2
simulate_arma11_garch11 <- function(n, burn, mu, phi, theta, omega, alpha, beta) {
  stopifnot(abs(phi) < 1, omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1)
  N <- n + burn
  z <- rnorm(N)
  x <- eps <- sigma2 <- numeric(N)
  
  # start at unconditional variance
  sigma2[1] <- omega / (1 - alpha - beta)
  eps[1]    <- sqrt(sigma2[1]) * z[1]
  x[1]      <- mu + eps[1]                       # no lags at t=1
  
  for (t in 2:N) {
    sigma2[t] <- omega + alpha * eps[t-1]^2 + beta * sigma2[t-1]
    eps[t]    <- sqrt(sigma2[t]) * z[t]
    x[t]      <- mu + phi * (x[t-1] - mu) + eps[t] + theta * eps[t-1]
  }
  tail(x, n)
}

# -------- three parameter sets (all stationary/invertible) -----------------
params <- list(
  list(mu = 0, phi =  0.30, theta = -0.40, omega = 0.05, alpha = 0.05, beta = 0.90),  # α+β=0.95
  list(mu = 0, phi = -0.50, theta =  0.50, omega = 0.05, alpha = 0.10, beta = 0.85),  # α+β=0.95
  list(mu = 0, phi =  0.70, theta = -0.20, omega = 0.05, alpha = 0.08, beta = 0.88)   # α+β=0.96
)

series <- lapply(params, function(p) do.call(simulate_arma11_garch11, c(list(n = n, burn = burn), p)))

# -------- helper: plot path + ACF + PACF for one series --------------------
triptych <- function(x, label, lags = 40) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "value",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(x, lag.max = lags, main = paste("ACF —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(x, lag.max = lags, main = paste("PACF —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
}

# -------- plot: 3×3 grid (Series | ACF | PACF) for the three models -------
labels <- sapply(params, function(p)
  sprintf("ARMA(1,1)–GARCH(1,1): μ=%.1f, φ=%.2f, θ=%.2f, ω=%.2f, α=%.2f, β=%.2f",
          p$mu, p$phi, p$theta, p$omega, p$alpha, p$beta))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
for (k in 1:3) triptych(series[[k]], labels[k], lags = lags)
par(op)

# (Optional) If you also want to see volatility clustering explicitly:
# op <- par(mfrow = c(3, 2), mar = c(4,4,3,1))
# for (k in 1:3) {
#   plot(series[[k]], type="l", main=paste(labels[k], "— path"), xlab="t", ylab="value"); grid(col="grey80", lty=1)
#   acf(series[[k]]^2, lag.max=lags, main=paste("ACF of squared — model", k)); grid(col="grey80", lty=1)
# }
# par(op)
