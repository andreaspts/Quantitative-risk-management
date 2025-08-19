set.seed(123)
n    <- 80      # length of each series
lags <- 40       # ACF/PACF max lag

# ---------- Simulators ----------
simulate_archq <- function(n, omega, alpha, burn = 500) {
  # ARCH(q): sigma_t^2 = omega + sum_{i=1..q} alpha_i * eps_{t-i}^2
  stopifnot(omega > 0, all(alpha >= 0), sum(alpha) < 1)
  q <- length(alpha)
  z <- rnorm(n + burn)
  eps <- numeric(n + burn)
  sigma2 <- rep(omega / (1 - sum(alpha)), n + burn)
  start <- if (q == 0) 2 else (q + 1)
  for (t in start:(n + burn)) {
    arch_term <- if (q > 0) sum(alpha * eps[(t - 1):(t - q)]^2) else 0
    sigma2[t] <- omega + arch_term
    eps[t]    <- sqrt(sigma2[t]) * z[t]
  }
  tail(eps, n)
}

simulate_garch_pq <- function(n, omega, alpha, beta, burn = 500) {
  # GARCH(p,q): sigma_t^2 = omega + sum alpha_i eps^2_{t-i} + sum beta_j sigma^2_{t-j}
  stopifnot(omega > 0, all(alpha >= 0), all(beta >= 0),
            sum(alpha) + sum(beta) < 1)
  q <- length(alpha); p <- length(beta)
  z <- rnorm(n + burn)
  eps <- numeric(n + burn)
  sigma2 <- rep(omega / (1 - sum(alpha) - sum(beta)), n + burn)
  start <- max(2, p + 1, q + 1)
  for (t in start:(n + burn)) {
    arch_term  <- if (q > 0) sum(alpha * eps[(t - 1):(t - q)]^2) else 0
    garch_term <- if (p > 0) sum(beta  * sigma2[(t - 1):(t - p)]) else 0
    sigma2[t] <- omega + arch_term + garch_term
    eps[t]    <- sqrt(sigma2[t]) * z[t]
  }
  tail(eps, n)
}

# ---------- Plot helper: series + ACF|PACF of squared ----------
triptych_vol <- function(x, label, lags = 40) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "return",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(x^2, lag.max = lags, main = paste("ACF of squared —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(x^2, lag.max = lags, main = paste("PACF of squared —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
}

##### =========================================================
#####  ARCH(1), ARCH(2), ARCH(3)
##### =========================================================
arch1 <- simulate_archq(n, omega = 0.10, alpha = c(0.50))
arch2 <- simulate_archq(n, omega = 0.10, alpha = c(0.40, 0.30))
arch3 <- simulate_archq(n, omega = 0.10, alpha = c(0.30, 0.25, 0.20))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
triptych_vol(arch1, "ARCH(1): ω=0.10, α1=0.50", lags)
triptych_vol(arch2, "ARCH(2): ω=0.10, α=(0.40,0.30)", lags)
triptych_vol(arch3, "ARCH(3): ω=0.10, α=(0.30,0.25,0.20)", lags)
par(op)

##### =========================================================
#####  GARCH(0,1), GARCH(0,2), GARCH(0,3)   (≡ ARCH(1/2/3))
##### =========================================================
g01 <- simulate_garch_pq(n, omega = 0.10, alpha = c(0.50),                  beta = numeric(0))
g02 <- simulate_garch_pq(n, omega = 0.10, alpha = c(0.40, 0.30),            beta = numeric(0))
g03 <- simulate_garch_pq(n, omega = 0.10, alpha = c(0.30, 0.25, 0.20),      beta = numeric(0))

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
triptych_vol(g01, "GARCH(0,1): ω=0.10, α1=0.50", lags)
triptych_vol(g02, "GARCH(0,2): ω=0.10, α=(0.40,0.30)", lags)
triptych_vol(g03, "GARCH(0,3): ω=0.10, α=(0.30,0.25,0.20)", lags)
par(op)

##### =========================================================
#####  GARCH(1,1), GARCH(2,1), GARCH(1,2)
##### =========================================================
g11 <- simulate_garch_pq(n, omega = 0.05, alpha = c(0.10),        beta = c(0.85))          # α+β=0.95
g21 <- simulate_garch_pq(n, omega = 0.05, alpha = c(0.08, 0.04),  beta = c(0.85))          # sum=0.97
g12 <- simulate_garch_pq(n, omega = 0.05, alpha = c(0.08),        beta = c(0.70, 0.20))    # sum=0.98

op <- par(mfrow = c(3, 3), mar = c(4,4,3,1))
triptych_vol(g11, "GARCH(1,1): ω=0.05, α=(0.10), β=(0.85)", lags)
triptych_vol(g21, "GARCH(2,1): ω=0.05, α=(0.08,0.04), β=(0.85)", lags)
triptych_vol(g12, "GARCH(1,2): ω=0.05, α=(0.08), β=(0.70,0.20)", lags)
par(op)
