set.seed(123)
n    <- 400        # length of each series
lags <- 40         # ACF/PACF max lag

# helper: plot series + ACF + PACF for one ts object
triptych <- function(x, label) {
  plot(x, type = "l", col = "black", xlab = "t", ylab = "value",
       main = paste(label, "— path"))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  acf(x, lag.max = lags, main = paste("ACF —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
  pacf(x, lag.max = lags, main = paste("PACF —", label))
  grid(col = "grey80", lty = 1, lwd = 0.8)
}

##### =========================
##### MA(1), MA(2), MA(3)
##### =========================
# Choose invertible MA coefficients
ma1 <- arima.sim(model = list(ma = 0.6),                 n = n, sd = 1)
ma2 <- arima.sim(model = list(ma = c(0.6, -0.4)),        n = n, sd = 1)
ma3 <- arima.sim(model = list(ma = c(0.6, -0.4, 0.3)),   n = n, sd = 1)

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
triptych(ma1, "MA(1): θ1=0.6")
triptych(ma2, "MA(2): θ1=0.6, θ2=-0.4")
triptych(ma3, "MA(3): θ1=0.6, θ2=-0.4, θ3=0.3")
par(op)

##### =========================
##### AR(1), AR(2), AR(3)
##### =========================
# Choose stationary AR coefficients (small magnitudes ⇒ safe)
ar1 <- arima.sim(model = list(ar = 0.6),                n = n, sd = 1)
ar2 <- arima.sim(model = list(ar = c(0.4, 0.2)),        n = n, sd = 1)
ar3 <- arima.sim(model = list(ar = c(0.4, 0.2, 0.1)),   n = n, sd = 1)

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
triptych(ar1, "AR(1): φ1=0.6")
triptych(ar2, "AR(2): φ1=0.4, φ2=0.2")
triptych(ar3, "AR(3): φ1=0.4, φ2=0.2, φ3=0.1")
par(op)

##### =========================
##### ARMA(1,1), ARMA(1,2), ARMA(2,1)
##### =========================
# Choose stationary/invertible combos
arma11 <- arima.sim(model = list(ar = 0.5,            ma = -0.4),             n = n, sd = 1)
arma12 <- arima.sim(model = list(ar = 0.5,            ma = c(-0.3, 0.2)),     n = n, sd = 1)
arma21 <- arima.sim(model = list(ar = c(0.5, -0.3),   ma = 0.2),               n = n, sd = 1)

op <- par(mfrow = c(3,3), mar = c(4,4,3,1))
triptych(arma11, "ARMA(1,1): φ1=0.5; θ1=-0.4")
triptych(arma12, "ARMA(1,2): φ1=0.5; θ1=-0.3, θ2=0.2")
triptych(arma21, "ARMA(2,1): φ1=0.5, φ2=-0.3; θ1=0.2")
par(op)
