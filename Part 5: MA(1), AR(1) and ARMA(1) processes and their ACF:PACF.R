set.seed(123)
n    <- 50
cols <- c("firebrick", "dodgerblue3", "darkgreen")

# ---------- MA(1) ----------
theta <- c(-0.8, 0.2, 0.8)
ma_list <- lapply(theta, function(th)
  arima.sim(model = list(ma = th), n = n, sd = 1))
MA <- do.call(cbind, ma_list); colnames(MA) <- paste0("theta=", theta)

par(mfrow = c(1,3), mar = c(4,4,3,1))
matplot(MA, type = "l", lty = 1, col = cols,
        xlab = "t", ylab = "value", main = "MA(1) paths")
grid(col = "grey80", lty = 1)
legend("topleft", bty = "n", lty = 1, col = cols, legend = colnames(MA))
acf(MA[,1], lag.max = 40, main = paste("ACF: MA(1),", colnames(MA)[1]))
grid(col = "grey80", lty = 1)
pacf(MA[,1], lag.max = 40, main = paste("PACF: MA(1),", colnames(MA)[1]))
grid(col = "grey80", lty = 1)

# ---------- AR(1) ----------
phi <- c(-0.7, 0.4, 0.8)
ar_list <- lapply(phi, function(ph)
  arima.sim(model = list(ar = ph), n = n, sd = 1))
AR <- do.call(cbind, ar_list); colnames(AR) <- paste0("phi=", phi)

par(mfrow = c(1,3), mar = c(4,4,3,1))
matplot(AR, type = "l", lty = 1, col = cols,
        xlab = "t", ylab = "value", main = "AR(1) paths")
grid(col = "grey80", lty = 1)
legend("topleft", bty = "n", lty = 1, col = cols, legend = colnames(AR))
acf(AR[,1], lag.max = 40, main = paste("ACF: AR(1),", colnames(AR)[1]))
grid(col = "grey80", lty = 1)
pacf(AR[,1], lag.max = 40, main = paste("PACF: AR(1),", colnames(AR)[1]))
grid(col = "grey80", lty = 1)

# ---------- ARMA(1,1) ----------
arma_params <- list(c(phi = 0.6,  theta = -0.5),
                    c(phi = -0.4, theta =  0.7),
                    c(phi = 0.8,  theta =  0.3))
arma_list <- lapply(arma_params, function(p)
  arima.sim(model = list(ar = p["phi"], ma = p["theta"]), n = n, sd = 1))
ARMA <- do.call(cbind, arma_list)
colnames(ARMA) <- sapply(arma_params, function(p)
  sprintf("phi=%.1f, theta=%.1f", p["phi"], p["theta"]))

par(mfrow = c(1,3), mar = c(4,4,3,1))
matplot(ARMA, type = "l", lty = 1, col = cols,
        xlab = "t", ylab = "value", main = "ARMA(1,1) paths")
grid(col = "grey80", lty = 1)
legend("topleft", bty = "n", lty = 1, col = cols, legend = colnames(ARMA))
acf(ARMA[,1], lag.max = 40, main = paste("ACF: ARMA(1,1),", colnames(ARMA)[1]))
grid(col = "grey80", lty = 1)
pacf(ARMA[,1], lag.max = 40, main = paste("PACF: ARMA(1,1),", colnames(ARMA)[1]))
grid(col = "grey80", lty = 1)

par(mfrow = c(1,1))
