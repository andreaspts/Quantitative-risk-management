# --- Data & marginals -------------------------------------------------------
#install.packages(c("quantmod","copula"))  # if needed
library(quantmod)
library(copula)

getSymbols("DBK.DE", from = "2007-01-01", src = "yahoo")
px  <- Ad(DBK.DE)                  # adjusted close
vol <- Vo(DBK.DE)                  # volume

r   <- diff(log(px))               # log-returns
lv  <- diff(log(vol))              # log-volume change
dat <- na.omit(merge(r, lv))
colnames(dat) <- c("ret", "lvol")

# --- Rank -> uniforms (pseudo-observations) --------------------------------
U <- pobs(as.matrix(dat))          # transforms each margin to ~U(0,1)
colnames(U) <- c("U_ret","U_lvol")

# --- Fit copulas ------------------------------------------------------------
gcop  <- normalCopula(dim = 2)
fit_g <- fitCopula(gcop, data = U, method = "ml")

tcop0 <- tCopula(dim = 2, df = 4, df.fixed = FALSE, dispstr = "un")
fit_t <- fitCopula(tcop0, data = U, method = "ml")

coef(fit_g)                        # Gaussian copula rho
coef(fit_t)                        # t-copula rho and df
AIC(fit_g); AIC(fit_t)             # compare fit (lower is better)

# Build fitted t-copula object for further use
tcop <- tCopula(param = coef(fit_t)["rho"], df = coef(fit_t)["df"],
                dim = 2, dispstr = "un")

# --- Dependence diagnostics -------------------------------------------------
kendall_tau <- cor(dat$ret, dat$lvol, method = "kendall")
tail_dep    <- lambda(tcop)        # lower/upper tail dependence (t-copula)

kendall_tau; tail_dep

# --- Plots ------------------------------------------------------------------
op <- par(mfrow = c(1,2), mar = c(4,4,2,1))

# Scatter in original margins
plot(as.numeric(dat$ret), as.numeric(dat$lvol), pch = 16, cex = 0.5,
     xlab = "DBK log-returns", ylab = "log-volume change",
     main = "Original margins")

# Normal-scores plot (rank transform) – dependence only
plot(qnorm(U[,1]), qnorm(U[,2]), pch = 16, cex = 0.5,
     xlab = "Normal score of returns", ylab = "Normal score of log-volume",
     main = "Dependence (copula scale)")
par(op)

# Optional: simulate from fitted copula, then map back via empirical quantiles
# Usim <- rCopula(5000, tcop)
# r_sim  <- quantile(dat$ret,  probs = Usim[,1], type = 8, na.rm = TRUE)
# lv_sim <- quantile(dat$lvol, probs = Usim[,2], type = 8, na.rm = TRUE)
