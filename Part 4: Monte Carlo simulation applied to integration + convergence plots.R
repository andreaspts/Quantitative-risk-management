set.seed(1)

## Target integral
f <- function(x) exp(-x^2)
a <- 0; b <- 1
I_true <- integrate(f, lower = a, upper = b)$value   # benchmark

## ---- Naive Monte Carlo -------------------------------------------------
n  <- 100
u  <- runif(n, a, b)
fx <- f(u)

k        <- 1:n
cs       <- cumsum(fx)
cs2      <- cumsum(fx^2)

Ihat     <- (b - a) * cs / k
var_fx   <- (cs2 - cs^2 / k) / pmax(k - 1, 1)        # running unbiased var
se       <- (b - a) * sqrt(var_fx / k)
lo       <- Ihat - 1.96 * se
hi       <- Ihat + 1.96 * se

## ---- Antithetic variates (variance reduction) --------------------------
m   <- ceiling(n / 2)                 # number of pairs
u1  <- runif(m)
u2  <- 1 - u1                         # antithetic mate
fxA <- (f(u1) + f(u2)) / 2            # paired average

kA   <- 1:m
csA  <- cumsum(fxA)
cs2A <- cumsum(fxA^2)

IhatA <- (b - a) * csA / kA
varA  <- (cs2A - csA^2 / kA) / pmax(kA - 1, 1)
seA   <- (b - a) * sqrt(varA / kA)
loA   <- IhatA - 1.96 * seA
hiA   <- IhatA + 1.96 * seA

## ---- Plots --------------------------------------------------------------
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# (A) Naive MC
plot(k, Ihat, type = "l",
     xlab = "samples n", ylab = "estimate",
     main = expression(paste("MC convergence for  ", integral(e^{-x^2} * dx, 0, 1))))
lines(k, lo, col = "grey70")
lines(k, hi, col = "grey70")
abline(h = I_true, col = "red", lwd = 2)
legend("topright", bty = "n",
       legend = c("MC estimate", "95% CI", "true value"),
       col = c("black", "grey70", "red"), lty = 1, lwd = c(1, 1, 2))

# (B) Antithetic variates (x-axis in total uniforms ≈ 2*kA)
plot(2 * kA, IhatA, type = "l",
     xlab = "samples n (counting both antithetic uniforms)",
     ylab = "estimate",
     main = "Antithetic variates")
lines(2 * kA, loA, col = "grey70")
lines(2 * kA, hiA, col = "grey70")
abline(h = I_true, col = "red", lwd = 2)
legend("topright", bty = "n",
       legend = c("Antithetic estimate", "95% CI", "true value"),
       col = c("black", "grey70", "red"), lty = 1, lwd = c(1, 1, 2))

par(op)

I_true    # for reference
tail(cbind(n = k, Ihat, se, lo, hi), 1)     # final estimate & CI (naive)
tail(cbind(n = 2*kA, IhatA, seA, loA, hiA), 1) # final estimate & CI (antithetic)
