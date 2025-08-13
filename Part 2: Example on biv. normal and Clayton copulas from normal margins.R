## Normal marginals + two copulas (Gaussian rho=0.7, Clayton theta=3)
## install.packages("copula")          # uncomment if needed
## install.packages("MASS")            # optional, for contour overlays
library(copula)

set.seed(123)
n <- 1000

## Marginal parameters (edit if you want non-standard normals)
mu    <- c(0, 0)
sigma <- c(1, 1)

## 1) Gaussian copula with correlation 0.7
gcop <- normalCopula(param = 0.7, dim = 2)
U_g  <- rCopula(n, gcop)  # uniforms on [0,1]^2 from Gaussian copula
X_g  <- cbind(qnorm(U_g[,1], mean = mu[1], sd = sigma[1]),
              qnorm(U_g[,2], mean = mu[2], sd = sigma[2]))  # Normal margins

## 2) Clayton copula with theta = 3 (lower-tail dependence)
ccop <- claytonCopula(param = 3, dim = 2)
U_c  <- rCopula(n, ccop)
X_c  <- cbind(qnorm(U_c[,1], mean = mu[1], sd = sigma[1]),
              qnorm(U_c[,2], mean = mu[2], sd = sigma[2]))

## ---- Plot helpers -------------------------------------------------------
plot_joint <- function(X, title) {
  plot(X[,1], X[,2],
       pch = 16, cex = 0.6, col = adjustcolor("black", 0.5),
       xlab = "X1 (Normal margin)", ylab = "X2 (Normal margin)", main = title)
  if (requireNamespace("MASS", quietly = TRUE)) {
    kd <- MASS::kde2d(X[,1], X[,2], n = 60)
    contour(kd, add = TRUE, drawlabels = FALSE)
  }
}

## ---- Make the plots -----------------------------------------------------
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot_joint(X_g, "Gaussian copula (ρ = 0.7) + Normal margins")
plot_joint(X_c, "Clayton copula (θ = 3) + Normal margins")
par(op)

## ---- Optional: print dependence measures --------------------------------
cat("Pearson r (Gaussian):", cor(X_g[,1], X_g[,2]), "\n")
cat("Pearson r (Clayton) :", cor(X_c[,1], X_c[,2]), "\n")
cat("Kendall tau (Gauss) :", cor(U_g[,1], U_g[,2], method = "kendall"), "\n")
cat("Kendall tau (Clayton):", cor(U_c[,1], U_c[,2], method = "kendall"),
    " (theory = theta/(theta+2) = ", 3/(3+2), ")\n", sep = "")
