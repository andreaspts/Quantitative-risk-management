# Six bivariate copulas: simulate and plot (n = 1000)
# install.packages("copula")  # uncomment if needed
library(copula)
set.seed(42)
n <- 1000

# --- Define copulas (choose parameters that show typical behavior) ----------
indep  <- indepCopula(dim = 2)

gauss  <- normalCopula(param = 0.7, dim = 2)       # rho = 0.7
tcop   <- tCopula(param = 0.7, df = 4, dim = 2)    # rho = 0.7, df = 4 (tail dep)

clay   <- claytonCopula(param = 2, dim = 2)        # lower-tail dependence
frank  <- frankCopula(param = 10, dim = 2)         # no tail dep, flexible center
gumbel <- gumbelCopula(param = 2, dim = 2)         # upper-tail dependence

# --- Simulate uniforms U ~ C(u1,u2) -----------------------------------------
U_ind <- rCopula(n, indep)
U_g   <- rCopula(n, gauss)
U_t   <- rCopula(n, tcop)
U_c   <- rCopula(n, clay)
U_f   <- rCopula(n, frank)
U_gu  <- rCopula(n, gumbel)

# --- Plot helper (unit square scatter) --------------------------------------
plotU <- function(U, title) {
  plot(U[,1], U[,2],
       pch = 16, cex = 0.6, col = adjustcolor("black", 0.5),
       xlab = expression(u[1]), ylab = expression(u[2]),
       xlim = c(0,1), ylim = c(0,1), main = title)
  rect(0, 0, 1, 1, border = "grey80")
}

# --- 2 x 3 panel of the six copulas ----------------------------------------
op <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plotU(U_ind, "Independence")
plotU(U_g,   "Gaussian (rho = 0.7)")
plotU(U_t,   "t (rho = 0.7, df = 4)")
plotU(U_c,   "Clayton (theta = 2)")
plotU(U_f,   "Frank (theta = 10)")
plotU(U_gu,  "Gumbel (theta = 2)")
par(op)

# --- Optional: view the same samples with standard normal margins -----------
# (purely cosmetic; dependence/coupling is unchanged)
# Z <- function(U) qnorm(U)
# op <- par(mfrow = c(2, 3), mar = c(4,4,2,1))
# plotU(Z(U_ind), "Independence (normal margins)")
# plotU(Z(U_g),   "Gaussian (normal margins)")
# plotU(Z(U_t),   "t (normal margins)")
# plotU(Z(U_c),   "Clayton (normal margins)")
# plotU(Z(U_f),   "Frank (normal margins)")
# plotU(Z(U_gu),  "Gumbel (normal margins)")
# par(op)
