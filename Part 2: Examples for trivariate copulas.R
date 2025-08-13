## Trivariate copulas: simulate + pairwise + 3D plots (n = 3000)
#install.packages(c("copula","scatterplot3d"))  # uncomment if needed
library(copula)
library(scatterplot3d)

set.seed(42)
n <- 3000

## ---- Define 3D copulas -------------------------------------------------
indep  <- indepCopula(dim = 3)

# For Gaussian/t, unconstrained correlation structure ("un"):
# param = (rho12, rho13, rho23); must form a valid PD correlation matrix.
rho_vec <- c(0.7, 0.5, 0.6)  # ρ12=0.7, ρ13=0.5, ρ23=0.6

gauss  <- normalCopula(param = rho_vec, dim = 3, dispstr = "un")
tcop   <- tCopula(param = rho_vec, df = 4, dim = 3, dispstr = "un")  # heavier tails

# Archimedean copulas
clay   <- claytonCopula(param = 2,  dim = 3)   # lower-tail dependence
frank  <- frankCopula(param = 10,  dim = 3)    # no tail dependence
gumbel <- gumbelCopula(param = 2,  dim = 3)    # upper-tail dependence (θ ≥ 1)

## ---- Simulate U ~ C(u1,u2,u3) -----------------------------------------
U_ind <- rCopula(n, indep)
U_g   <- rCopula(n, gauss)
U_t   <- rCopula(n, tcop)
U_c   <- rCopula(n, clay)
U_f   <- rCopula(n, frank)
U_gu  <- rCopula(n, gumbel)

## ---- Helper: pairwise scatter matrix on [0,1]^3 -----------------------
plotPairs <- function(U, title) {
  colnames(U) <- c("u1","u2","u3")
  pairs(U,
        pch = 16, cex = 0.35, col = adjustcolor("black", 0.5),
        labels = c(expression(u[1]), expression(u[2]), expression(u[3])),
        main = title)
}

## ---- Helper: static 3D scatter using scatterplot3d --------------------
plot3D <- function(U, title, angle = 55) {
  scatterplot3d(U[,1], U[,2], U[,3],
                pch = 16, cex.symbols = 0.35,
                color = adjustcolor("black", 0.5),
                xlab = "u1", ylab = "u2", zlab = "u3",
                main = title, box = TRUE, angle = angle)
}

## ---- 2 x 3 grid: pairwise matrices ------------------------------------
op <- par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2, 1))
plotPairs(U_ind, "Independence")
plotPairs(U_g,   "Gaussian (ρ12=0.7, ρ13=0.5, ρ23=0.6)")
plotPairs(U_t,   "t (df = 4; ρ as above)")
plotPairs(U_c,   "Clayton (θ = 2)")
plotPairs(U_f,   "Frank (θ = 10)")
plotPairs(U_gu,  "Gumbel (θ = 2)")
par(op)

## ---- 2 x 3 grid: static 3D plots --------------------------------------
op <- par(mfrow = c(2, 3), mar = c(3.5, 3.5, 2, 1))
plot3D(U_ind, "Independence")
plot3D(U_g,   "Gaussian")
plot3D(U_t,   "t (df = 4)")
plot3D(U_c,   "Clayton (θ = 2)")
plot3D(U_f,   "Frank (θ = 10)")
plot3D(U_gu,  "Gumbel (θ = 2)")
par(op)

## ---- OPTIONAL: interactive 3D (rgl) -----------------------------------
## install.packages("rgl")  # if you want interactive, uncomment and run
# if (requireNamespace("rgl", quietly = TRUE)) {
#   library(rgl)
#   mfrow3d(2, 3, sharedMouse = TRUE)
#   plotRgl <- function(U, title) {
#     next3d(); bg3d("white")
#     points3d(U[,1], U[,2], U[,3], size = 2)
#     axes3d(); title3d(main = title, xlab = "u1", ylab = "u2", zlab = "u3")
#   }
#   plotRgl(U_ind, "Independence")
#   plotRgl(U_g,   "Gaussian")
#   plotRgl(U_t,   "t (df = 4)")
#   plotRgl(U_c,   "Clayton (θ = 2)")
#   plotRgl(U_f,   "Frank (θ = 10)")
#   plotRgl(U_gu,  "Gumbel (θ = 2)")
# }
