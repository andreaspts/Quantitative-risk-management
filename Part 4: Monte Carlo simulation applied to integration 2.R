# another example for MC integration: ∫_0^1 (x^3 - 0.5 x^2) dx

f <- function(x) {x^3-0.5*x^2} #Defining a function for integration

I<-rep (0,200) #Initializing vector of calculated integrals 
#for n=1,...，200

set.seed(7)
for(j in 1:length(I)) {
    x=runif(j,0,1) #simulate data from uniform distribution
    for (k in 1:j) {
      I[j]=I[j]+f(x[k])
    }
  I[j]=I[j]/j    #calculate arithmetic mean
}
  
#plotting results
plot((1:length(I)),I,type="l",xlab="number of observations",
       ylab="calculated integral")
abline (h=1/12, col="blue")


#improved version thereof
f <- function(x) x^3 - 0.5 * x^2          # define integrand

I <- numeric(200)                         # store estimates for n = 1..200

set.seed(7)
for (j in seq_along(I)) {
  x <- runif(j, 0, 1)                     # sample U(0,1)
  I[j] <- mean(f(x))                       # Monte Carlo estimate (average)
}

# plot
plot(seq_along(I), I, type = "l",
     xlab = "number of observations", ylab = "calculated integral")
abline(h = 1/12, col = "blue")            # true value

#higher dimensional scenario: 3d case
set.seed(7)

# Integrand on [0,1]^3
f <- function(x, y, z) x^3 * y^2 + 2 * x * y * z + z^4

# ---- Monte Carlo sampling ----
N <- 100000
U <- matrix(runif(3 * N), ncol = 3)
vals <- f(U[,1], U[,2], U[,3])         # f(X,Y,Z)

# Running estimate + 95% CI
k    <- 1:N
cs   <- cumsum(vals)
cs2  <- cumsum(vals^2)
Ihat <- cs / k                          # volume = 1 on [0,1]^3
var_hat <- (cs2 - cs^2 / k) / pmax(k - 1, 1)
se   <- sqrt(var_hat / k)
lo   <- Ihat - 1.96 * se
hi   <- Ihat + 1.96 * se

I_true <- 1/12 + 1/4 + 1/5              # = 8/15 ≈ 0.5333333

# ---- Plots ----
op <- par(mfrow = c(1, 2), mar = c(4,4,3,1))

## (A) Convergence plot
plot(k, Ihat, type = "l",
     xlab = "samples n", ylab = "estimate",
     main = expression(paste("MC convergence for  ",
                             tripleintegral(f(x,y,z) * dx * dy * dz, 0, 1, 0, 1, 0, 1))))
lines(k, lo, col = "grey70")
lines(k, hi, col = "grey70")
abline(h = I_true, col = "red", lwd = 2)
legend("topright", bty = "n",
       legend = c("estimate", "95% CI", "true value (= 8/15)"),
       col = c("black", "grey70", "red"), lty = 1, lwd = c(1,1,2))

## (B) 3D scatter colored by f(x,y,z)
# install.packages("scatterplot3d")  # uncomment if needed
library(scatterplot3d)
n_plot <- min(N, 3000)                 # keep the 3D plot light
idx <- sample.int(N, n_plot)
vals_plot <- vals[idx]
U_plot <- U[idx, , drop = FALSE]

# simple color ramp by function value
vals_norm <- (vals_plot - min(vals_plot)) / (max(vals_plot) - min(vals_plot) + 1e-12)
pal <- colorRampPalette(c("navy","skyblue","orange","red"))(100)
cols <- pal[pmax(1, pmin(100, floor(vals_norm * 99) + 1))]

scatterplot3d(U_plot[,1], U_plot[,2], U_plot[,3],
              pch = 16, cex.symbols = 0.5, color = cols,
              xlab = "x", ylab = "y", zlab = "z",
              main = expression(paste("Samples on [0,1]^3 colored by  f(x,y,z)")))

par(op)

# Quick printout
cat(sprintf("True integral: %.6f\n", I_true))
cat(sprintf("MC (N=%d):     %.6f  (SE %.6f)\n", N, Ihat[N], se[N]))
