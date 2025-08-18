#another example for using MC simulation for integration, now in 10d
#note that work for MC grows linearly with N , regardless of dimension; 
#the error shrinks like O(N^−1/2).
#Work for tensor grids grows as m^d (exponential in d); even modest 
#m explodes in high d.
#The integrand is not separable, 
#so product rules do not simplify it;
#MC stays simple and accurate.

set.seed(7)

# ----- Problem setup ---------------------------------------------------------
d <- 10
f_sum3 <- function(U) {                    # U is an N x d matrix in [0,1]
  S <- rowSums(U)
  S^3
}

I_true <- d/4 + d*(d-1)/2 + d*(d-1)*(d-2)/8   # = 137.5 for d = 10

# ----- Monte Carlo -----------------------------------------------------------
N <- 200000
U <- matrix(runif(N * d), ncol = d)
vals <- f_sum3(U)

Ihat <- mean(vals)
se   <- sd(vals) / sqrt(N)
ci95 <- Ihat + c(-1, 1) * 1.96 * se

cat(sprintf("Monte Carlo (N=%d): Î = %.6f, SE = %.6f, 95%% CI = [%.6f, %.6f]\n",
            N, Ihat, se, ci95[1], ci95[2]))
cat(sprintf("Exact:              I  = %.6f\n", I_true))

# Convergence (running mean & 95% CI)
k   <- 1:N
cs  <- cumsum(vals)
cs2 <- cumsum(vals^2)
Irun <- cs / k
var_hat <- (cs2 - cs^2 / k) / pmax(k - 1, 1)
se_run  <- sqrt(var_hat / k)
lo <- Irun - 1.96 * se_run
hi <- Irun + 1.96 * se_run

op <- par(mfrow = c(1, 2), mar = c(4,4,3,1))

plot(k, Irun, type = "l",
     xlab = "samples n", ylab = "estimate",
     main = "MC convergence in 10D")
lines(k, lo, col = "grey70")
lines(k, hi, col = "grey70")
abline(h = I_true, col = "red", lwd = 2)
legend("topright", bty = "n",
       legend = c("estimate", "95% CI", "true value"),
       col = c("black", "grey70", "red"), lty = 1, lwd = c(1,1,2))

# ----- Tiny tensor mid-point grid (illustration) -----------------------------
m <- 3                                  # 3 points per axis -> 3^10 = 59,049 evals
mid <- (seq_len(m) - 0.5) / m
# Build all combinations (careful: grows as m^d)
G <- expand.grid(rep(list(mid), d))
vals_grid <- f_sum3(as.matrix(G))
I_grid <- mean(vals_grid)               # mid-point rule on [0,1]^d (volume = 1)

plot(density(vals), main = "Distribution of f(U) under MC",
     xlab = expression((sum(U[i]))^3)); abline(v = I_true, col = "red", lwd = 2)
par(op)

cat(sprintf("Tensor mid-point (m=%d per axis, %d pts): I_grid ≈ %.6f\n",
            m, m^d, I_grid))

# Note: increasing m quickly becomes infeasible: m=5 -> 9,765,625 points; m=7 -> 282,475,249.
