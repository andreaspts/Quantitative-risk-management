# MC application: binomial “up/down” stock over 12 months
set.seed(1)

# --- Inputs ---
S0 <- 100          # starting price
N  <- 12           # months
p  <- 0.2          # prob of an "up" month
u  <- 1.02         # +3% month
d  <- 0.98         # -3% month
n  <- 1000000        # Monte Carlo paths

# --- Simulate ---
k  <- rbinom(n, size = N, prob = p)        # number of up months
ST <- S0 * (u^k) * (d^(N - k))             # terminal price per path

# --- Simple stats you care about ---
mean_ST   <- mean(ST)                       # expected terminal price
prob_loss <- mean(ST < S0)                  # chance you end below S0
ci_mean   <- mean_ST + c(-1,1) * 1.96 * sd(ST)/sqrt(n)  # 95% CI for E[ST]

cat(sprintf("E[ST] ≈ %.2f (95%% CI [%.2f, %.2f])\n", mean_ST, ci_mean[1], ci_mean[2]))
cat(sprintf("P(ST < %g) ≈ %.3f\n", S0, prob_loss))

# --- Quick visualization ---
hist(ST, breaks = 50, probability = TRUE,
     main = "1-year terminal price (binomial up/down)",
     xlab = "ST")
abline(v = S0, col = "red", lwd = 2, lty = 2)  # starting price
