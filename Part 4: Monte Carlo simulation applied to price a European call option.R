set.seed(1)

# --- Inputs ---
S0    <- 100      # spot
K     <- 100      # strike
r     <- 0.02     # risk-free rate
sigma <- 0.20     # vol
T     <- 1        # years
n     <- 1e5      # simulations

# --- Antithetic variates for lower variance ---
m  <- n %/% 2
Z  <- rnorm(m)
Z  <- c(Z, -Z)    # antithetic pair
ST <- S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * Z)

payoff <- pmax(ST - K, 0)
disc   <- exp(-r * T)
price_hat <- disc * mean(payoff)
se        <- disc * sd(payoff) / sqrt(n)
ci95      <- price_hat + c(-1, 1) * 1.96 * se

# --- Closed-form Black–Scholes for a check ---
bs_call <- function(S0, K, r, sigma, T) {
  d1 <- (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
}
price_bs <- bs_call(S0, K, r, sigma, T)

cat(sprintf("MC price: %.4f  (SE %.4f)  95%% CI [%.4f, %.4f]\n",
            price_hat, se, ci95[1], ci95[2]))
cat(sprintf("BS price: %.4f\n", price_bs))

# --- Quick visualization ---
hist(ST, breaks = 80, probability = TRUE,
     main = "Terminal price distribution ST",
     xlab = "ST")
abline(v = K, col = "red", lwd = 2, lty = 2)  # strike
