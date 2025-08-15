#============================================================
# Function simulates paths of jump diffusion processes.=====
#============================================================
# Inputs:
#   S0      starting value
#   mu      drift
#   sigma   volatility
#   lambda  jump intensity (per unit time)
#   T       time horizon
#   NSteps  number of time steps
#   NRepl   number of paths
#   supr    0 to suppress plot
#============================================================

jumpdiffusionPaths <- function(S0, mu, sigma, lambda, T,
                               NSteps, NRepl, supr) {
  # Initialize paths
  JPaths <- matrix(0, nrow = NRepl, ncol = 1 + NSteps)
  JPaths[, 1] <- S0
  
  # Discretization
  dt   <- T / NSteps
  nudt <- (mu - 0.5 * sigma^2) * dt        # drift term per step (GBM in log)
  sidt <- sigma * sqrt(dt)                 # diffusion term per step
  
  # --- Poisson jumps: number of jumps per step dN[i, j] ~ Poisson(lambda * dt)
  dN <- matrix(rpois(NRepl * NSteps, lambda = lambda * dt),
               nrow = NRepl, ncol = NSteps)
  
  # (Jump sizes) If a step has dN = k > 0, add sum of k i.i.d. N(0,1) jumps
  # (Keeps your interface; change rnorm(...) to control jump distribution)
  # Example alternative (Merton): rnorm(k, mean = muJ, sd = sigmaJ)
  
  # Simulate the paths
  for (i in 1:NRepl) {           # NRepl paths
    for (j in 1:NSteps) {        # NSteps points
      k  <- dN[i, j]                          # number of jumps this step
      J  <- if (k > 0) sum(rnorm(k)) else 0   # aggregate jump size
      Z  <- rnorm(1)                          # diffusion shock
      JPaths[i, j + 1] <- JPaths[i, j] * exp(nudt + sidt * Z + J)
    }
  }
  
  # Graphic output
  if (supr != 0) {
    # empty canvas close to your original plotting style
    plot(0:NSteps,
         seq(min(JPaths), max(JPaths), length.out = NSteps + 1),
         type = "n", xlab = "step", ylab = "", main = "Jump Diffusion Process")
    for (j in 1:nrow(JPaths)) {
      lines(0:NSteps, JPaths[j, ], col = j)
    }
  }
  
  return(JPaths)
}

# --- Example run ---
# set.seed(1)
# JP <- jumpdiffusionPaths(S0 = 100, mu = 0.05, sigma = 0.2,
#                          lambda = 0.3, T = 1, NSteps = 252, NRepl = 5, supr = 1)

JP <- jumpdiffusionPaths(S0 = 50, mu = 0.1, sigma = 0.3,
                         lambda = 0.03, T = 1, NSteps = 365, NRepl = 3, supr = 1)