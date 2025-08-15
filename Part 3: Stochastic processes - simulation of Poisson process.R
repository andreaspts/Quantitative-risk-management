# ==========================================================
# Simulate paths of a Poisson counting process
# Inputs:
#   lambda  > 0   intensity
#   ymax         number of jumps per path (levels = S0, ..., S0 + ymax)
#   S0            starting count level (integer, typically 0)
#   NRepl         number of simulated paths
#   supr          0 = no plot, nonzero = plot
# ==========================================================

poissonsimulation <- function(lambda, ymax, S0 = 0, NRepl = 1, supr = 1) {
  stopifnot(lambda > 0, ymax >= 1, NRepl >= 1)
  
  # Jump times matrix: each row j holds the (ymax) ordered jump times
  x <- matrix(0, nrow = NRepl, ncol = ymax)
  x[, 1] <- 0  # first jump time measured from 0
  
  # y-levels of the counting process (jumps of height 1)
  y <- S0:(S0 + ymax - 1)
  
  if (ymax > 1) {
    for (j in seq_len(NRepl)) {
      for (i in 1:(ymax - 1)) {
        x[j, i + 1] <- x[j, i] + rexp(1, rate = lambda)  # exponential waiting times
      }
    }
  }
  
  # Plot (step functions of a Poisson counting process)
  if (supr != 0) {
    plot(NA,
         xlim = c(0, max(x)),
         ylim = c(S0, S0 + ymax),
         xlab = "t",
         ylab = "N(t)",
         main = sprintf("Poisson process with intensity λ = %.3g", lambda))
    for (j in seq_len(NRepl)) {
      lines(stats::stepfun(x[j, ], c(y, S0 + ymax)), col = "blue")
    }
  }
  
  invisible(list(times = x, levels = y))
}

# Example:
# set.seed(1)
# poissonsimulation(lambda = 2, ymax = 10, S0 = 0, NRepl = 5, supr = 1)
poissonsimulation(lambda = 2, ymax = 12, S0 = 0, NRepl = 3, supr = 1)
