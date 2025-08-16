#Metripolis MC applied to 2d Ising model of statistical physics
set.seed(1)

# ----- helpers --------------------------------------------------------------
idxp <- function(i, L) if (i < L) i + 1 else 1
idxm <- function(i, L) if (i > 1) i - 1 else L

local_dE <- function(S, i, j, L) {
  s  <- S[i, j]
  nb <- S[idxp(i,L), j] + S[idxm(i,L), j] + S[i, idxp(j,L)] + S[i, idxm(j,L)]
  2 * s * nb  # J=1
}

energy_per_spin <- function(S) {
  L <- nrow(S)
  up    <- rbind(S[L, ], S[1:(L-1), ])
  right <- cbind(S[, 2:L], S[, 1])
  E <- -sum(S * (up + right))  # count each bond once (up and right)
  E / (L * L)
}

magnetization <- function(S) mean(S)

# ----- one temperature run --------------------------------------------------
ising_run <- function(T, L = 32, sweeps = 3000, burn = 1000) {
  S <- matrix(sample(c(-1, 1), L * L, replace = TRUE), L, L)
  N <- L * L
  
  # burn-in
  for (s in 1:burn) {
    for (step in 1:N) {
      i <- sample(L, 1); j <- sample(L, 1)
      dE <- local_dE(S, i, j, L)
      if (dE <= 0 || runif(1) < exp(-dE / T)) S[i, j] <- -S[i, j]
    }
  }
  
  # measurement sweeps
  m_vals <- numeric(sweeps)
  e_vals <- numeric(sweeps)
  for (s in 1:sweeps) {
    for (step in 1:N) {
      i <- sample(L, 1); j <- sample(L, 1)
      dE <- local_dE(S, i, j, L)
      if (dE <= 0 || runif(1) < exp(-dE / T)) S[i, j] <- -S[i, j]
    }
    m_vals[s] <- abs(magnetization(S))   # abs to avoid symmetry sign flips
    e_vals[s] <- energy_per_spin(S)
  }
  
  list(m = mean(m_vals), e = mean(e_vals))
}

# ----- sweep temperatures and plot ------------------------------------------
Ts <- seq(1.5, 3.5, length.out = 14)
res <- lapply(Ts, function(T) ising_run(T, L = 64, sweeps = 6000, burn = 2000))
mT  <- sapply(res, `[[`, "m")
eT  <- sapply(res, `[[`, "e")

op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
plot(Ts, mT, type = "b", pch = 16, xlab = "Temperature T",
     ylab = "Magnetization per spin", main = "2D Ising: m(T)")
abline(v = 2.269, lty = 2, col = "red")  # Tc (theory)

plot(Ts, eT, type = "b", pch = 16, xlab = "Temperature T",
     ylab = "Energy per spin", main = "2D Ising: e(T)")
abline(v = 2.269, lty = 2, col = "red")
par(op)
