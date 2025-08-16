# Monte Carlo simulation applied to 2D Lennard–Jones fluid

# 2D Lennard–Jones fluid (reduced units: epsilon = sigma = 1)
# Metropolis MC in NVT, periodic boundary conditions

set.seed(1)

run_lj_mc <- function(N = 64, rho = 0.6, Tstar = 1.0,
                      sweeps = 3000, burn = 1000,
                      max_disp = 0.15, rc = 2.5, nbins = 120) {
  
  stopifnot(N > 3, rho > 0, Tstar > 0, sweeps > burn, nbins > 10)
  L   <- sqrt(N / rho)          # box length (2D)
  beta <- 1 / Tstar
  rc2  <- rc^2
  
  # minimum-image displacement
  mic <- function(dx) dx - L * round(dx / L)
  
  # LJ pair energy with potential shift at rc (so U(rc) = 0)
  lj_pair_shifted <- function(r2) {
    inv_r2 <- 1 / r2
    inv_r6 <- inv_r2^3
    u <- 4 * (inv_r6^2 - inv_r6)     # 4*(1/r^12 - 1/r^6)
    u
  }
  # shift value at cutoff
  inv_rc2 <- 1 / rc2
  inv_rc6 <- inv_rc2^3
  u_shift <- 4 * (inv_rc6^2 - inv_rc6)
  
  # lattice-ish initial positions
  nx <- ceiling(sqrt(N))
  ny <- ceiling(N / nx)
  xs <- (0:(nx - 1)) * (L / nx) + 0.5 * (L / nx)
  ys <- (0:(ny - 1)) * (L / ny) + 0.5 * (L / ny)
  grid <- as.matrix(expand.grid(xs, ys))[1:N, ]
  pos  <- matrix(grid, ncol = 2)
  
  # total potential energy (i<j)
  total_energy <- function(pos) {
    E <- 0
    for (i in 1:(N - 1)) {
      dx <- mic(pos[i, 1] - pos[(i + 1):N, 1])
      dy <- mic(pos[i, 2] - pos[(i + 1):N, 2])
      r2 <- dx * dx + dy * dy
      ok <- r2 < rc2
      if (any(ok)) {
        u <- lj_pair_shifted(r2[ok]) - u_shift
        E <- E + sum(u)
      }
    }
    E
  }
  
  # energy contribution of particle i with all j != i
  energy_i <- function(i, pos) {
    jset <- setdiff(1:N, i)
    dx <- mic(pos[i, 1] - pos[jset, 1])
    dy <- mic(pos[i, 2] - pos[jset, 2])
    r2 <- dx * dx + dy * dy
    ok <- r2 < rc2
    if (!any(ok)) return(0)
    sum(lj_pair_shifted(r2[ok]) - u_shift)
  }
  
  # initialize energy
  E <- total_energy(pos)
  
  # g(r) accumulators
  rmax <- L / 2
  dr   <- rmax / nbins
  counts <- numeric(nbins)
  samples <- 0L
  E_traj <- numeric(sweeps - burn)
  
  # MC sweeps
  for (s in 1:sweeps) {
    order <- sample.int(N, N, replace = FALSE)
    for (idx in order) {
      i <- idx
      # old energy contribution
      Ei_old <- energy_i(i, pos)
      
      # propose move
      trial <- pos[i, ] + runif(2, -max_disp, max_disp)
      # periodic wrap
      trial <- trial %% L
      
      # compute new contribution
      pos_i_old <- pos[i, ]
      pos[i, ] <- trial
      Ei_new <- energy_i(i, pos)
      
      dE <- Ei_new - Ei_old
      if (dE <= 0 || runif(1) < exp(-beta * dE)) {
        # accept: E already updated via dE (pos kept as trial)
        E <- E + dE
      } else {
        # reject
        pos[i, ] <- pos_i_old
      }
    }
    
    # measurements
    if (s > burn) {
      # energy per particle
      E_traj[s - burn] <- E / N
      
      # radial distribution counts (unordered pairs i<j)
      for (i in 1:(N - 1)) {
        dx <- mic(pos[i, 1] - pos[(i + 1):N, 1])
        dy <- mic(pos[i, 2] - pos[(i + 1):N, 2])
        r  <- sqrt(dx * dx + dy * dy)
        ok <- r < rmax
        if (any(ok)) {
          bin <- pmax(1L, pmin(nbins, floor(r[ok] / dr) + 1L))
          for (b in bin) counts[b] <- counts[b] + 1L
        }
      }
      samples <- samples + 1L
    }
  }
  
  # build g(r) normalization (2D):
  # counts are over unordered pairs per configuration.
  # For ideal gas: expected pairs in shell = (N * rho * 2*pi*r*dr) / 2
  r_mid <- (1:nbins - 0.5) * dr
  rho_  <- N / (L * L)
  norm  <- (samples * (N) * rho_ * 2 * pi * r_mid * dr) / 2
  g     <- counts / norm
  
  list(r = r_mid, g = g, E_traj = E_traj, L = L, rho = rho_, Tstar = Tstar)
}

# ---- Run and plot ---------------------------------------------------------
res <- run_lj_mc(N = 64, rho = 0.6, Tstar = 1.0,
                 sweeps = 2500, burn = 800,
                 max_disp = 0.15, rc = 2.5, nbins = 120)

op <- par(mfrow = c(1,2), mar = c(4,4,3,1))
plot(res$r, res$g, type = "l", xlab = "r", ylab = "g(r)",
     main = "2D Lennard–Jones: radial distribution")
abline(h = 1, lty = 2, col = "grey70")

plot(res$E_traj, type = "l",
     xlab = "measurement sweep", ylab = "E/N",
     main = "Potential energy per particle")
par(op)
