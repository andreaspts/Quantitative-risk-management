# Simulating paths of a random walk
set.seed(7)
T <- 200       # number of steps
m <- 10        # number of paths
y <- matrix(0, nrow = m, ncol = T)

for (ii in 1:m) {
  for (j in 2:T) {
    y[ii, j] <- y[ii, j - 1] + rnorm(1, mean = 0, sd = 1)
  }
}

# --- choose which path to analyze ---
sel <- 10                   # pick any index from 1 to m
stopifnot(sel >= 1, sel <= m)

# Re-plot all paths, highlight the chosen one
cols <- rep("grey70", m); cols[sel] <- "red"
plot(1:T, y[1, ], type = "l", lwd = ifelse(1 == sel, 2.2, 1.2),
     col = cols[1], ylim = range(y), xlab = "", ylab = "",
     main = sprintf("Random walk paths (highlight: path %d)", sel))
for (k in 2:m) {
  lines(1:T, y[k, ], col = cols[k], lwd = ifelse(k == sel, 2.2, 1.2))
}
legend("topleft", bty = "n", lwd = c(2.2, 1.2),
       col = c("red", "grey70"), legend = c("selected", "others"))

# --- Diagnostics for the selected path: ACF, PACF, QQ -----------------
x  <- y[sel, ]      # levels of the chosen random walk
dx <- diff(x)       # increments

op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Levels (nonstationary)
acf(x,  lag.max = 40, main = sprintf("ACF: path %d (levels)", sel))
pacf(x, lag.max = 40, main = sprintf("PACF: path %d (levels)", sel))
qqnorm(x, main = sprintf("QQ plot: path %d (levels)", sel)); qqline(x, col = "red")

# Increments (white noise if shocks ~ N(0,1))
acf(dx,  lag.max = 40, main = sprintf("ACF: path %d (increments)", sel))
pacf(dx, lag.max = 40, main = sprintf("PACF: path %d (increments)", sel))
qqnorm(dx, main = sprintf("QQ plot: path %d (increments)", sel)); qqline(dx, col = "red")

par(op)
