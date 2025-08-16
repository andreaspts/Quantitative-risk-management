#Monte Carlo simulation application to integration

set.seed(1)

mc_integral <- function(f, a, b, n = 1e5, antithetic = FALSE) {
  if (!antithetic) {
    x  <- runif(n, a, b)
    fx <- f(x)
    Ihat <- (b - a) * mean(fx)
    se   <- (b - a) * sd(fx) / sqrt(n)
    ci   <- Ihat + c(-1, 1) * 1.96 * se
    return(list(estimate = Ihat, se = se, ci95 = ci))
  } else {
    m <- ceiling(n / 2)
    u <- runif(m)
    x1 <- a + (b - a) * u
    x2 <- a + (b - a) * (1 - u)        # antithetic mate
    fx <- (f(x1) + f(x2)) / 2
    Ihat <- (b - a) * mean(fx)
    se   <- (b - a) * sd(fx) / sqrt(m) # m pairs
    ci   <- Ihat + c(-1, 1) * 1.96 * se
    return(list(estimate = Ihat, se = se, ci95 = ci))
  }
}

f <- function(x) exp(-x^2)

base   <- mc_integral(f, 0, 1, n = 1e5)
anti   <- mc_integral(f, 0, 1, n = 1e5, antithetic = TRUE)

base
anti     # usually narrower CI (lower variance)

#Importance sampling on [0,\infty)
set.seed(1)
n <- 1e5
x  <- rexp(n, rate = 1)            # importance proposal matches e^{-x}
w  <- 1 / (1 + x^2)                # integrand under the proposal
Jhat <- mean(w)
se   <- sd(w) / sqrt(n)
ci   <- Jhat + c(-1, 1) * 1.96 * se
list(estimate = Jhat, se = se, ci95 = ci)


# MC application to differentiation of an integral

set.seed(1)

mc_deriv_integral <- function(theta, n = 1e5) {
  u <- runif(n)
  f  <- exp(-theta * u^2)             # integrand
  df <- -u^2 * exp(-theta * u^2)      # pathwise derivative
  
  Ihat  <- mean(f)
  dIhat <- mean(df)
  
  se_I  <- sd(f)  / sqrt(n)
  se_dI <- sd(df) / sqrt(n)
  
  ci_I   <- Ihat  + c(-1, 1) * 1.96 * se_I
  ci_dI  <- dIhat + c(-1, 1) * 1.96 * se_dI
  
  list(I = Ihat, I_ci = ci_I, dI = dIhat, dI_ci = ci_dI)
}

res <- mc_deriv_integral(theta = 1, n = 1e5)
res

