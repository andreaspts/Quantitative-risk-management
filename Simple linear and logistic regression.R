# linear and logistic regression in R

#linear regression
#simulate data from uniform distribution
#on the interval [-5,5]
x<- runif(n=50, min=-50, max=50)

#simulare error terms from the standard
#normal distribution
e<-rnorm(50)

y<- 0.2*x+e

#Estimation of the linear model:
#y=beta_hat*x+epsilon

library(stats)

#linear ols model
model<-lm(y~x)

summary(model)

lm(formula = y ~ x)

plot(x,y, main="Simple linear regression model")
abline(model$coefficients,col="blue")
abline(a=0,b=0.2, col="red")
legend(-4.5,2,legend=c("True parameters",
                       "Estimated parameters"), 
       col=c("red","blue"), lty=1:1, cex=1)

#logistic regression
# 1) Simulate predictor on [-5, 5]
n <- 300
x <- runif(n, min = -5, max = 5)

# 2) True logistic model: logit P(Y=1|x) = alpha + beta*x
alpha <- -0.5
beta  <-  0.8
eta   <- alpha + beta * x
p     <- 1 / (1 + exp(-eta))      # = plogis(eta)

# 3) Generate binary outcome
y <- rbinom(n, size = 1, prob = p)

# 4) Fit simple logistic regression (OLS analogue)
model <- glm(y ~ x, family = binomial(link = "logit"))
summary(model)

# 5) Plot points (jittered) + estimated and true probability curves
plot(x, jitter(y, amount = 0.06),
     pch = 16, col = adjustcolor("black", 0.4),
     xlab = "x", ylab = "y (0/1)",
     main = "Simple Logistic Regression")

xs   <- seq(min(x), max(x), length.out = 300)
phat <- predict(model, newdata = data.frame(x = xs), type = "response")
lines(xs, phat, col = "blue", lwd = 2)                     # estimated P(Y=1|x)

ptru <- plogis(alpha + beta * xs)
lines(xs, ptru, col = "red", lwd = 2, lty = 2)             # true P(Y=1|x)

legend("topleft",
       legend = c("Estimated P(Y=1|x)", "True P(Y=1|x)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

# 6) Quick interpretation helpers
coef(model)                 # log-odds coefficients (alpha_hat, beta_hat)
exp(coef(model))            # odds ratios
# e.g., exp(beta_hat): multiplicative change in odds per +1 in x

# 7) Optional: confusion table at cutoff 0.5
phat_obs <- fitted(model)
yhat     <- ifelse(phat_obs >= 0.5, 1, 0)
table(Observed = y, Predicted = yhat)
