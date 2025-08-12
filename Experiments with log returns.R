#plot DBAG log returns
library("tseries")

#downloading stock price data from Yahoo Finance
price=get.hist.quote(instrument="DBK.DE",start="2007-01-01",
                     quote="AdjClose") 

#calculate log returns
returns=diff(log(price))

#plot the calculated log returns
plot(returns, xlab="year", ylab="log returns", 
     main="DB AG log returns")

#simple statistics to get a feeling for the data
mean(returns)
sd(returns)
min(returns)
max(returns)
round(summary(as.vector(returns)),4)
skewness(returns)
kurtosis(returns)
# test for normal distribution H_0 
# (if p smaller than threshold of 5% then normality rejected).
jarque.bera.test(returns)

#understand autocorrelation of logreturns and logreturns^2
acf(as.vector(returns), lag.max=5, plot=FALSE)
acf(as.vector(returns^2), lag.max=5, plot=FALSE)

#boxplot and normal QQ plot
boxplot(returns,
        main = "DBK.DE log returns",
        ylab = "r_t",
        notch = TRUE,          # shows median CI notch
        horizontal = TRUE)     # nicer for one variable

qqnorm(returns, main = "Normal Q-Q plot")
qqline(returns, col = "red", lwd = 2)

d<-density(returns)
plot(d, main = "Density of log-returns", xlab = "r_t")
# Overlay Normal with same mean/sd (for tail check)
m <- mean(returns); s <- sd(returns)
lines(d$x, dnorm(d$x, mean = m, sd = s), lty = 4)




#correlation with dax30 logreturns
# DBK.DE vs DAX30: scatter with linear fit
library(tseries)
library(zoo)

# 1) Download adjusted close prices from Yahoo
dbk <- get.hist.quote(instrument = "DBK.DE",
                      start = "2007-01-01",
                      quote = "AdjClose",
                      retclass = "zoo")

dax <- get.hist.quote(instrument = "^GDAXI",
                      start = "2007-01-01",
                      quote = "AdjClose",
                      retclass = "zoo")

# 2) Compute daily log returns
r_dbk <- diff(log(dbk))
r_dax <- diff(log(dax))

# 3) Align by date (intersection) and drop NAs
both <- na.omit(merge(r_dax, r_dbk, all = FALSE))
colnames(both) <- c("DAX", "DBK")

# 4) Scatter plot (DBK on y, DAX on x) + linear fit
plot(both$DAX, both$DBK,
     xlab = "DAX30 daily log return",
     ylab = "DBK.DE daily log return",
     main = "DBK.DE vs DAX30: Daily Log Returns")

fit <- lm(DBK ~ DAX, data = as.data.frame(both))
abline(fit, col = "red", lwd = 2)

# Optional: show regression equation and R^2
co <- coef(fit); r2 <- summary(fit)$r.squared
legend("topleft", bty = "n",
       legend = c(
         sprintf("DBK = %.4f + %.2f × DAX", co[1], co[2]),
         sprintf("R² = %.3f", r2)
       ))

# Pearson correlation (DBK vs DAX)
df <- na.omit(data.frame(
  DAX = as.numeric(both$DAX),
  DBK = as.numeric(both$DBK)
))

# point estimate
r <- cor(df$DAX, df$DBK, method = "pearson")
r

#a 92-day rolling Pearson r
library(zoo)
roll_r <- rollapply(both, width = 92,
                    FUN = function(z) cor(z[,1], z[,2]),
                    by.column = FALSE, align = "right")
plot(roll_r, main = "92-day rolling Pearson r (DBK vs DAX)",
     ylab = "r", xlab = "Date", type = "l")

# hypothesis test with 95% CI
ct <- cor.test(df$DAX, df$DBK, method = "pearson")
ct$estimate         # Pearson's r
ct$p.value          # p-value
ct$conf.int         # 95% CI

#check for volatility clustering of the DBAG logreturn data

#1. compute rolling volatility
#idea: Compute a rolling (e.g., 30‑day) annualized stdev and plot it. 
#Spikes and long high-vol stretches = clusters.

library(zoo)

r <- as.numeric(returns)
d <- as.Date(time(returns))

vol30 <- sqrt(252) * rollapply(r, 30, sd, align = "right", fill = NA)

op <- par(mfrow = c(2,1), mar = c(3,4,2,1))
plot(d, r, type = "h", main = "DBK.DE log returns", xlab = "", ylab = "r_t")
abline(h = 0, col = "grey")

plot(d, vol30, type = "l", main = "Rolling volatility (30-day, annualized)",
     xlab = "year", ylab = "σ̂_t")
par(op)

#2. compute ACF of squared returns
#idea: Clustering → significant autocorrelation in r_t^2. 
#The ARCH‑LM test checks this formally.

acf(r^2, main = "ACF of squared returns")

# ARCH-LM test
# install.packages("FinTS")  # if needed
library(FinTS)
ArchTest(r, lags = 12)

#3. compute GARCH(1,1) conditional volatility
#idea: Fit a simple GARCH.
#The fitted conditional sigma(t) shows clusters very clearly.

# install.packages("rugarch")  # if needed
library(rugarch)

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "norm"
)

fit <- ugarchfit(spec, data = r)
sig <- sigma(fit)

plot(d, sig, type = "l", xlab = "year", ylab = "σ_t",
     main = "GARCH(1,1) conditional volatility")

#4. identify high-volume regimes
#idea: Highlight “high‑vol” regimes.

thr <- quantile(vol30, 0.8, na.rm = TRUE)
hi <- vol30 >= thr
plot(d, vol30, type = "l", xlab = "year", ylab = "σ̂_t",
     main = "Rolling vol with high-vol regimes")
lines(d[hi], vol30[hi], type = "h")  # spikes where vol is high
abline(h = thr, lty = 2)