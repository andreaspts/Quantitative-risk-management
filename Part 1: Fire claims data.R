#plot danish fire data

#load QRM library
library(QRM) 

plot(danish,
     xlab = "Year",
     ylab = "Insurance claims in million DKK")


#flag unusually large claims in the Danish fire losses data

#1. Threshold-based detection
#idea: Picks a high quantile or 
#an absolute value cutoff and marks points above it.

data(danish)
y <- as.numeric(danish)
dates <- time(danish)

# Example: threshold at 99th percentile
thr <- quantile(y, 0.99)
high_idx <- y > thr

plot(dates, y, type = "h",
     xlab = "Year", ylab = "Insurance claims (million DKK)",
     main = "Danish Fire Losses")
points(dates[high_idx], y[high_idx], col = "red", pch = 19)
abline(h = thr, col = "blue", lty = 2)

#2. Median Absolute Deviation (MAD) method
#idea: Identifies outliers.

med <- median(y)
mad_val <- mad(y)
high_idx <- y > med + 5 * mad_val  # "5" = sensitivity factor

plot(dates, y, type = "h",
     main = "Outliers via MAD",
     xlab = "Year", ylab = "Insurance claims (million DKK)")
points(dates[high_idx], y[high_idx], col = "red", pch = 19)

#3. EVT-based approach (fits a Generalized Pareto Distribution to tail)
#idea: Gives statistical parameters for the tail 
#and helps to define "unusually high".

library(evir)
gpd_fit <- gpd(y, threshold = thr)  # Fit GPD above 99% quantile
gpd_fit

