## =========================================
## GARCH workflow with robust interactive menus
## =========================================

graphics.off()

suppressPackageStartupMessages({
  library(tseries)
  library(zoo)
  # install.packages("fGarch")  # <- uncomment if needed
  library(fGarch)
})

download_data <- function() {
  pricel <- get.hist.quote("msft", start="2007-06-01", end="2025-12-31",
                           quote="AdjClose", quiet=TRUE)
  price2 <- get.hist.quote("ms",   start="2007-06-01", end="2025-12-31",
                           quote="AdjClose", quiet=TRUE)
  price3 <- get.hist.quote("GS",   start="2007-06-01", end="2025-12-31",
                           quote="AdjClose", quiet=TRUE)
  p <- cbind(pricel, price2, price3)
  colnames(p) <- c("MSFT","MS","GS")
  y <- na.omit(diff(log(p)))   # daily log returns
  list(p = p, y = y)
}

overview_plots <- function(p, y) {
  cols <- c("steelblue","firebrick","darkgreen")
  
  # Prices
  op <- par(mfrow = c(1,1), mar = c(4,4,3,1))
  plot.zoo(p, plot.type="single", col=cols, lwd=1.6,
           xlab="", ylab="Adjusted close",
           main="MSFT, MS, GS — Adjusted Close")
  grid(col="grey80", lty=1, lwd=0.8)
  legend("topleft", bty="n", lwd=1.6, col=cols, legend=colnames(p))
  par(op)
  
  # Returns (3 panels)
  op <- par(mfrow = c(3,1), mar = c(4,4,3,1))
  plot.zoo(y[,"MSFT"], col=cols[1], lwd=1.4, xlab="", ylab="log return",
           main="MSFT — Daily Log Returns"); grid(col="grey80", lty=1, lwd=0.8)
  plot.zoo(y[,"MS"],   col=cols[2], lwd=1.4, xlab="", ylab="log return",
           main="MS — Daily Log Returns");   grid(col="grey80", lty=1, lwd=0.8)
  plot.zoo(y[,"GS"],   col=cols[3], lwd=1.4, xlab="", ylab="log return",
           main="GS — Daily Log Returns");   grid(col="grey80", lty=1, lwd=0.8)
  par(op)
  
  invisible(NULL)
}

fit_one_series <- function(y, choice_idx) {
  lab  <- colnames(y)[choice_idx]
  yvec <- as.numeric(na.omit(y[, choice_idx]))
  
  fitN <- garchFit(~ garch(1,1), data=yvec, include.mean=FALSE,
                   cond.dist="norm", trace=FALSE)
  fitT <- garchFit(~ garch(1,1), data=yvec, include.mean=FALSE,
                   cond.dist="std",  trace=FALSE)
  
  aicN <- fitN@fit$ics["AIC"]; aicT <- fitT@fit$ics["AIC"]
  best <- if (aicT < aicN) fitT else fitN
  
  cat("\n====", lab, "GARCH(1,1) ====\n")
  cat("Normal AIC:", round(aicN, 4), " | Student-t AIC:", round(aicT, 4), "\n")
  cat("Chosen dist:", if (identical(best, fitT)) "Student-t" else "Normal", "\n\n")
  show(summary(best))
  
  list(best = best, label = lab)
}

diagnostics_menu <- function(best, label) {
  diag_labels <- c(
    "Time Series",
    "Conditional SD",
    "Series with 2 Conditional SD Superimposed",
    "ACF of Observations",
    "ACF of Squared Observations",
    "Cross Correlation",
    "Residuals",
    "Conditional SDs",
    "Standardized Residuals",
    "ACF of Standardized Residuals",
    "ACF of Squared Standardized Residuals",
    "Cross Correlation between r^2 and r",
    "QQ-Plot of Standardized Residuals",
    "Series with -VaR Superimposed",
    "Series with -ES Superimposed",
    "Series with -VaR & -ES Superimposed"
  )
  
  repeat {
    pick <- menu(diag_labels,
                 title = paste0("Diagnostics for ", label, " — pick 1..16; 0 to exit"))
    if (pick == 0) {
      message("Done with diagnostics.")
      break
    }
    plot(best, which = pick)  # show exactly the chosen diagnostic
  }
}


run_workflow <- function() {
  dat <- download_data()
  p <- dat$p; y <- dat$y
  
  overview_plots(p, y)
  
  choice <- menu(c("MSFT","MS","GS"),
                 title = "\nWhich series to fit with GARCH?")
  if (choice == 0) stop("Selection cancelled.")
  
  fit <- fit_one_series(y, choice)
  diagnostics_menu(fit$best, fit$label)
}

## --------- RUN (nothing comes after this call) ---------
run_workflow()
