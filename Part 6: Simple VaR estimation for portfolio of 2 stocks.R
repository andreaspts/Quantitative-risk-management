#VaR estimation - practical example

library ("tseries")

#retrieve stock data
p1 = get.hist.quote(instrument="msft",
                    start="2000-01-01", end="2009-12-31", quote="AdjClose")
p2 = get.hist.quote(instrument="ibm",
                      start="2000-01-01", end="2009-12-31", quote="AdjClose")

#convert to log returns
y1=coredata(diff(log(p1)))
y2=coredata(diff(log(p2)))
y1=tail(y1, T-14)
y2=tail(y2, T-14)
T=length(y1)

#Portfolio value 
value=1000
y=cbind(y1, y2)
#Confidence level 1-p
p=0.01

#evaluating a multivariate ES using historical sumulation
ys=sort(y1) # sort returns 
op=T*p # p % smallest
VaR1=-ys[op]*value # Calculating the VaR level

#portfolio weights 
w=matrix(c(0.3,0.7))
#portfolio returns
yp=y%*%w
yps=sort(yp)
VaR2=-yps[op] *value

#calculating a multivariate ES using historical simulatio
ES1 = -mean (ys [1: op]) *value

#univariate VaR estimation based on a GARCH-normal model
library (fGarch)
g=garchFit(~garch (1,1),y1, cond.dist="norm", include.mean=FALSE, trace=FALSE)
omega=g@fit$matcoef[1, 1]
alpha=g@fit$matcoef[2,1]
beta=g@fit$matcoef[3,1]
# Calculate sigma2 for t+1 (forecast)
sigma2=omega + alpha*y[T]^2 + beta*g@h.t[T]
#VaR forecast
VaR9=-sqrt(sigma2)*qnorm(p)*value



