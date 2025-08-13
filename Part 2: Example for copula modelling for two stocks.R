#Motivation: Using the theory of copulas, we will apply it 
#to the return series of two large US stocks 
#and model their joint behavior in R. 
#We therefore download stock data from Yahoo Finance, 
#calculate pseudo observations from the log returns, 
#determine a suitable copula model,
#fit the copula model to the pseudo observations,
#model the marginal distributions 
#and simulate new observations from the fitted model.

# Download stock data from Yahoo finance

library(tseries) #to download stock timeseries

price1 = get.hist.quote(instrument="msft" ,start="2010-01-01",
                        end="2020-04-01", quote="AdjClose" )
price2 = get.hist.quote(instrument="ms",start="2010-01-01",
                        end="2020-04-01", quote="AdjClose" )
# Calculate log returns
return1<-diff (log(price1)) 
return2<-diff (log (price2))

plot(return1, return2,pch=16, cex=0.7,xlab="Microsoft",
      ylab="Morgan Stanley")
abline(lm (return2~return1),col='blue')

cor(return1, return2, method='pearson')
## [1] 0.518635
cor(return1, return2, method='kendall')
## [1] 0.3061241
cor(return1, return2, method='spearman' )
## [1] 0.4355747

library(copula) #for measuring tail dependence 
pseudoObservationMatrix<-as.matrix(cbind(pobs(return1),pobs(return2)))
fitLambda (u=pseudoObservationMatrix, lower.tail=TRUE) [2, 1]
## [1] 0.4191571
fitLambda (u=pseudoObservationMatrix, lower.tail=FALSE) [2,1]
## [1] 0.25483

#Determining a suitable copula model
#As some copulas are better suited for modelling 
#an asymmetrical dependence and others put more 
#emphasis on tail dependence, it is very important 
#to choose a suitable copula model that can 
#adequatly capture the dependence structure.
#The VineCopula R-package offers the function 
#BiCopSelect to perform copula selection based on 
#BIC and AlC for bivariate data.
#Estimates are obtained by maximum likelihood estimation.
#The parameter familyset of the BiCopSelect function 
#allows for specifying a subset of copula families 
#to be considered.

#install.packages("VineCopula") 
library (VineCopula)
#creating pseudo observations
u<-pobs (return1)
v<-pobs (return2)

#selecting a suitable copula model
# The familyset parameter corresponds to:
# 0 = independence copula
# 1 = Gaussian copula
# 2 = Student t copula (t-copula)
# 3 = Clayton copula
# 4 = Gumbel copula
# 5 = Frank copula
# 6 = Joe copula
copulaModel<-BiCopSelect (u,v, familyset=0: 6)

summary(copulaModel)

#fitting a t-copula to the data as suggested
library (copula)

tCop<-tCopula (dim=2)

pseudoobservationMatrix<-as.matrix(cbind(pobs(return1),
                                         pobs(return2)))

fittedModel<-fitCopula (tCop,pseudoObservationMatrix,
                        method="ml")

#density plot of the fitted copula
rho=coef(fittedModel) [1]
df=coef(fittedModel) [2]
persp(tCopula (dim=2,param=rho, df=df),dCopula)

#Sampling data points from the estimated t-copula
#• We sample pseudo observations from the estimated t-copula.
#However, the pseudo observations have 
#the same correlation as the original returns.
#• The t-copula is particularly well suited 
#for modeling extreme events especially 
#when there is high correlation 
#in the tails of the distribution.
#•As the t-copula is symmetric it is not suited 
#to capture asymmetries in the data. 
#However, we are neglecting this issue in 
#our little application.
#•Before pseudo observations can be compared 
#to the actual returns they have to be 
#re-transformed based on a marginal distribution. 
#This is done in the following.

nSim=length(return1) #number of pseudo observations
pseudoObservations<-rCopula(nSim,
                             tCopula(dim=2, param=rho,df=df))
plot(pseudoObservations[,1],pseudoObservations[,2],
      pch=16, cex=0.7,xlab="Microsoft",ylab="Morgan Stanley")

#• We model the marginals based on a Student-t distribution.
#•This is sufficient in our setting 
#as our emphasis is on modeling the joint distribution 
#of Microsoft and Morgan Stanley stock returns 
#and not on modeling the univariate marginal distributions.

#Fitting univariate distributions to the marginals
#install.packages("MASS") 
#install.packages("metRology") 
library (MASS) #For fitting Student-t-distribution
library (metRology) #scaled t quantile function

fit1<-fitdistr(return1, "t" )
fit2<-fitdistr(return2, "t")

simReturns1<-qt.scaled(mean=coef (fit1) [1],
                       sd=coef(fit1) [2],df=coef(fit1) [3],
                       p=pseudoObservations [,1])
simReturns2<-qt.scaled(mean=coef(fit2) [1],
                       sd=coef(fit2) [2] ,df=coef(fit2) [3] ,
                       p=pseudoObservations [,2])

#plotting the actual and simulated returns
plot(return1, return2, xlab="Microsoft", ylab="Morgan Stanley")
points (simReturns1, simReturns2,col="red" )
legend("topleft",legend=c ("observed", "simulated") ,
        col=c("black", "red"),pch=1)

#one would argue that this model is suitable to model 
#the joint behavior of the two stocks
