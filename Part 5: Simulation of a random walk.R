#Simulating paths of a random walk

set.seed(7)
T=200     #number of steps
m= 5      #number of paths to simulate
y=matrix(0,nrow=m, ncol=T) #initializing matrix to store

#simulated paths
for (i in 1:m) {
  for (j in 2:T) {
    #calculate paths
    y[i,j]=y[i, j-1]+rnorm(n=1, mean=0, sd=1)
  }
}

#plotting
plot(1:T, y[1,],xlim=c(1,T),type="l", lwd=1.5,
              ylim=c(min(y),max(y)),xlab="",ylab="",main="Simulated
paths of a random walk")
  for(i in 2:m) {
    lines (1:T, y[i,],col=i,lwd=1.5) 
    }