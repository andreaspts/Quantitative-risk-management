#simulate one path/trajectory of a standard Wiener process in discrete time

set.seed (10)

T=1; N=500; dt=T/N 

dW=rep(0, N) #change in the Wiender processs
W=rep(0,N) # Wiener process

dW[1]=sqrt(dt) *rnorm(1) #W_t=\sqrt{t}Z_t
W[1] =dW[1]

for (j in 2:N) {
  dW[j]=sqrt(dt) *rnorm(1) #increment evaluation
  W[j]=W[j-1]+dW[j]        #add increment to previous value of W
}

plot(seq(0,T,by=dt), c(0,W),type="l",xlab="t",
      ylab="W(t)",col="blue")