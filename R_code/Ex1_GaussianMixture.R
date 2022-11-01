## Practical exercise 1: Gaussian Mixture model 
rm(list=ls())

# Generate data
set.seed(231171)
K = 3
mu = c(-1,1,3)
sig.mu = 1
tau.mu = 1/sig.mu^2
sig2 = c(1,1,1) 
n1 = 100 # For simplicity, 100 observations from each component
n = K*n1
x = rep(NA,n)
eps = 0.001
for(k in 1:K){
  x[(k-1)*n1+1:n1] = rnorm(n1,mu[k],sqrt(sig2[k]))
}
# Look at density
plot(density(x))

# Find mean-field estimates using CAVI
phi = matrix(1/K,nrow=n,ncol=K); m = mu; s2 = rep(1,K); more = TRUE
Elbo = 0
while(more){
  for(i in 1:n){
    phi[i,] = exp(m*x[i]-0.5*s2-0.5*m^2)
    phi[i,] = phi[i,]/sum(phi[i,])
  }
  for(k in 1:K){
    m[k] = sum(phi[,k]*x)/(tau.mu+sum(phi[,k]))
    s2[k] = 1/(tau.mu+sum(phi[,k]))
  }
  elbo = -0.5*tau.mu*sum(s2+m^2)-sum(rowSums(phi*log(phi)))+0.5*sum(log(s2))
  for(k in 1:K){
    elbo = elbo + sum(phi[,k]*(x*m[k]-0.5*s2[k]-0.5*m[k]^2))
  }
  more = abs(tail(Elbo,n=1)-elbo)>eps
  Elbo = c(Elbo,elbo)
}

# Final parameter estimates
phi[nrow(phi),]

# Look at convergence of ELBO
plot.ts(Elbo[-1])



