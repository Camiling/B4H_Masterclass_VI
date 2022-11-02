## Practical exercise 1: Gaussian Mixture model 
rm(list=ls())
library(matrixcalc)

# Generate data
set.seed(123)
beta=c(-1,2)
phi=0.5
n=50
p=2
X = cbind(1, replicate(p-1, rnorm(n)))         # Generate x
y = X %*% beta  + rnorm(n, sd = sqrt(1/phi)) # Generate y
# Look at density
plot(density(y))

# Find mean-field estimates using CAVI
a0=b0=0.001
E_kappa = a0/b0 # Initial value
a = a0+p/2
more=TRUE
XX = crossprod(X,X)
Xy = crossprod(X,y)
yy = c(crossprod(y))
eps=0.0001
m_old = c(1000,1000)

while(more){
  S = solve(diag(E_kappa, p) + phi * XX)
  m = phi * S %*% Xy
  E_betabeta =  as.numeric(crossprod(m) + matrixcalc::matrix.trace(S))
  b = b0 + 0.5 * E_betabeta
  E_kappa = a / b
  # Assessing convergence by variational factors (can also use ELBO)
  if(sum((m-m_old)^2)<eps){
    more =FALSE
  }
  m_old = m
}

# Visualise bivariate Gaussian approximation for coefficients
library(mnormt)
set.seed(0)
beta1 = seq(-1.5, 2.5, 0.1)
beta2 = seq(-1.5, 2.5, 0.1)
f = function(beta1, beta2) dmnorm(cbind(beta1, beta2), c(m), S)
y.plot = outer(beta1, beta2, f)
# Create surface plot
persp(beta1, beta2, y.plot, theta=-20, phi=20, col = 'blue',
      expand=0.8, ticktype='detailed')

