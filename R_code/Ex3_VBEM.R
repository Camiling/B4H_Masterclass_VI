## Practical exercise 3: Linear regression with empirical Bayes estimation for hyperparameters
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

# Initialize
a0=b0=0.001
E_kappa = a0/b0 # Initial value
# Initial value
phi_hat =0.2
a = a0+p/2
more=TRUE
more_phi = TRUE
XX = crossprod(X,X)
Xy = crossprod(X,y)
yy = c(crossprod(y))
eps=0.0001
m_old = c(1000,1000)
phi_old = 1000

while(more_phi){
  # E-step
  while(more){
    S = solve(diag(E_kappa, p) + phi_hat * XX)
    m = phi_hat * S %*% Xy
    E_betabeta =  as.numeric(crossprod(m) + matrixcalc::matrix.trace(S))
    b = b0 + 0.5 * E_betabeta
    E_kappa = a / b
    # Assessing convergence by variational factors (can also use ELBO)
    if(sum((m-m_old)^2)<eps){
      more =FALSE
    }
    m_old = m
  }
  # M-step
  phi_hat = as.numeric(n/(yy + matrixcalc::matrix.trace(XX%*%(m%*%t(m) + S)) - 2*t(m)%*%Xy))
  # Check for convergence
  more_phi = abs(phi_old-phi_hat) >= eps
  phi_old = phi_hat
}

# Resulting estimate for phi
phi_hat

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

