" Mixture distribution: weighted combination of probability distribtuions. 
  
  Example: Exponential distribution with mean 1 and normal distribution with mean 3 
  and variance 1, with weights 0.4 and 0.6 respectively. 

  p(y) = (0.4).exp(-y).I(y>=0) + (0.6).(1/√2π).exp(-((y-3)^2)/2)

  
  General form of discrete mixture of distributions is:

  p(y) = nΣi=1 wi.fi(y)

  Where wi is weight and fi(y) is a PDF.
  The weights are probabilities and their sum is 1.

"

# PDF of distribution
curve( 0.4*dexp(x, 1.0) + 0.6*dnorm(x, 3.0, 1.0), from=-2.0, to=7.0, ylab="density", 
       xlab="y", main="40/60 mixture of exponential and normal distributions", lwd=2)

# PDF for each population
curve( 0.4*dexp(x, 1.0) + 0.6*dnorm(x, 3.0, 1.0), from=-2.0, to=7.0, ylab="density",
       xlab="y", main="40/60 mixture of exponential and normal distributions", lwd=2)
curve( 0.4*dexp(x, 1.0), from=-2.0, to=7.0, col="red", lty=2, add=TRUE)
curve( 0.6*dnorm(x, 3.0, 1.0), from=-2.0, to=7.0, col="blue", lty=2, add=TRUE)

# simulate mixture distribution
set.seed(117)

n = 1000
z = numeric(n)
y = numeric(n)

for (i in 1:n) {
  z[i] = sample.int(2, 1, prob=c(0.4, 0.6)) # returns a 1 with probability 0.4, or a 2 with probability 0.6
  if (z[i] == 1) {
    y[i] = rexp(1, rate=1.0)
  } else if (z[i] == 2) {
    y[i] = rnorm(1, mean=3.0, sd=1.0)
  }
}

hist(y, breaks=30)

" This can be written as:

  p(y) = 2Σj=1 p(y,z=j) = 2Σj=1 p(z=j).p(y|z=j) = 2Σj=1 wj.fj(y) 

  
  When a mixture model is fitted, we only know about the y variables, and not the
  population they belong to i.e z is unknown. Since the z variables are not observed,
  they are called latent variables. They can be treated as parameters in a hierarchical
  model and bayesian inference can be performed for them.

  Hierarchical model:
  
  yi|zi,θ ind~ fzi(y|θ) ; i = 1,...,n
  Pr(zi=j|w) = wj ; j = 1,...,J
  
  w ~ p(w)
  θ ~ p(θ)

  A Dirichlet prior may be used for weight and population specific parameters in θ.

"

# load dataset
dat = read.csv("mixture.csv", header=FALSE)
y = dat$V1
(n = length(y))

hist(y, breaks=20)

plot(density(y))

" The data appears to be from two distributions. A mixture model of two normal
  distributions can be used to learn more. "

# JAGS Model
library("rjags")

mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[z[i]], prec)
z[i] ~ dcat(omega)
}

mu[1] ~ dnorm(-1.0, 1.0/100.0)
mu[2] ~ dnorm(1.0, 1.0/100.0) T(mu[1],) # ensures mu[1] < mu[2]

prec ~ dgamma(1.0/2.0, 1.0*1.0/2.0)
sig = sqrt(1.0/prec)

omega ~ ddirich(c(1.0, 1.0))
} "

set.seed(11)

data_jags = list(y=y)

params = c("mu", "sig", "omega", "z[1]", "z[31]", "z[49]", "z[6]") # Select some z's to monitor

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim, ask=TRUE)

autocorr.diag(mod_sim)
effectiveSize(mod_sim)

# model summary
summary(mod_sim)

# posterior densities for the population parameters and the mixing weights
par(mfrow=c(3,2))
densplot(mod_csim[,c("mu[1]", "mu[2]", "omega[1]", "omega[2]", "sig")])

# for the z's
par(mfrow=c(2,2))
densplot(mod_csim[,c("z[1]", "z[31]", "z[49]", "z[6]")])


# posterior probabilities for z[1], the membership of y[1]
table(mod_csim[,"z[1]"]) / nrow(mod_csim) 

# posterior probabilities for z[31], the membership of y[31]
table(mod_csim[,"z[31]"]) / nrow(mod_csim) 

# posterior probabilities for z[49], the membership of y[49]
table(mod_csim[,"z[49]"]) / nrow(mod_csim) 

# posterior probabilities for z[6], the membership of y[6]
table(mod_csim[,"z[6]"]) / nrow(mod_csim) 


y[c(1, 31, 49, 6)]

" y1 is clearly in Population 1's territory, y31 is ambiguous, y49 is ambiguous but 
  is closer to Population 2's territory, and y6 is clearly in Population 2's territory.
  The posterior distributions for the z variables closely reflect our assessment. "
