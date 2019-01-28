" Model:

  posterior:
  p(μ|y) ~ p(y|μ).p(μ)

  likelihood: normal with known variance
  p(yi|μ) iid~ N(μ,1) for i=1,...,N
  
  prior:  t distribution
  p(μ) ~ t(0,1,1)
"


#install.packages("rjags")
library("rjags")

# Model specification for JAGS

mod_string = " model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu, 1.0/sig2)
  }
  mu ~ dt(0.0, 1.0/1.0, 1.0) # location, inverse scale, degrees of freedom
  sig2 = 1.0
} "


# Model setup

set.seed(50)
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9) # given sample
n = length(y)

data_jags = list(y=y, n=n) # data supplied to JAGS model, in the form of a list
params = c("mu")           # parameter to be observed

# initialise values - optional (and fixed)
inits = function() {
  inits = list("mu"=0.0)
} 

# create a JAGS model object
mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits) 

# Run MCMC sampler
update(mod, 500) # burn-in - run sampler for 500 iterations to reach stationary distribution (does not store values)

# run updated model for 1000 iterations
mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=1000)

summary(mod_sim)

# plot the Markov chain
library("coda")
plot(mod_sim)
