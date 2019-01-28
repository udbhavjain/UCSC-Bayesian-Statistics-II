" posterior distribution:
  
  p(μ,σ^2∣yi) ∝ p(yi|μ,σ^2).p(μ).p(σ^2)

  Likelihood:
  p(yi|μ,σ^2) iid~ N(μ,σ^2) ; i = 1,...,n

  Here, we consider both mean and variance to be unknown. 
  In Gibbs sampling, the values for both will be simulated and used as inputs for generating 
  full conditionals of each other. 

  Priors:
  μ ~ N(μ0,σ0^2)
  σ^2 ~ IG(ν0,β0)

  When variance is known, normal is conjugate prior for mean.
  When mean is known, inverse gamma is conjugate prior for variance.


  p(μ|σ^2,yi) ∝ p(μ,σ^2∣yi)

  Expanding p(yi|μ,σ^2).p(μ).p(σ^2) and dropping terms without μ gives:

  p(μ|σ^2,y1,...,yn) ∝  N (  μ |    n.Ybar/σ^2 + μ0/σ0^2   ,        1        )
                                     -------------------       -------------
                                       n/σ^2 + 1/σ^2           n/σ^2 + 1/σ^2



  Similarly,

  p(σ^2|μ,yi) ∝ p(μ,σ^2∣yi)
  
  Expanding p(yi|μ,σ^2).p(μ).p(σ^2) and dropping terms without σ^2 gives:

  p(σ^2|μ,y1,...,yn) ∝ IG ( σ^2 | ν0 + n/2 , β0 + (Σ(yi-μ)^2)/2 ) ; i = 1,...,n

"

# full conditional for mean
update_mu = function(n, ybar, sig2, mu_0, sig2_0) 
{
  "n = size of given sample
   ybar = mean of given sample
   mu_0 = hyperparameter - mean for normal prior distribution
   sig2_0 = hyperparameter - variance for normal prior distribution"
  
  # variance
  sig2_1 = 1.0 / (n / sig2 + 1.0 / sig2_0)
  # mean
  mu_1 = sig2_1 * (n * ybar / sig2 + mu_0 / sig2_0)
  # simulate a value from full conditional distribution defined by the parameters above
  rnorm(n=1, mean=mu_1, sd=sqrt(sig2_1))
}

# full conditional for variance
update_sig2 = function(n, y, mu, nu_0, beta_0) 
{
  "n = size of given sample
   y = given sample values
   mu = value of mean for current iteration
   nu_0 = hyperparameter - shape for prior inverse-gamma distribution
   beta_0 = hyperparameter - rate for prior inverse-gamma distribution"
  
  # shape
  nu_1 = nu_0 + n / 2.0
  # rate
  sumsq = sum( (y - mu)^2 ) 
  beta_1 = beta_0 + sumsq / 2.0
  
  # simulate a value from gamma distribution defined by the parameters above
  out_gamma = rgamma(n=1, shape=nu_1, rate=beta_1)
  
  # convert to inverse gamma
  1.0 / out_gamma 
}

# function for Gibbs sampling
gibbs = function(y, n_iter, init, prior) 
{
  "y = given sample values
   n_iter = number of iterations for Markov chain
   init = list with initial state for mu
   prior = list with hyperparameters for prior normal and inverse-gamma distributions
  "
  ybar = mean(y)
  n = length(y)
  
  ## initialize
  mu_out = numeric(n_iter)
  sig2_out = numeric(n_iter)
  
  mu_now = init$mu # initialise mean with supplied initial state
  
  ## Gibbs sampler
  for (i in 1:n_iter) {
    # simulate new values for mean and variance using current values of each other
    sig2_now = update_sig2(n=n, y=y, mu=mu_now, nu_0=prior$nu_0, beta_0=prior$beta_0)
    mu_now = update_mu(n=n, ybar=ybar, sig2=sig2_now, mu_0=prior$mu_0, sig2_0=prior$sig2_0)
    
    # store generated values
    sig2_out[i] = sig2_now
    mu_out[i] = mu_now
  }
  
  cbind(mu=mu_out, sig2=sig2_out)
}

# given sample
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)

# hyperparameters for priors 
prior = list()

# hyperparameters for normal prior
prior$mu_0 = 0.0
prior$sig2_0 = 1.0

# hyperparameters for inverse-gamma prior
prior$n_0 = 2.0 # prior effective sample size for sig2
prior$s2_0 = 1.0 # prior point estimate for sig2
prior$nu_0 = prior$n_0 / 2.0 # prior parameter for inverse-gamma
prior$beta_0 = prior$n_0 * prior$s2_0 / 2.0 # prior parameter for inverse-gamma

# compare prior for mu and distribution of data
hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dnorm(x=x, mean=prior$mu_0, sd=sqrt(prior$sig2_0)), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean


# run the sampler
set.seed(53)

init = list()
init$mu = 0.0 # initial state for mu

# run Gibbs sampler for 1000 iterations
post = gibbs(y=y, n_iter=1e3, init=init, prior=prior)
head(post)

library("coda")
plot(as.mcmc(post))

summary(as.mcmc(post))
