" Consider the posterior distribution:
  
  p(μ∣y1,...,yn) ∝  exp[n( Ybar*μ − (μ^2/2) )]
                    -------------------------
                            1+μ^2

  The expression on right is the product of normal likelihood with sd=1 and t-distribution prior 
  with mean = 0, scale = 1 and df = 1, with the constants removed.

  As g(μ) can take very small values and cause problems when calculating alpha, we will use log(g(μ))

  log(g(μ))=n(Ybar*μ− (μ^2/2))−log(1+μ^2)

"

# function to calculate log of g()
lg = function(mu, n, ybar) {
  mu2 = mu^2
  n * (ybar * mu - mu2 / 2.0) - log(1 + mu2)
}

# Random-Walk Metropolis-Hastings algorithm

mh = function(n, ybar, n_iter, mu_init, cand_sd) {
  
  "
    n = given sample size
    ybar = given sample mean
    n_iter = number of simulations
    mu_init = initial Markov chain state
    cand_sd = standard deviation for proposal distribution
    
  "
    ## step 1, initialize
  mu_out = numeric(n_iter) # for storing (accepted) generated draws 
  accpt = 0                # for storing number of accepted draws
  mu_now = mu_init         # current head of Markov chain, initialised with given sample mean
  lg_now = lg(mu=mu_now, n=n, ybar=ybar) # current log(g(mu)), initialised with supplied Markov chain state, 
                                         # given sample size, and given sample mean
  
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## step 2a
    
    # draw a candidate from proposal distribution 
    # (normal here, with mean equal to current markov chain state and supplied proposal standard deviation) 
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd) 
    
    ## step 2b
    lg_cand = lg(mu=mu_cand, n=n, ybar=ybar) # evaluate log of g with the candidate
    lalpha = lg_cand - lg_now # log of acceptance ratio i.e log(g(MUi)/g(MUi-1))
    alpha = exp(lalpha) # convert log of alpha to real alpha
    
    ## step 2c
    u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha)
    if (u < alpha) { # then accept the candidate
      mu_now = mu_cand
      accpt = accpt + 1 # to keep track of acceptance
      lg_now = lg_cand
    }
    
    ## collect results
    mu_out[i] = mu_now # save this iteration's value of mu
  }
  
  ## return a list of output
  list(mu=mu_out, accpt=accpt/n_iter)
}


# initial data for the function

# given sample
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)

hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean


# run the sampler for 1000 iterations, proposal standard deviation = 3.0 and initial Markov chain state = 0.0

set.seed(43)

post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=3.0)
str(post)

# plot Markov chain generation history (simulated values vs number of iterations)
library("coda")

traceplot(as.mcmc(post$mu))

" Acceptance rate is 12.2%, which is below 23%. Step sizes need to be reduced. "

# draw samples with proposal sd = 0.05
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.05)
post$accpt

traceplot(as.mcmc(post$mu))

" Acceptance rate is 93.7%, which is higher than 50%. Step sizes need to be increased. "

post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.9)
post$accpt

traceplot(as.mcmc(post$mu))

" Acceptance rate is 36.7%, which is between 23% and 50%. "


# initialise Markov chain with a value different from sample mean
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=30.0, cand_sd=0.9)
post$accpt

traceplot(as.mcmc(post$mu))

" The chain eventually hits stationary distribution. "

# compare posterior and prior densities

post$mu_keep = post$mu[-c(1:100)] # discard the first 100 samples

# plot density estimate of the posterior
plot(density(post$mu_keep, adjust=2.0), main="", xlim=c(-1.0, 3.0), xlab=expression(mu)) 
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu
points(ybar, 0, pch=19) # sample mean
curve(0.017*exp(lg(mu=x, n=n, ybar=ybar)), from=-1.0, to=3.0, add=TRUE, col="blue") # approximation to the true posterior in blue