# use Metropolis-Hastings random walk example
source("3_metropolis_hastings.R", echo = FALSE)

# converged chain
set.seed(61)
post0 = mh(n=n, ybar=ybar, n_iter=10e3, mu_init=0.0, cand_sd=0.9)
coda::traceplot(as.mcmc(post0$mu[-c(1:500)]))

# wandering chain
set.seed(61)
post1 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post1$mu[-c(1:500)]))

" A trend can be observed, so the number of iterations needs to be increased to hit stationary distribution. "
set.seed(61)
post2 = mh(n=n, ybar=ybar, n_iter=100e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post2$mu))

# autocorrelation

# chain 1
coda::autocorr.plot(as.mcmc(post0$mu))
coda::autocorr.diag(as.mcmc(post0$mu))

# chain 2
coda::autocorr.plot(as.mcmc(post1$mu))
coda::autocorr.diag(as.mcmc(post1$mu))

" Autocorrelation reduction is low in the second chain. "


# effective sample size

str(post2)
coda::effectiveSize(as.mcmc(post2$mu)) 
" Equivalent to a sample of 374 after 100k iterations. "

# autocorrelation for chain 3
coda::autocorr.plot(as.mcmc(post2$mu), lag.max=500)

" Autocorrelation is 0 after every ~400 iterations. "

# Pick and store every 400th element in Markov Chain
thin_interval = 400 # how far apart the iterations are for autocorrelation to be essentially 0.
thin_indx = seq(from=thin_interval, to=length(post2$mu), by=thin_interval)
head(thin_indx)

post2mu_thin = post2$mu[thin_indx]

# compared filtered Markov chain with original 
traceplot(as.mcmc(post2$mu))
traceplot(as.mcmc(post2mu_thin))

# autocorrelation for filtered chain
coda::autocorr.plot(as.mcmc(post2mu_thin), lag.max=10)

# effective sample size vs number of iterations for filtered chain
effectiveSize(as.mcmc(post2mu_thin))
length(post2mu_thin)

# effective sample size vs number of iterations for converged chain
str(post0) 
coda::effectiveSize(as.mcmc(post0$mu)) 

?effectiveSize
" ~2500 independent samples for every 10k iterations. "

# number of samples required for confidence intervals

# Raftery and Lewis diagnostic
raftery.diag(as.mcmc(post0$mu))

" Total 13218 samples and effective sample size of 3746 required for reliable 95% interval. "

# sample size for 99% posterior interval with 0.001 margin of error
raftery.diag(as.mcmc(post0$mu), q=0.005, r=0.001, s=0.95)

?raftery.diag

" Minimum 19112 samples required for reliable 99% interval. "



# burn-in
set.seed(62)
post3 = mh(n=n, ybar=ybar, n_iter=500, mu_init=10.0, cand_sd=0.3)

coda::traceplot(as.mcmc(post3$mu))

" First ~100 iterations are not from a stationary distribution and should be discarded. "

# simulate multiple chains with different starting states
set.seed(61)

nsim = 500
post1 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=15.0, cand_sd=0.4)
post1$accpt

post2 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-5.0, cand_sd=0.4)
post2$accpt

post3 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=7.0, cand_sd=0.1)
post3$accpt

post4 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=23.0, cand_sd=0.5)
post4$accpt

post5 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-17.0, cand_sd=0.4)
post5$accpt

# store MCMC chains in a single list
pmc = mcmc.list(as.mcmc(post1$mu), as.mcmc(post2$mu), 
                as.mcmc(post3$mu), as.mcmc(post4$mu), as.mcmc(post5$mu))
str(pmc)

# plot the 5 chains
coda::traceplot(pmc)

" All chains appear to have converged to stationary distribution after ~200 iterations. "

# variability between chains

# Gelman and Rubin diagnostic
coda::gelman.diag(pmc)

" Scale reduction factor is almost 1, so there is no variability and all of them have hit the stationary distribution. "

# plot iterations vs shrink factor
coda::gelman.plot(pmc)

" The chains seem to converge after 300 iterations. "

# Monte Carlo estimation

" Once the Markov chain has converged, it can be used as a Monte Carlo sample from the posterior distribution. "

nburn = 1000 # discard early iterations
post0$mu_keep = post0$mu[-c(1:1000)]
summary(as.mcmc(post0$mu_keep))

mean(post$mu_keep > 1.0) # posterior probability that mu  > 1.0
