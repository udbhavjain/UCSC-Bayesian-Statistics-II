# load cookies dataset
dat = read.table(file="cookies.dat", header=TRUE)
head(dat)

table(dat$location)

# distribution of chips
hist(dat$chips)

# distribution by location
boxplot(chips ~ location, data=dat)

# prior predictive checks

" Priors for λ, are drawn from a gamma distribution, the hyperparameters for which
  are alpha and beta. Every location has its own λ, and the alpha and beta values
  control the distribution λ between locations. Mean of λ distribution is the mean
  of the overall number of chips. "

# simulate alpha and beta
set.seed(112)
n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri
sig_pri = sqrt(alpha_pri/beta_pri^2)

summary(mu_pri)
summary(sig_pri)

# simulate lambda
lam_pri = rgamma(n=n_sim, shape=alpha_pri, rate=beta_pri)
summary(lam_pri)

# prior predictive reconstruction of dataset
(lam_pri = rgamma(n=5, shape=alpha_pri[1:5], rate=beta_pri[1:5]))
(y_pri = rpois(n=150, lambda=rep(lam_pri, each=30)))


# JAGS model
library("rjags")

"
  Model:
  yi|li,λi ind~ Pois(λli) ; li = {1,2,3,4,5}, i = {1,...,150}
  λl|α,β iid~ Gamma(α,β) ; l = {1,2,3,4,5}
  
  α ~ p(α)
  β ~ p(β)
"

mod_string = " model {
for (i in 1:length(chips)) {
chips[i] ~ dpois(lam[location[i]])
}

for (j in 1:max(location)) {
lam[j] ~ dgamma(alpha, beta)
}

alpha = mu^2 / sig^2
beta = mu / sig^2

mu ~ dgamma(2.0, 1.0/5.0)
sig ~ dexp(1.0)

} "

set.seed(113)

data_jags = as.list(dat)

params = c("lam", "mu", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

# residuals
(pm_params = colMeans(mod_csim))

yhat = rep(pm_params[1:5], each=30)
resid = dat$chips - yhat

# vs index
plot(resid)

# vs predicted value
plot(jitter(yhat), resid)

# residual variance for different values of y
var(resid[yhat<7])
var(resid[yhat>11])

# vs lambda mean 
lam_resid = pm_params[1:5] - pm_params["mu"]
plot(lam_resid)
abline(h=0, lty=2)


# model summary
summary(mod_sim)


# posterior predictive simulation

# use draws for mean and sigma to simulate new values
(n_sim = nrow(mod_csim))

lam_pred = rgamma(n=n_sim, shape=mod_csim[,"mu"]^2/mod_csim[,"sig"]^2, 
                  rate=mod_csim[,"mu"]/mod_csim[,"sig"]^2)
hist(lam_pred)

mean(lam_pred > 15)

# simulate number of chips per cookie using draws for lambda
y_pred = rpois(n=n_sim, lambda=lam_pred)
hist(y_pred)

mean(y_pred > 15)
hist(dat$chips)

# probability that a cookie from location 1 will have <7 chips
y_pred1 = rpois(n=n_sim, lambda=mod_csim[,"lam[1]"])
hist(y_pred1)

mean(y_pred1 < 7)
