# load infant mortality dataset
library("car")
data("Leinhardt")

?Leinhardt
str(Leinhardt)

pairs(Leinhardt)
head(Leinhardt)

dat = na.omit(Leinhardt)
dat$logincome = log(dat$income)
dat$loginfant = log(dat$infant)
str(dat)

# JAGS model for linear regression, with different intercepts for different regions

library("rjags")

mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = a[region[i]] + b[1]*log_income[i] + b[2]*is_oil[i]
}

for (j in 1:max(region)) {
a[j] ~ dnorm(a0, prec_a)
}

a0 ~ dnorm(0.0, 1.0/1.0e6)
prec_a ~ dgamma(1/2.0, 1*10.0/2.0)
tau = sqrt( 1.0 / prec_a )

for (j in 1:2) {
b[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(116)
data_jags = list(y=dat$loginfant, log_income=dat$logincome,
                 is_oil=as.numeric(dat$oil=="yes"), region=as.numeric(dat$region))
data_jags$is_oil
table(data_jags$is_oil, data_jags$region)

params = c("a0", "a", "b", "sig", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combine multiple chains

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

# DIC
dic.samples(mod, n.iter=1e3)

" Better than the non-hierarchical model. "

# posterior summary
summary(mod_sim)
