library("COUNT")
library("rjags")

data("badhealth")

# run model and create posterior distribution

mod_string = " model {
for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
}

int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/1e4)
b_age ~ dnorm(0.0, 1.0/1e4)
b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

# posterior density for bad health coefficient
plot(density(mod_csim[,"b_badh"]))

# using near-zero prior for bad health and interaction coefficients
mod2_string = " model {
    for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
}

int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/0.2^2)
b_age ~ dnorm(0.0, 1.0/1e4)
b_intx ~ dnorm(0.0, 1.0/0.01^2)
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

# compare posterior densities for bad health coefficient
curve(dnorm(x, mean=0.0, sd=sqrt(1e4)), from=-3.0, to=3.0, ylim=c(0.0, 3.0), lty=2,
      main="b_badh", ylab="density", xlab="b_badh")
curve(dnorm(x, mean=0.0, sd=0.2), from=-3.0, to=3.0, col="red", lty=2, add=TRUE)
lines(density(mod_csim[,"b_badh"]))
lines(density(mod2_csim[,"b_badh"]), col="red")
legend("topleft", legend=c("noninformative prior", "posterior", "skeptical prior", "posterior"),
       lty=c(2,1,2,1), col=rep(c("black", "red"), each=2), bty="n")

" Although the posterior for bad health coefficient has dropped, it is still above 0. "

# compare posterior densities for interaction coefficient
curve(dnorm(x, mean=0.0, sd=sqrt(1e4)), from=-0.05, to=0.05, ylim=c(0.0, 140.0), lty=2,
      main="b_intx", ylab="density", xlab="b_intx")
curve(dnorm(x, mean=0.0, sd=0.01), from=-0.05, to=0.05, col="red", lty=2, add=TRUE)
lines(density(mod_csim[,"b_intx"]))
lines(density(mod2_csim[,"b_intx"]), col="red")
legend("topleft", legend=c("noninformative prior", "posterior", "skeptical prior", "posterior"),
       lty=c(2,1,2,1), col=rep(c("black", "red"), each=2), bty="n")

mean(mod2_csim[,"b_intx"] > 0)

" Using a skeptical prior makes the interaction coefficient positive. "
