# load warpbreaks dataset
data("warpbreaks")

?warpbreaks
head(warpbreaks)
table(warpbreaks$wool, warpbreaks$tension)

# plot distribution of breaks
boxplot(breaks ~ wool + tension, data=warpbreaks)

# plot distribution of log of breaks
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)

# one-way model using tension variable
library("rjags")

mod1_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[tensGrp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*2.0/2.0)
sig = sqrt(1.0 / prec)
} "

set.seed(83)
str(warpbreaks)

data1_jags = list(y=log(warpbreaks$breaks), tensGrp=as.numeric(warpbreaks$tension))

params1 = c("mu", "sig")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)

## convergence diagnostics
plot(mod1_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

" Chains converged, autocorrelation is ~0, shrink factor is 1, 
  effective sample sizes are close to number of iterations. "

summary(mod1_sim)

" 95% posterior interval for medium tension overlaps with both high and low tension intervals.
  High and low tension intervals barely overlap with each other. "

# DIC 
(dic1 = dic.samples(mod1, n.iter=1e3))


# two-way additive model

# design matrix
X = model.matrix( ~ wool + tension, data=warpbreaks)
head(X)
tail(X)

" Wool A and low tension are intercepts in reference model. "

# model in JAGS

mod2_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = int + alpha*isWoolB[i] + beta[1]*isTensionM[i] + beta[2]*isTensionH[i]
}

int ~ dnorm(0.0, 1.0/1.0e6)
alpha ~ dnorm(0.0, 1.0/1.0e6)
for (j in 1:2) {
beta[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(3/2.0, 3*1.0/2.0)
sig = sqrt(1.0 / prec)
} "

data2_jags = list(y=log(warpbreaks$breaks), isWoolB=X[,"woolB"], 
                  isTensionM=X[,"tensionM"], isTensionH=X[,"tensionH"])

params2 = c("int", "alpha", "beta", "sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

## convergene diagnostics
plot(mod2_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

summary(mod2_sim)

# DIC
(dic2 = dic.samples(mod2, n.iter=1e3))

# compare models
dic2
dic1

" Model 2 is better. "

# visualise two-way model
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)

" Wool B has a negative effect on number of breaks, according to model summary.
  In the box plot, wool B has lower number of breakages for low and high tension, 
  but higher for medium tension. Therefore, the effect is not consistent across 
  tension levels and an interaction term should be added."

# reference model with interacton terms
lmod2 = lm(log(breaks) ~ .^2, data=warpbreaks)
summary(lmod2)

" We now have the effect of wool B with medium and high tension. "

# JAGS model
mod3_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec)
    }
    
    for (j in 1:max(woolGrp)) {
        for (k in 1:max(tensGrp)) {
            mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
        }
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

str(warpbreaks)

data3_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params3 = c("mu", "sig")

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)
raftery.diag(mod3_sim)

# DIC
(dic3 = dic.samples(mod3, n.iter=1e3))

# DIC of older models
dic2
dic1

" Model that includes interaction between wool type and tension is the best. "

# model summary
summary(mod3_sim)

# 95% HPD interval
HPDinterval(mod3_csim)

# posterior density plots
par(mfrow=c(3,2)) # arrange frame for plots
densplot(mod3_csim[,1:6], xlim=c(2.0, 4.5))

# group with smallest number of breaks
prop.table( table( apply(mod3_csim[,1:6], 1, which.min) ) )

" Wool B with high tension. "