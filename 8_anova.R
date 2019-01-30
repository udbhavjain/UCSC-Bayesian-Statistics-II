# load plant growth dataset
data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)

# plot distribution for each group
boxplot(weight ~ group, data=PlantGrowth)

# reference analysis with linear model
lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)

anova(lmod)

" p-value is small i.e at least one of the means is significantly different from the others. "

# fit model in JAGS
library("rjags")

mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*1.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
  inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

# check for convergence
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)

" The chains appear to have converged, shrink factor is 1, autocorrelation is almost 0,
  and effective sample sizes are close to number of iterations for all groups. "

# calculate residuals
(pm_params = colMeans(mod_csim))
yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat

# residuals vs data index
plot(resid)

# residuals vs predicted values
plot(yhat, resid)

" Residual variance varies between groups. "


# posterior summary of parameters
summary(mod_sim)

# highest posterior density intervals
HPDinterval(mod_csim)

# probability of treatment 2 increasing mean yield
mean(mod_csim[,3] > mod_csim[,1])

# probability of treatment 2 increasing mean yield by at least 10%
mean(mod_csim[,3] > 1.1*mod_csim[,1])

" Approximately 50/50 odds of at least 10% improvement. "


# model with different variance for each group

mod2_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
sig[j] = sqrt( 1.0 / prec[j] )
}


} "

set.seed(32)
str(PlantGrowth)
data2_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params2 = c("mu", "sig")

inits2 = function() {
  inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                       variable.names=params2,
                       n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combined chains

# check for convergence
plot(mod2_sim)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
effectiveSize(mod2_sim)

" The chains appear to have converged, shrink factor is 1, autocorrelation is almost 0,
  and effective sample sizes are close to number of iterations for all groups. "

# calculate residuals
(pm2_params = colMeans(mod2_csim))
y2hat = pm2_params[1:3][data2_jags$grp]
resid2 = data2_jags$y - y2hat

# residuals vs data index
plot(resid2)

# residuals vs predicted values
plot(y2hat, resid2) # model 2
plot(yhat, resid)   # model 1
