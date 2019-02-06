#install.packages("COUNT")
library("COUNT")

# load badhealth dataset
data("badhealth")

?badhealth
head(badhealth)

# check NAs
any(is.na(badhealth))

# plot number of visits to doctor
hist(badhealth$numvisit, breaks=20)

# plot number of visits against age for bad and good health
plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")

# JAGS model with interaction term
"
Model:

likelihood:
  y ~ Pois(λi)
  
  where λi is defined by the link function:
  log(λi) = b0 + b1.x1 + b2.x2 + bintx.x1.x2
  
"

library("rjags")

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

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

" Autocorrelation is somewhat high and effective sample size is small. "

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

# model matrix
X = as.matrix(badhealth[,-1])
X = cbind(X, with(badhealth, badh*age))
head(X)

# medians of coefficients
(pmed_coef = apply(mod_csim, 2, median))

# calculate lambdas
llam_hat = pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)

hist(lam_hat)

# residuals
resid = badhealth$numvisit - lam_hat
plot(resid) # dataset was ordered

# lambda vs number of visits
plot(lam_hat, badhealth$numvisit)
abline(0.0, 1.0)

# lambda vs residuals
plot(lam_hat[which(badhealth$badh==0)], resid[which(badhealth$badh==0)], xlim=c(0, 8), ylab="residuals", xlab=expression(hat(lambda)), ylim=range(resid))
points(lam_hat[which(badhealth$badh==1)], resid[which(badhealth$badh==1)], col="red")

" Variability increases with lambda, as expected. "

# variances in data
var(resid[which(badhealth$badh==0)])
var(resid[which(badhealth$badh==1)])

" Residual variances are a lot higher than the lambdas. "

# posterior summary
summary(mod_sim)

" Age has a positive association with number of doctor visits. Bad health is associated 
  with an increase in expected number of visits. The interaction coefficient is interpreted 
  as an adjustment to the age coefficient for people in bad health. Hence, for people with 
  bad health, age is essentially unassociated with number of visits."

# predictive distributions

x1 = c(0, 35, 0) # good health
x2 = c(1, 35, 35) # bad health

# posterior samples
head(mod_csim)

# calculate log lambdas
loglam1 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x1
loglam2 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x2

# inverse link
lam1 = exp(loglam1)
lam2 = exp(loglam2)

# simulate number of visits using these lambdas
n_sim = length(lam1)

y1 = rpois(n=n_sim, lambda=lam1)
y2 = rpois(n=n_sim, lambda=lam2)

plot(table(factor(y1, levels=0:18))/n_sim, pch=2, ylab="posterior prob.", xlab="visits")
points(table(y2+0.1)/n_sim, col="red")

# probability that person with poor health will have more visits
mean(y2 > y1)
