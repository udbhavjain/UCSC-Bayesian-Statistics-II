#install.packages("car")
library("car")

# load Leinhardt dataset from 'car' package
data("Leinhardt")
?Leinhardt

head(Leinhardt)
str(Leinhardt)

pairs(Leinhardt)

# correlation between infant mortality and per capita income
plot(infant ~ income, data=Leinhardt)

" There does not seem to be a linear relationship. "

# distribution of infant mortality
hist(Leinhardt$infant)

# distribution of per capita income
hist(Leinhardt$income)

" Both are highly right skewed. "

# correlation between log of infant mortality and log of income
Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)

plot(loginfant ~ logincome, data=Leinhardt)

" Negative, linear relationship. "

# linear model of the logs
lmod = lm(loginfant ~ logincome, data=Leinhardt)
summary(lmod)


# Fit model in JAGS

" Equation for linear regression:
  y = B0 + B1X1 + ... + BnXn + Ei; Ei iid~ N(0, σ^2) for i=1,...,n 
  
  Where B0 + B1X1 + ... + BnXn can be considered the mean for y.

  Model:
  
  Likelihood:
  yi iid~ N(μi,σ^2)
  μi ind~ b1 + b2.log(income)

  Prior:
  bj ~ N(0,1) ; j = 1,2
  σ^2 ~ inv-gamma(5/2,25)

"

# remove NA values
dat = na.omit(Leinhardt)

library("rjags")

# specify model
mod1_string = " model {
    for (i in 1:n) {
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = b[1] + b[2]*log_income[i] 
    }
    
    for (i in 1:2) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*10.0/2.0)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

# data for JAGS
set.seed(72)
data1_jags = list(y=dat$loginfant, n=nrow(dat), 
                  log_income=dat$logincome)

# parameters to be monitored - combined betas and standard deviation
params1 = c("b", "sig")

# initial values
inits1 = function() {
  inits = list("b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

# create model with specifications above and 3 chains
mod1 = jags.model(textConnection(mod1_string), data=data1_jags, inits=inits1, n.chains=3)

# run model
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains


# Convergence testing
plot(mod1_sim)

# shrink factors
gelman.diag(mod1_sim)

# autocorrelation
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)

" High autocorrelation for betas. "

# check effective sample sizes
effectiveSize(mod1_sim)

" Very small effective sample sizes for betas. "

# posterior summary of parameters
summary(mod1_sim)


# residual checks

# un-transformed variables

# check for independence
lmod0 = lm(infant ~ income, data=Leinhardt)
plot(resid(lmod0))

" Seems independent. "

# check for linearity and constant variance
plot(predict(lmod0), resid(lmod0))

" Downward trend visible, variance not constant at higher predictions. "

# check for normality
qqnorm(resid(lmod0))

" Curved plot. "

# log-transformed variables

# log-transformed income 
X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)

# posterior means of simulated samples
(pm_params1 = colMeans(mod1_csim))

# multiply values of income with beta 2 and add beta 1
" y = b[1] + b[2].x "

yhat1 = drop(X %*% pm_params1[1:2])

# calculate residuals 
resid1 = data1_jags$y - yhat1

# plot residuals against data index
plot(resid1) 
" No patterns. "

# plot residuals against predicted values
plot(yhat1, resid1)
" Centred around 0, but variance seems to increase for higher predicted values. A couple of outliers can be seen. "

# normality check for residuals
qqnorm(resid1)

" Linear, except for the outliers. "

# predicted values vs residuals for reference linear model
plot(predict(lmod), resid(lmod)) 

" Identical to the simulated predictions vs residuals plot. "

# check the outliers
rownames(dat)[order(resid1, decreasing=TRUE)[1:5]] 


# check for additional covariates

# Add oil export variable to the model
mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
    }
    
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*10.0/2.0)
    sig = sqrt( 1.0 / prec )
} "


set.seed(73)
data2_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"))
data2_jags$is_oil

params2 = c("b", "sig")

inits2 = function() {
  inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
update(mod2, 1e3) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combine multiple chains

# convergence diagnostics

# view chains
plot(mod2_sim)

# shrink factor
gelman.diag(mod2_sim)

# autocorrelation
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)

# effective sample size
effectiveSize(mod2_sim)

# posterior summary
summary(mod2_sim)

" Oil production and infant mortality are positively correlated. "

# calculate residuals
X2 = cbind(rep(1.0, data1_jags$n), data2_jags$log_income, data2_jags$is_oil)
head(X2)

(pm_params2 = colMeans(mod2_csim)) # posterior mean

yhat2 = drop(X2 %*% pm_params2[1:3])
resid2 = data2_jags$y - yhat2

# residuals vs data index
plot(resid2)

# residuals vs predicted values 
plot(yhat2, resid2) # model 2
plot(yhat1, resid1) # model 1

# standard deviation of residuals
sd(resid2)

" Although deviation is lower, Saudi Arabia and Libya are still outliers. "


# model with t-distribution likelihood for y instead of normal

mod3_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dt( mu[i], tau, df )
        mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
    }
    
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
    nu ~ dexp(1.0)
    
    tau ~ dgamma(5/2.0, 5*10.0/2.0) # tau is close to, but not equal to the precision
    sig = sqrt( 1.0 / tau * df / (df - 2.0) ) # standard deviation of errors
} "


set.seed(28)
data3_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"))

params3 = c("b", "sig")

inits3 = function() {
  inits = list("b"=rnorm(3,0.0,100.0), "tau"=rgamma(1,1.0,1.0))
}

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, inits = inits3, n.chains=3)
update(mod3, 1e3) # burn-in

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)

mod3_csim = as.mcmc(do.call(rbind, mod3_sim)) # combine multiple chains

# convergence diagnostics

# view chains
plot(mod3_sim)

# shrink factor
gelman.diag(mod3_sim)

# autocorrelation
autocorr.diag(mod3_sim)
autocorr.plot(mod3_sim)

# effective sample size
effectiveSize(mod3_sim)

# posterior summary
summary(mod3_sim)

# calculate residuals
X3 = cbind(rep(1.0, data1_jags$n), data3_jags$log_income, data3_jags$is_oil)
head(X3)

(pm_params3 = colMeans(mod3_csim)) # posterior mean

yhat3 = drop(X3 %*% pm_params3[1:3])
resid3 = data3_jags$y - yhat3

# residuals vs data index
plot(resid3)

# residuals vs predicted values 
plot(yhat3, resid3) # model 3
plot(yhat2, resid2) # model 2

# standard deviation of residuals
sd(resid3)


# compare models

# deviance information criterion
dic.samples(mod1, n.iter=1e3)
dic.samples(mod2, n.iter=1e3)
dic.samples(mod3, n.iter=1e3)

" Model 2 has the lowest penalised deviance. "