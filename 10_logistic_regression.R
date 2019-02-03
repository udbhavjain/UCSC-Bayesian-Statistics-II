library("boot")

# load urine dataset
data("urine")

?urine
head(urine)

dat = na.omit(urine)

# pairwise scatterplots
pairs(dat)

# check correlation
#install.packages("corrplot")
library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

# variable selection

" To select the best variables, priors with values near 0 should be used for the
  beta-coefficients. The association between explanatory and response variables will 
  be established by the data. "

# scale the covariates to be around 0
X = scale(dat[,-1], center=TRUE, scale=TRUE)

head(X[,"gravity"])
colMeans(X)

# standard deviation for explanatory variables after standardisation
apply(X, 2, sd)

# model

" Priors for the betas will be drawn from a Laplace/double exponential distribution. "

# Laplace distribution example
ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

# define model for JAGS
library("rjags")

" Logit of probability is calculated using a regression model. 
  
  Log(p/(1-p)) = b0 + b1x1 + b2x2 + ... + bnxn
  
  Response variable is drawn from a Bernoulli trial.

"

mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

" Chains appear to have converged, autocorrelation is close to 0, shrink factor is almost 1,
  and effective sample sizes are good, except for 'osmo'. "

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)

# summary of model
summary(mod1_sim)

# plot marginal posterior probabilities for beta coefficients
par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))

" We will drop the variables that have posterior distributions close to 0. 
  The variable 'urea' is a borderline case, but it is correlated with 'gravity',
  so it will be dropped. "


# Model 2
mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)

# compare models
dic1
dic2

# summary of second model
summary(mod2_sim)

" 'specific gravity' and 'calcium concentration' are associated with increased probability
  of observing crystals. 'conductivity' is inversely related to probability of crystals. "

# 95% highest posterior density interval
HPDinterval(mod2_csim)

# plot marginal posterior probabilities of beta coefficients
par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))


# prediction

# posterior means
(pm_coef = colMeans(mod2_csim))

" log(p/1-p) = b1x1 + b2x2 + b3x3 + int

  Substitute Z for b1x1 + b2x2 + b3x3 + int 
  
  p/1-p = e^Z
  p = e^Z - p.e^Z
  p(1 + e^Z) = e^Z
  p = e^Z/(1+e^Z)

  dividing numerator and denominator on rhs of equation with e^Z:
  p = 1/(1 + e^(-Z))

  Replacing Z with original value:

  p = 1/(1 + e^-(b1x1 + b2x2 + b3x3 + int))
  
  This is called the link function.
  
  As all covariates were normalised, 0 is the average value for them.
  Therefore, at average values of covariates, the probability is:

  p = 1/(1 + e^-(-0.15)) ; value of intercept taken from simulated distribution.

  If we want to make a prediction for a new specimen whose value of gravity is average, 
  cond is one standard deviation below the mean, calc is one standard deviation above 
  the mean, the point estimate for the probability of calcium oxalate crystals is:

  p = 1/(1 + e^-(-0.15 + 1.4*0 -1.3*(-1) + 1.9*1))

"

# predicted probability for each data point

# multiply explanatory variables with beta coefficients
pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]

# calculate probability using link function
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

# compare predicted values against actual data
plot(phat, jitter(dat$r))

# set cutoff as >0.5 = 1 and check correctness of classifications
(tab0.5 = table(phat > 0.5, data_jags$y))

# accuracy of model = (true positives + true negatives) / total observations 
sum(diag(tab0.5)) / sum(tab0.5)

" Suppose it is considered really bad to predict no calcium oxalate crystal when 
  there in fact is one (i.e false negatives are bad). We might then choose to lower 
  the threshold for classifying data points as 1s. "

# set cutoff as >0.3 = 1 and check correctness of classifications
(tab0.3 = table(phat > 0.3, data_jags$y))

# accuracy of model
sum(diag(tab0.3)) / sum(tab0.3)
