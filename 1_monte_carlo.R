set.seed(32) 

# number of draws
m = 100

# shape parameter for gamma
a = 2.0

# rate parameter for gamma
b = 1.0 / 3.0

# simulate m samples
theta = rgamma(n=m, shape=a, rate=b)

# plot the monte carlo sample distribution and the true gamma curve 
hist(theta, freq=FALSE)
curve(dgamma(x=x, shape=a, rate=b), col="blue", add=TRUE)

# expected value/mean of monte carlo sample
sum(theta)/m
mean(theta)

# true expected value/mean of gamma
a/b

# mean with 10k samples
m = 1e4
theta = rgamma(n=m, shape=a, rate=b)
mean(theta)

" Much closer to true value. "

# variance of sample
var(theta)

# true variance of gamma
a/b^2

# probability of theta < 5 in sample
ind = theta < 5.0
mean(ind)

# true probability of theta < 5
pgamma(q=5.0, shape=a, rate=b) 

# 90th quantile of sample
quantile(x=theta, probs=0.9)

# true 90th quantile
qgamma(p=0.9, shape=a, rate=b) 



# calculation of monte carlo error

# calculate standard error and margin of error for 95% confidence interval for estimate of mean
se = sd(theta) / sqrt(m)
2.0 * se

# margin of error for 95% confidence interval for estimate of p(theta < 0.5)
ind = theta < 5.0
se = sd(ind) / sqrt(m)
2.0 * se 



# simulation of hierarchical model

" 
 y|phi ~ Bin(10,phi)
 phi ~ Beta(2,2)
"

# number of simulations
m = 10e4

# vectors for storing phi and y
y = numeric(m) 
phi = numeric(m)

# simulate 100k samples for Beta(2,2)
phi = rbeta(n=m, shape1=2.0, shape2=2.0)

# simulate 100k draws for Bin(10,phi) using the simulated values for phi 
y = rbinom(n=m, size=10, prob=phi)

# marginal distribution for draws of y
mean(y)

plot(prop.table(table(y)), ylab="P(y)", main="Marginal distribution of y")
