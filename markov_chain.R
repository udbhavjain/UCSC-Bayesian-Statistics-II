# Continuous Markov chain - random walk

" Transition model: p(Xt+1|Xt=xt)=N(xt,1)
  Probability distribution for the next state is Normal with variance 1 and mean equal to the current state.
  Markov chain, since transition to next state depends only on current state. "

# build a Markov chain of 100 elements
set.seed(34)

n = 100
x = numeric(n)

for (i in 2:n) {
  x[i] = rnorm(1, mean=x[i-1], sd=1.0)
}

x11()
plot.ts(x)

# transition matrix

" Experiment: Flip a coin.

    1. If the coin turns up heads, then increase your secret number by one (5 increases to 1).
    2. If the coin turns up tails, then decrease your secret number by one (1 decreases to 5).

  Transition matrix for this is defined below, where the rows represent the current state and the
  the numbers represent the probabilities of moving to the state equivalent to the column number. 
"

Q = matrix(c(0.0, 0.5, 0.0, 0.0, 0.5,
             0.5, 0.0, 0.5, 0.0, 0.0,
             0.0, 0.5, 0.0, 0.5, 0.0,
             0.0, 0.0, 0.5, 0.0, 0.5,
             0.5, 0.0, 0.0, 0.5, 0.0), 
           nrow=5, byrow=TRUE)

" p(Xt+1=5∣Xt=4) can be found in the fourth row, fifth column. "
" To see the probabilities of transition after two steps, the matrix is squared. "

Q %*% Q

" p(Xt+2=1∣Xt=4) is 0.25. "


# stationary distribution - discrete

# transition matrix for 5 steps ahead
Q5 = Q %*% Q %*% Q %*% Q %*% Q 
round(Q5, 3)

# transition matrix for 10 steps ahead
Q10 = Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q
round(Q10, 3)


# transition matrix for 30 steps ahead
Q30 = Q
for (i in 2:30) {
  Q30 = Q30 %*% Q
}
round(Q30, 3)

" Transition distributions appear to be converging. "

# perform number experiment 5000 times, with probabilities drawn from matrix Q

n = 5000
x = numeric(n)
x[1] = 1 # start from state 1

for (i in 2:n) 
{
  # draw any number from 1 to 5, based on probabilities selected from rows of Q, depending on current state
  x[i] = sample.int(5, size=1, prob=Q[x[i-1],]) 
}

# calculate probability of each state based on the experiment above
table(x) / n

" Approximately equal to stationary distribution. "



# stationary distribution - continuous
" Modified random walk so that it has a stationary distribution: 
   p(Xt+1|Xt=xt)=N(PHI*xt,1) where −1<PHI<1
"

# simulation for PHI = -0.6

set.seed(38)

n = 1500
x = numeric(n)
phi = -0.6

for (i in 2:n) {
  x[i] = rnorm(1, mean=phi*x[i-1], sd=1.0)
}

x11()
plot.ts(x)

" The theoretical stationary distribution for this chain is normal with mean 0 and variance 1/(1−phi^2), 
  which approximately equals 1.562 for phi = -0.6. "

# compare simulated and theoretical distributions
x11()
hist(x, freq=FALSE)
curve(dnorm(x, mean=0.0, sd=sqrt(1.0/(1.0-phi^2))), col="red", add=TRUE)
legend("topright", legend="theoretical stationary\ndistribution", col="red", lty=1, bty="n")

" The simulated distribution is similar to the theoretical. Therefore, the chain has reached stationary distribution
  and the simulation can be treated as a Monte Carlo sample. "