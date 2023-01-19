#Lab #X: Bayesian Inference and Markov chain Monte Carlo Basics (Kelsey L. Ruckert, Tony E. Wong, Benjamin Seiyon Lee, Yawen Guan, and Murali Haran)
#Edited by Sitara Baboolal

#Copyright 2017 by the Authors

#This file is part of Risk Analysis in the Earth Sciences: A Lab Manual with Exercises in R. 

#This e-textbook is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. -->

#This e-textbook is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. -->

#You should have received a copy of the GNU General Public License along with this e-textbook.  If not, see <http://www.gnu.org/licenses/>. -->

###### This script folows the workflow outlined in MCMC Chapter 1 #######

# Set up the known probabilities based on prior knowledge.
# Define the transition matrix.
# Rows represent the current state.
# Columns are the next state.
#       FOGGY   SUNNY  RAINY
FOGGY <- c( 0.5, 0.25,  0.25)  
SUNNY <- c( 0.5,    0,   0.5)
RAINY <- c(0.25, 0.25,   0.5)
# Define the parameter space.
parameter <- c("Foggy", "Sunny", "Rainy")
P <- matrix(c(FOGGY, SUNNY, RAINY), nrow=3, ncol=3, byrow = TRUE, 
            dimnames = list(parameter, parameter))
print(P)

# PREDICTING THE WEATHER WITH A MARKOV CHAIN EXAMPLE
#================================================
# Set the prediction length and the initial parameter value.
Prediction_Days <- 25
CurrentStateProb <- mat.or.vec(Prediction_Days, length(P[1,]))
CurrentStateProb[1,] <- c(1, 0, 0) # Today is foggy
# Run the Markov chain.
for(i in 2:Prediction_Days){
  # Current parameter value times probability
  Prob <- CurrentStateProb[i-1, ] * P
  CurrentStateProb[i, ] <- c(sum(Prob[,1]), sum(Prob[,2]), sum(Prob[,3]))
}
colnames(CurrentStateProb) <- parameter
print(CurrentStateProb)
# Print weather predictions.
print(paste("p(", parameter, ") in 1 day = ", round(CurrentStateProb[2, ], 2)))
print(paste("p(", parameter, ") in 5 days = ", round(CurrentStateProb[6, ], 2)))
print(paste("p(", parameter, ") in 24 days = ", round(CurrentStateProb[25, ], 2)))


# SAMPLING THE PROBABILITY DISTRIBUTION; A MARKOV CHAIN MONTE CARLO EXAMPLE
#================================================
# Set the initial parameter value probability.
CurrentState <- sample(parameter, 1, prob = c(0.5, 0.25, 0.25)) # Today is foggy
# Start sampling (iteratively) from the distribution.
weather <- c()
for (i in 1:1e4){
  NextState <- sample(parameter, 1, prob = P[CurrentState,])
  weather[i] <- NextState
  CurrentState <- NextState
}
# Throw away the first 1% of the data (Burn-in)
burnin <- seq(from = 1, to = 1e4*0.1, by = 1)
weatherDraws <- weather[-burnin]
# DISPLAY 
#================================================
par(mfrow = c(1,2))
# Display the results from the Markov Chain
plot(1:Prediction_Days, CurrentStateProb[ ,1], xlab = "Days", ylab = "P(Weather)", 
     ylim = c(0,1), type = "l", lwd = 2, col = "gray")
lines(1:Prediction_Days, CurrentStateProb[ ,2], lwd = 2, col = "gold")
lines(1:Prediction_Days, CurrentStateProb[ ,3], lwd = 2, col = "blue")
legend("right", c("Foggy", "Sunny", "Rainy"), lwd = 2, bty = "n", col = c("gray", "gold", 
                                                                          "blue"))
# Display the results from MCMC
barplot(prop.table(table(weatherDraws)), ylim = c(0,1), 
        sub = "Equilibrium distribution\nof the weather", col = c("gray", "blue", "gold"))
abline(h = 0.4, lty = 2); abline(h = 0.2, lty = 2)

###Tutorial
# Clear away any existing variables or figures and set the seed for random sampling.
rm(list = ls())
graphics.off()
set.seed(1234)
# Define the true model parameters and set up the time scale.
alpha.true <- 2 # Arbitrary choices.
beta.true <- -5
time <- 1:10
# Generate some observations with noise (measurement error) and the model.
y.true = alpha.true*time + beta.true
sigma = 1
meas_err <- rnorm(length(time), mean = 0, sd = sigma)
y = y.true + meas_err
# Plot the true model and observations.
par(mfrow = c(1,1))
plot(time, y.true, type="l", main="True observations and model",
     xlab = "Time", ylab = "Observations")
points(time, y, pch = 20)

####Define the log posterior
plot(0:4, 0:4, type="n", xlab="x", ylab="f(x)", xaxt = "n", yaxt = "n")
segments(1, -1, 1, 2.5)
segments(3, -1, 3, 2.5)
segments(1, 2.5, 3, 2.5, lty=2)
axis(1, c(0, 1,3), labels=c(0,"a", "b"))
axis(2, 2.5, labels=expression(frac(1,b-a)), las=1)

# Define the unnormalized log posterior.
logp = function(theta){
  N = length(time)
  
  # Calulate model simulations.
  alpha = theta[1]; beta = theta[2]
  model = alpha*time + beta
  
  # Estimate the residuals (i.e., the deviation of the observations from the 
  # model simulation).
  resid = 	y - model
  
  # Get the log of the likelihood function.
  log.likelihood = -N/2*log(2*pi) - N*log(sigma) - 1/2*sum(resid^2)/sigma^2
  
  # Use an improper uniform "relatively uninformative" prior.
  log.prior = 0 # log(1)
  
  # Bayesian updating: update the probability estimates based on the observations.
  log.posterior = log.likelihood + log.prior
  
  # Return the unnormalized log posterior value. 
  return(log.posterior)
}

#####Constructing a Markov Chain
# Set the number of MCMC iterations.
NI = 30000
# Start with some initial state of parameter estimates, theta^initial
alpha.init = runif(1, -50, 50) # Arbitrary choices.
beta.init = runif(1, -50, 50)  
theta = c(alpha.init, beta.init)
```

Using the initial parameter values $\theta^{initial}$, you then evaluate the unnormalized posterior of that value, $p(\theta^{initial} \mid y)$, using the likelihood function ($L(y \mid \theta)$) and the prior distribution.
```{r eval=FALSE}
# Evaluate the unnormalized posterior of the parameter values
# P(theta^initial | y)
lp = logp(theta)

# Setup some variables and arrays to keep track of:
theta.best = theta         # the best parameter estimates
lp.max = lp                # the maximum of the log posterior
theta.new = rep(NA,2)      # proposed new parameters (theta^new)
accepts = 0                # how many times the proposed new parameters are accepted
mcmc.chains = array(dim=c(NI,2)) # and a chain of accepted parameters
# Set the step size for the MCMC.
step = c(0.1, 1)
# Metropolis-Hastings algorithm MCMC; the proposal distribution proposes the next 
# point to which the random walk might move. For this algorithm, this proposal  
# distribution is symmetric, that is P(x to x`) = P(x` to x).
for(i in 1:NI) {
  # Propose a new state (theta^new) based on the current parameter values 
  # theta and the transition probability / step size
  theta.new = rnorm(2, theta, sd = step)
  
  # Evaluate the new unnormalized posterior of the parameter values
  # and compare the proposed value to the current state
  lp.new = logp(theta.new)
  lq = lp.new - lp
  # Metropolis test; compute the acceptance ratio
  # Draw some uniformly distributed random number 'lr' from [0,1]; 
  lr = log(runif(1))
  
  # If lr < the new proposed value, then accept the parameters setting the 
  # proposed new theta (theta^new) to the current state (theta).
  if(lr < lq) {
    # Update the current theta and log posterior to the new state.
    theta = theta.new
    lp = lp.new
    
    # If this proposed new parameter value is “better” (higher posterior probability)
    # than the current state, then accept with a probability of 1. Hence, increase
    # the number of acceptions by 1.
    accepts = accepts + 1
    
    # Check if the current state is the best, thus far and save if so.
    if(lp > lp.max) {
      theta.best = theta
      lp.max = lp	
    }
  }
  # Append the parameter estimate to the chain. This will generate a series of parameter
  # values (theta_0, theta_1, ...). 
  mcmc.chains[i,] = theta
}


#####Checking the acceptance rate
# Calculate the parameter acceptance rate; it should be higher than 23.4%.
accept.rate <- (accepts/NI) * 100
print(accept.rate)



