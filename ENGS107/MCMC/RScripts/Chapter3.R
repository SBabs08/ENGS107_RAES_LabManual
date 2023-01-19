#Lab #X: Applying Markov chain Monte Carlo to sea-level data (Kelsey L. Ruckert, Tony E. Wong, Yawen Guan, and Murali Haran)
#Edited by Sitara Baboolal 

#Copyright 2017 by the Authors

#This file is part of Risk Analysis in the Earth Sciences: A Lab Manual with Exercises in R. 

#This e-textbook is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. -->

#This e-textbook is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. -->

#You should have received a copy of the GNU General Public License along with this e-textbook.  If not, see <http://www.gnu.org/licenses/>. -->

###### This script follows the workflow outlined in MCMC Chapter 3 #######

#Estimate the log likelihood of the AR1 process
logl.ar1 <- function(r, sigma1, rho1, eps1 = y.meas.err){
  n <- length(r) # r is the residuals
  logl <- 0
  if(n>1) {
    # Whiten the residuals (w)
    w <- r[2:n] - rho1*r[1:(n-1)]
    logl <- logl + sum(dnorm(w, sd=sqrt((sigma1)^2 + (eps1[c(-1)])^2), log=TRUE))
  }
  return(logl)
}

## Tutorial
# Clear away any existing variables or figures.  
rm(list = ls())
graphics.off()
# Install and read in packages.
# install.packages("adaptMCMC")
library(adaptMCMC)
library(compiler)
enableJIT(3)
# Set the seed.
set.seed(111)
# Make directory for this lab
dir.create("lab_X")
dir.path <- paste(getwd(), "/lab_X/", sep="")
scrim.git.path <- "https://raw.githubusercontent.com/scrim-network/"
# Global mean sea level from Church and White (2011)
# year, time in years from 1880 to 2013 at the half year mark
# slr, global mean sea level in mm
# err.obs, global mean sea level measurement error in mm
url <- paste(scrim.git.path, "BRICK/BRICKms/data/GMSL_ChurchWhite2011_yr_2015.txt", sep="")
download.file(url, paste(dir.path, "GMSL_ChurchWhite2011_yr_2015.txt", sep=""))
church.data <- read.table("lab_X/GMSL_ChurchWhite2011_yr_2015.txt")
year <- church.data[ ,1] - 0.5 
slr <- church.data[ ,2]
err.obs <- church.data[ ,3]
# Calculate sea level -/+ observation errors.
err_pos <- slr + err.obs
err_neg <- slr - err.obs
# Read in the information from the Smith (2008) temperature file.
# project.years, time in years from 1880 to 2100.
# hist.temp, historical global mean temperatures from 1880 to 2013 in C
# rcp85, merged historical + rcp 8.5 temperatures from 1880 to 2100 in C 
url <- paste(scrim.git.path, 
             "Ruckertetal_SLR2016/master/RFILES/Data/NOAA_IPCC_RCPtempsscenarios.csv", 
             sep="")
download.file(url, paste(dir.path, "NOAA_IPCC_RCPtempsscenarios.csv", sep=""))
temp.data <- read.csv("lab_X/NOAA_IPCC_RCPtempsscenarios.csv")
project.years <- temp.data[1:match(2100, temp.data[,1]), 1]
hist.temp <- temp.data[1:match(2013, temp.data[,1]), 2]
rcp85 <- temp.data[1:match(2100, temp.data[,1]), 4]
# The sea level and temperature data have a yearly timestep.
tstep = 1
# Load the sea-level model (Rahmstorf, 2007)
url <- paste(scrim.git.path, "BRICK/BRICKms/R/gmsl_r07.R", sep="")
download.file(url, paste(dir.path, "gmsl_r07.R", sep=""))
source("lab_X/gmsl_r07.R")

### Characterizing the error structure
par(mfrow = c(1,1))
plot(year, err.obs, type = "l", ylab = "Observation errors [mm]", xlab = "Time")
points(year, err.obs, pch = 20)

# Generate fit to sea-level observations using General-purpose Optimization.
# by calculating the root mean squared error.
fn <- function(pars, tstep, Temp, obs){
  a   <- pars[1]
  Teq <- pars[2]
  SL0 <- pars[3]
  np      <- length(Temp)
  gmsl    <- rep(NA, np)
  gmsl[1] <- SL0
  
  # Use Forward Euler to estimate sea level over time.
  for (i in 2:np) {
    gmsl[i] <- gmsl[i-1]  +  tstep * a * ( Temp[i-1]-Teq )
  }
  resid <- obs - gmsl
  # Return the root mean square error
  rmse <- sqrt(mean(resid^2))
  return(rmse) 
}
# Plug in random values for the parameters.
parameter_guess <- c(3.4, -0.5, slr[1])
# Optimize the model to find initial starting values for parameters.
result <- optim(parameter_guess, fn, gr=NULL, tstep, Temp = hist.temp, obs = slr)
start_alpha <- result$par[1]
start_Teq <- result$par[2]
start_SL0 <- result$par[3]
# Make a plot of observations and fit.
slr.est <- gmsl_r07(a = start_alpha, Teq = start_Teq, SL0 = start_SL0, tstep, 
                    Tg = hist.temp)
plot(year, slr, pch = 20, ylab = "Sea-level Anomaly [mm]", xlab = "Year")
lines(year, slr.est, col = "blue", lwd = 2)
# Calculate residuals and initial sigma value.
res <- slr - slr.est
start_sigma <- sd(res)

# Apply the auto-correlation function to determine correlation coefficients.
rho <- rep(NA,3)
ac <- acf(res, lag.max = 5, plot = TRUE, main = "", xlab = "Time lag", 
          ylab = "Autocorrelation function")
rho[1] <- ac$acf[1]
rho[2] <- ac$acf[2]
rho[3] <- ac$acf[3]
rho[4] <- ac$acf[4]
rho[5] <- ac$acf[5]
start_rho <- rho[2]

### Setting up the prior information
# Set up MCMC calibration.
# Set up priors.
bound.lower <- c( 0, -3, err_neg[1], 0, -0.99)
bound.upper <- c(20,  2, err_pos[1], 10,  0.99)
# Name the parameters and specify the number of physical model parameters.
# Sigma and rho are statistical parameters and will not be counted in the number.
parnames <- c("alpha", "Teq", "SL0", "sigma", "rho")
model.p <- 3
# Set the measurement errors.
y.meas.err <- err.obs
# Load the likelihood model accounting for the error structure.
source("Rahm_obs_likelihood_AR.R")
# Setup initial parameters.
p <- c(start_alpha, start_Teq, start_SL0, start_sigma, start_rho) 
# Set p0 to a vector of values to be used in optimizing the log-posterior function.
# We optimize the negative log-posterior function because “optim” only knows how to 
# minimize a function, while we would like the log-posterior to be maximized.
p0 <- c(3.4, -0.5, slr[1], 6, 0.5)
p0 <- optim(p0, function(p) -log.post(p))$par
print(round(p0, 4))

### Running MCMC
# Set optimal acceptance rate as # parameters->infinity 
# (Gelman et al, 1996; Roberts et al, 1997)
# Set number of iterations and rate of adaptation 
# (gamma.mcmc, between 0.5 and 1, lower is faster adaptation)
accept.mcmc = 0.234
gamma.mcmc = 0.5
NI <- 4e5
step <- c(0.2, 0.02, 1, 0.1, 0.01)
# Run MCMC calibration. Note: Runs for several minutes.
mcmc.out <- MCMC(log.post, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc, 
                 gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains <- mcmc.out$samples

### Determine the burn-in period
## Gelman and Rubin's convergence diagnostic:
set.seed(1708)
p0 <- c(1.9, -0.9, -145, 4, 0.7) # Arbitrary choice.
mcmc.out2 <- MCMC(log.post, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc, 
                  gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains2 <- mcmc.out2$samples
set.seed(1234)
p0 <- c(2.9, 0, -160, 5, 0.8) # Arbitrary choice.
mcmc.out3 <- MCMC(log.post, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc, 
                  gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains3 <- mcmc.out3$samples
# Turn the chains into an mcmc object.
# This is a requirement for the gelman.diag and gelman.plot commmand.
mcmc1 <- as.mcmc(mcmc.chains) 
mcmc2 <- as.mcmc(mcmc.chains2) 
mcmc3 <- as.mcmc(mcmc.chains3)
# Test for convergence.
set.seed(111) # revert back to original seed
mcmc_chain_list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
gelman.diag(mcmc_chain_list)
gelman.plot(mcmc_chain_list)
# Create a vector to test the statistic at several spots throughout the chain.
# Test from 5000 to 4e5 in increments of 5000.
niter.test <- seq(from = 5000, to = NI, by = 5000)
gr.test <- rep(NA, length(niter.test))
# Once the statistic is close to 1 (adopting the standard of < 1.1), the between-chain 
# variability is indistinguishable from the within-chain variability, and hence converged.
for (i in 1:length(niter.test)){
  mcmc1 = as.mcmc(mcmc.chains[1:niter.test[i],])
  mcmc2 = as.mcmc(mcmc.chains2[1:niter.test[i],])
  mcmc3 = as.mcmc(mcmc.chains3[1:niter.test[i],])
  mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2, mcmc3))
  gr.test[i] = gelman.diag(mcmc_chain_list)[2]
}
# Plot Gelman and Rubin statistics as a function of iterations. The iterations prior  
# to the point in which convergence happens is what we set as the burn-in.
plot(niter.test, gr.test, pch=20)
abline(h = 1.1, col = "red")
gr.test < 1.1; first.stat = max(which((gr.test < 1.1) == FALSE))+1; print(first.stat)
# Set the burn-in to the 1st statistic that remains under 1.1 and remove it.
burnin <- seq(1, niter.test[first.stat], 1)
slr.burnin.chain <- mcmc.chains[-burnin, ]

### Hindcasting and projecting sea-level rise
# Extract parameter vectors from the chain to enhance code readability.
alpha.chain <- slr.burnin.chain[ ,1]
Teq.chain <-   slr.burnin.chain[ ,2]
SL_0.chain <-  slr.burnin.chain[ ,3]
sigma.chain <- slr.burnin.chain[ ,4]
rho.chain <-   slr.burnin.chain[ ,5]
# Set up empty matrices for sea-level output.
# Loop over the sea level model to generate a distribution of sea level simulations.
NI_length <- length(slr.burnin.chain[ ,1])
SLR.model.sim <- mat.or.vec(NI_length, length(project.years))
for(n in 1:NI_length) {
  SLR.model.sim[n, ] = gmsl_r07(a = alpha.chain[n], Teq = Teq.chain[n], 
                                SL0 = SL_0.chain[n], tstep, Tg = rcp85)
}
# Estimate the residuals with the AR(1) coefficient and sigma.
SLR.residuals <- mat.or.vec(NI_length, length(project.years)) #(nr,nc)
for(n in 1:NI_length) {
  for(i in 2:length(project.years)) {
    SLR.residuals[n,i] <- rho.chain[n]*SLR.residuals[n,i-1] + 
      rnorm(1, mean = 0, sd = sigma.chain[n])
  }
}
# Estimate the hindcasts and projections: add the residuals onto the model simulations.
SLR.projections <- mat.or.vec(NI_length, length(project.years)) #(nr,nc)
for(i in 1:NI_length) {
  SLR.projections[i,] <- SLR.model.sim[i,] + SLR.residuals[i,]
}


