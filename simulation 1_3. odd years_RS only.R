rm(list=ls())
library(boot)
library(nimble)
library(rstan)

#setwd('')

nsim <- 100
#==============
# Basic values
#==============
nsite_imbcr <- 100 # number of sites with IMBCR data
nsite       <- 100 # total number of sites
nyear       <- 20  # number of years

### Process
beta_lambda0 <- c(log(10), .2, -.8) # initial abundance, log mean,                     effect of temperature, effect of square temperature
beta_rho <- c(log(1), -.4, .2, -.6) # population growth, log mean, density dependence, effect of temperature, effect of square temperature
delta <- 1.5                        # immigration, expectation

### IMBCR
nreps <- 4    # number of points
ndist <- 5    # number of distance bins
ntime <- 6    # number of time bins
cutoff <- 200 # maximum distance for distance sampling
nobsv <- 30   # number of observers

alpha_sigma <- c(log(60)  , .6, -.4) # distance sampling, log   mean, effect of w1, effect of w2
alpha_theta <- c(logit(.4), .8, -.6) # removal  sampling, logit mean, effect of w1, effect of w2

#===============
# Simulate data
#===============
years_imbcr <- rep(c(1,0), times=nyear/2)

for (sim in 1:nsim) {

### Process covariates
x <- array(, dim=c(nsite, nyear, 2))
x[,,1] <- rnorm(nsite*nyear, 0, 1)
x[,,2] <- x[,,1] ^ 2

### Abundance
N <- matrix(, nsite, nyear)
lambda0 <- exp(cbind(1, x[,1,]) %*% beta_lambda0)
N[,1] <- rpois(nsite, lambda0)

rho <- matrix(, nsite, nyear-1)
for (t in 2:nyear) {
  rho[,t-1] <- exp(cbind(1, (N[,t-1] - lambda0) / lambda0, x[,t,]) %*% beta_rho)
  N[,t] <- rpois(nsite, N[,t-1] * rho[,t-1] + delta)
} # t

### IMBCR data
w <- array(rnorm(nsite*nyear*nreps*2, 0, 1), dim=c(nsite, nyear, nreps, 2))
sigma <- theta <- array(, dim=c(nsite, nyear, nreps))
for (t in 1:nyear) {
  for (j in 1:nreps) {
    sigma[,t,j] <- exp      (cbind(1, w[,t,j,]) %*% alpha_sigma)
    theta[,t,j] <- inv.logit(cbind(1, w[,t,j,]) %*% alpha_theta)
  } # j
} # t

breaks <- seq(0, cutoff, length.out=ndist + 1)
psi <- array(, dim=c(nsite, nyear, nreps, ndist))
for(d in 1:ndist) {
  psi[,,,d] <- 2 * (sigma^2 * (exp(-1*(breaks[d]^2)/(2*sigma^2)) - exp(-1*(breaks[d+1]^2)/(2*sigma^2)))) / (cutoff ^ 2)
} # k
psi_sum <- apply(psi, 1:3, sum)
psi_prop <- array(, dim=c(nsite, nyear, nreps, ndist))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (j in 1:nreps) {
      psi_prop[i,t,j,] <- psi[i,t,j,] / psi_sum[i,t,j]
    } # j
  } # t
} # i

phi <- array(, dim=c(nsite, nyear, nreps, ntime))
for (r in 1:ntime) {
  phi[,,,r] <- (1 - theta) ^ (r - 1) * theta
} # r
phi_sum <- apply(phi, 1:3, sum)
phi_prop <- array(, dim=c(nsite, nyear, nreps, ntime))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (j in 1:nreps) {
      phi_prop[i,t,j,] <- phi[i,t,j,] / phi_sum[i,t,j]
    } # j
  } # t
} # i

pi <- array(, dim=c(nsite, nyear, nreps, ndist * ntime))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (j in 1:nreps) {
      for (d in 1:ndist) {
        for (r in 1:ntime) {
          pi[i,t,j,(d-1)*ntime+r] <- psi_prop[i,t,j,d] * phi_prop[i,t,j,r]
        } # r
      } # d
    } # j
  } # t
} # i

imbcr_cnt_sum <- array(, dim=c(nsite, nyear, nreps))
imbcr_cnt <- array(, dim=c(nsite, nyear, nreps, ndist*ntime))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    if (years_imbcr[t] == 1) {
      for (j in 1:nreps) {
        imbcr_cnt_sum[i,t,j] <- rbinom(1, N[i,t], psi_sum[i,t,j] * phi_sum[i,t,j])
        imbcr_cnt[i,t,j,] <- rmultinom(1, imbcr_cnt_sum[i,t,j], pi[i,t,j,])
      } # j
    }
  } # t
} # i

#========================
# Define model in Nimble
#========================
code <- nimbleCode({
      
  # Priors
  for (k in 1:3) {
    beta_lambda0[k] ~ dnorm(0, sd=10)
  } # k
  for (k in 1:4) {
    beta_rho[k] ~ dnorm(0, sd=10)
  } # k
  delta ~ dgamma(.01, .01)
  for (k in 1:3) {
    alpha_sigma[k] ~ dnorm(0, sd=10)
  } # k
  for (k in 1:3) {
    alpha_theta[k] ~ dnorm(0, sd=10)
  } # k

  # Process model
  for (i in 1:nsite) {
    lambda0[i] <- exp(beta_lambda0[1] + beta_lambda0[2] * x[i,1,1] + beta_lambda0[3] * x[i,1,2])
    N[i,1] ~ dpois(lambda0[i])
    for (t in 2:nyear) {
      rho[i,t-1] <- exp(beta_rho[1] + beta_rho[2] * (N[i,t-1] - lambda0[i]) / lambda0[i] + 
                        beta_rho[3] * x[i,t,1] + beta_rho[4] * x[i,t,2])
      N[i,t] ~ dpois(N[i,t-1] * rho[i,t-1] + delta)
    } # t
  } # i

  # Observation model, IMBCR
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        sigma[i,t,j] <- exp   (alpha_sigma[1] + alpha_sigma[2] * w[i,t,j,1] + alpha_sigma[3] * w[i,t,j,2])
        theta[i,t,j] <- ilogit(alpha_theta[1] + alpha_theta[2] * w[i,t,j,1] + alpha_theta[3] * w[i,t,j,2])
      } # j
    } # t
  } # i

  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        for (d in 1:ndist) {
          intval[i,t,j,d] <- sigma[i,t,j]^2 * (exp(-1*(breaks[d]^2)/(2*sigma[i,t,j]^2)) - exp(-1*(breaks[d+1]^2)/(2*sigma[i,t,j]^2)))
          psi[i,t,j,d] <- 2 * intval[i,t,j,d] / (cutoff ^ 2)
        } # d
        psi_sum[i,t,j] <- sum(psi[i,t,j,1:ndist])
      } # j
    } # t
  } # i
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        for(d in 1:ndist) {
          psi_prop[i,t,j,d] <- psi[i,t,j,d] / psi_sum[i,t,j]
        } # d
      } # j
    } # t
  } # i

  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        for (r in 1:ntime) {
          phi[i,t,j,r] <- (1 - theta[i,t,j]) ^ (r - 1) * theta[i,t,j]
        } # r
        phi_sum[i,t,j] <- sum(phi[i,t,j,1:ntime])
      } # j
    } # t
  } # i
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        for (r in 1:ntime) {
          phi_prop[i,t,j,r] <- phi[i,t,j,r] / phi_sum[i,t,j]
        } # r
      } # j
    } # t
  } # i

  for (i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        for (d in 1:ndist) {
          for (r in 1:ntime) {
            pi[i,t,j,(d-1)*ntime+r] <- psi_prop[i,t,j,d] * phi_prop[i,t,j,r]
          } # r
        } # d
      } # j
    } # t
  } # i

  for(i in 1:nsite) {
    for (t in 1:nyear) {
      for (j in 1:nreps) {
        imbcr_cnt_sum[i,t,j] ~ dbinom(psi_sum[i,t,j] * phi_sum[i,t,j], N[i,t])
        imbcr_cnt[i,t,j,1:(ndist*ntime)] ~ dmultinom(pi[i,t,j,1:(ndist*ntime)], imbcr_cnt_sum[i,t,j])
      } # j
    } # t
  } # i

}) # nimbleCode

#============
# Run Nimble
#============
# Data
constants <- list(
  nsite=nsite, nyear=nyear, 
  nreps=nreps, ndist=ndist, ntime=ntime, breaks=breaks, cutoff=cutoff, 
  x=x, w=w
)

data <- list(
  imbcr_cnt_sum=imbcr_cnt_sum, imbcr_cnt=imbcr_cnt
)

# Initial values
Ni <- matrix(apply(apply(imbcr_cnt, 1:3, sum), 1:2, max)*2+20, nsite, nyear)
inits <- list(
  beta_lambda0=rep(0,3), 
  beta_rho=rep(0,4), delta=1, 
  alpha_sigma=c(4,0,0), alpha_theta=rep(0,3),  
  N=Ni)

model <- nimbleModel(code, constants=constants, data=data, inits=inits)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

mcmcConf$removeSamplers(c('beta_lambda0', 'beta_rho', 'alpha_sigma', 'alpha_theta'))
mcmcConf$printSamplers()

mcmcConf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmcConf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_sigma'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_theta'), type = "RW_block")
mcmcConf$printSamplers()

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

fit <- runMCMC(compiled$mcmc, nchains = 3, niter = 5000, nburnin = 0)

} # sim


