rm(list=ls())
library(truncnorm)
library(boot)
library(nimble)
library(rstan)

#setwd('')

x_quantile <- .9 # quantile of standard normal to set lower bound for x, 0.5, 0.75 or 0.9

nsim <- 100
#==============
# Basic values
#==============
nsite_imbcr <- 100 # number of sites with IMBCR data
nsite_bbs   <- 150 # number of sites with BBS data
nsite_ebird <- 200 # number of sites with ebird data
nsite       <- 300 # total number of sites
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

### BBS
bbs_chi_mean <- 1.5  # mean scaling factor for BBS observation # 1.5 or 0.5
bbs_chi_sd <- .3     # SD of scaling factors for BBS observation
bbs_chi_eta <- 1.2   # new observer effect of scaling factors for BBS observation
alpha_bbs_chi <- c(.2, -.4)  # effect of z for scaling factors for BBS observation

### eBird
alpha_ebird_chi <- c(log(5), -.5, .4, .6) # scaling factor for eBird observation, log mean, effects of stationary survey, survey distance and duration

#===============
# Simulate data
#===============
### Sites
sites_imbcr <- rep(c(0,0,1), times=nsite/3)
sites_bbs <- rep(c(0,1), times=nsite/2)
sites_ebird <- rep(0, nsite)
sites_ebird[which(sites_imbcr == 0 & sites_bbs == 0)] <- 1
sites_ebird[which(sites_imbcr == 0 & sites_bbs == 1)] <- rep(c(0,1), times=length(which(sites_imbcr == 0 & sites_bbs == 1))/2)
sites_ebird[which(sites_imbcr == 1 & sites_bbs == 0)] <- rep(c(0,1), times=length(which(sites_imbcr == 1 & sites_bbs == 0))/2)
sites_ebird[which(sites_imbcr == 1 & sites_bbs == 1)] <- rep(c(0,1), times=length(which(sites_imbcr == 1 & sites_bbs == 1))/2)
sites <- data.frame(sites_imbcr, sites_bbs, sites_ebird)
#table(data.frame(sites_imbcr, sites_bbs, sites_ebird))

years_imbcr <- c(rep(0,nyear/2), rep(1,nyear/2))
years_ebird <- c(rep(0,nyear/4), rep(1,nyear*3/4))

for (sim in 1:nsim) {

### Process covariates
imbcr_lower <- qnorm(x_quantile) # lower bound of x for imbcr data
x <- array(, dim=c(nsite, nyear, 2))
x[which(sites_imbcr==1),,1] <- rtruncnorm(n=nsite_imbcr*nyear, a=imbcr_lower, b=Inf, mean=0, sd=1)
x[which(sites_imbcr==0),,1] <- rnorm((nsite-nsite_imbcr)*nyear, 0, 1)
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
  if (sites_imbcr[i] == 1) {
    for (t in 1:nyear) {
      if (years_imbcr[t] == 1) {
        for (j in 1:nreps) {
          imbcr_cnt_sum[i,t,j] <- rbinom(1, N[i,t], psi_sum[i,t,j] * phi_sum[i,t,j])
          imbcr_cnt[i,t,j,] <- rmultinom(1, imbcr_cnt_sum[i,t,j], pi[i,t,j,])
        } # j
      }
    } # t
  }
} # i

### BBS data
obs <- numeric(nsite)
obs[which(sites_bbs==1)] <- sample(1:(nobsv/2), nsite_bbs, replace=T)
obs <- matrix(obs, nsite, nyear, byrow=F)
for (i in 1:nsite) {
  if (sites_bbs[i] == 1) {
    tt <- sample(5:15, 1)
    obs[i, tt:nyear] <- sample((nobsv/2+1):nobsv, 1)
  }
} # i

new_obs <- sort(sample(1:nobsv, nobsv/2, replace=F))
new <- matrix(0, nsite, nyear)
for (i in 1:nsite) {
  if (obs[i,1] %in% new_obs) {
    new[i,1] <- 1
  }
  for (t in 2:nyear) {
    if (obs[i,t] != obs[i,t-1]) {
      new[i,t] <- 1
    }
  } # t
} # i

obs[which(obs == 0)] <- nobsv + 1

z <- array(rnorm(nsite*nyear*2, 0, 1), dim=c(nsite, nyear, 2))

bbs_chi <- matrix(, nsite, nyear)
bbs_chi_obs <- exp(rnorm(nobsv+1, log(bbs_chi_mean), bbs_chi_sd))
bbs_chi_eta_mat <- matrix(, nsite, nyear)
for (i in 1:nsite) {
  for (t in 1:nyear) {
    bbs_chi_eta_mat[i,t] <- ifelse(new[i,t]==1, bbs_chi_eta, 1)
    bbs_chi[i,t] <- exp(log(bbs_chi_obs[obs[i,t]] * bbs_chi_eta_mat[i,t]) + z[i,t,] %*% alpha_bbs_chi)
  } # t
} # i

bbs_cnt <- matrix(, nsite, nyear)
for (i in 1:nsite) {
  if (sites_bbs[i] == 1) {
    for (t in 1:nyear) {
      bbs_cnt[i,t] <- rpois(1, N[i,t] * bbs_chi[i,t])
    } # t
  }
} # i

### eBird data
ebird_stat <- matrix(sample(c(rep(1,nsite*nyear*.2),rep(0,nsite*nyear*.8)), nsite*nyear, replace=F), nsite, nyear)
ebird_dist <- exp(matrix(rnorm(nsite*nyear,0,1), nsite, nyear))
ebird_dist[which(ebird_stat == 1)] <- 0
ebird_dist <- log(ebird_dist + 1)
ebird_dist <- (ebird_dist - mean(ebird_dist)) / sd(ebird_dist)
ebird_time <- matrix(rnorm(nsite*nyear,0,1), nsite, nyear)
v <- array(, dim=c(nsite, nyear, 3))
v[,,1] <- ebird_stat
v[,,2] <- ebird_dist
v[,,3] <- ebird_time

ebird_chi <- matrix(, nsite, nyear)
for (t in 1:nyear) {
  ebird_chi[,t] <- exp(cbind(1, v[,t,]) %*% alpha_ebird_chi)
} # t

ebird_cnt <- matrix(, nsite, nyear)
for (i in 1:nsite) {
  if (sites_ebird[i] == 1) {
    for (t in 1:nyear) {
      if (years_ebird[t] == 1) {
        ebird_cnt[i,t] <- rpois(1, N[i,t] * ebird_chi[i,t])
      }
    } # t
  }
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
  bbs_chi_mean ~ dgamma(.01, .01)
  bbs_chi_sd ~ dgamma(.01, .01)
  bbs_chi_eta ~ dgamma(.01, .01)
  for (k in 1:2) {
    alpha_bbs_chi[k] ~ dnorm(0, sd=10)
  } # k
  for (k in 1:4) {
    alpha_ebird_chi[k] ~ dnorm(0, sd=10)
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

  # Observation model, BBS
  for (k in 1:(nobsv+1)) {
    log_bbs_chi_obs[k] ~ dnorm(log(bbs_chi_mean), sd=bbs_chi_sd)
    bbs_chi_obs[k] <- exp(log_bbs_chi_obs[k])
  } # k

  bbs_chi_eta_vec[1] <- 1
  bbs_chi_eta_vec[2] <- bbs_chi_eta
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      bbs_chi[i,t] <- exp(log(bbs_chi_obs[obs[i,t]] * bbs_chi_eta_vec[new[i,t]+1]) + 
                          alpha_bbs_chi[1] * z[i,t,1] + alpha_bbs_chi[2] * z[i,t,2])
      bbs_cnt[i,t] ~ dpois(N[i,t] * bbs_chi[i,t])
    } # t
  } # i


  # Observation model, eBird
  for (i in 1:nsite) {
    for (t in 1:nyear) {
      ebird_chi[i,t] <- exp(alpha_ebird_chi[1] + alpha_ebird_chi[2] * v[i,t,1] + 
                            alpha_ebird_chi[3] * v[i,t,2] + alpha_ebird_chi[4] * v[i,t,3])
      ebird_cnt[i,t] ~ dpois(N[i,t] * ebird_chi[i,t])
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
  nobsv=nobsv, obs=obs, new=new, 
  x=x, w=w, z=z, v=v
)

data <- list(
  imbcr_cnt_sum=imbcr_cnt_sum, imbcr_cnt=imbcr_cnt, bbs_cnt=bbs_cnt, ebird_cnt=ebird_cnt
)

# Initial values
Ni <- matrix(apply(apply(imbcr_cnt, 1:3, sum), 1:2, max)*2+20, nsite, nyear)
Ni[which(is.na(Ni))] <- max(Ni, na.rm=T)
inits <- list(
  beta_lambda0=rep(0,3), 
  beta_rho=rep(0,4), delta=1, 
  alpha_sigma=c(4,0,0), alpha_theta=rep(0,3),  
  bbs_chi_mean=1, bbs_chi_sd=1, log_bbs_chi_obs=rep(0,nobsv+1), 
  bbs_chi_eta=1, alpha_bbs_chi=rep(0,2), 
  alpha_ebird_chi=rep(0,4), 
  N=Ni)

model <- nimbleModel(code, constants=constants, data=data, inits=inits)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

mcmcConf$removeSamplers(c('beta_lambda0', 'beta_rho', 'alpha_sigma', 'alpha_theta', 'alpha_bbs_chi', 'alpha_ebird_chi'))
mcmcConf$printSamplers()

mcmcConf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmcConf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_sigma'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_theta'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_bbs_chi'), type = "RW_block")
mcmcConf$addSampler(target = c('alpha_ebird_chi'), type = "RW_block")
mcmcConf$printSamplers()

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

fit <- runMCMC(compiled$mcmc, nchains = 3, niter = 5000, nburnin = 0)

} # sim


