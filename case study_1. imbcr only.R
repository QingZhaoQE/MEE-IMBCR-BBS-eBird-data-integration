rm(list=ls())
library(nimble)
library(rstan)

#setwd('')

#==============
# Read in data
#==============
imbcr <- read.csv('data/ready/imbcr.csv')
gra <- read.csv('data/ready/grass.csv')
tem <- read.csv('data/ready/temperature.csv')

### Standardize covariates
x1 <- as.matrix(gra[,-c(1,2)])
x2 <- as.matrix(tem[,-c(1,2)])
x1 <- (x1 - 0.5) * 4
x2 <- (x2 - mean(x2)) / sd(x2)
gra[,-c(1,2)] <- x1
tem[,-c(1,2)] <- x2

gra <- gra[which(gra$fish_id %in% imbcr$fish_id),]
tem <- tem[which(tem$fish_id %in% imbcr$fish_id),]

info <- data.frame(fish_id = gra$fish_id, site = 1:dim(gra)[1])
gra <- as.matrix(gra[,-c(1:10)])
tem <- as.matrix(tem[,-c(1:10)])
nsite <- dim(gra)[1]
nyear <- dim(gra)[2]

imbcr_t <- merge(info, imbcr)
imbcr_t <- imbcr_t[order(imbcr_t$year, imbcr_t$site),]
nobs_imbcr <- dim(imbcr_t)[1]
years_imbcr <- imbcr_t$year
sites_imbcr <- imbcr_t$site
imbcr_cnt <- as.matrix(imbcr_t[,-c(1:4)])
imbcr_cnt_sum <- rowSums(imbcr_cnt)
ndist <- 20    # number of distance bins
ntime <- 6    # number of time bins
cutoff <- 200 # maximum distance for distance sampling
breaks <- seq(0, cutoff, length.out=ndist + 1)

year0 <- min(years_imbcr)
years_imbcr <- years_imbcr - year0 + 1

#========================
# Define model in Nimble
#========================
code <- nimbleCode({

  # Priors
  zeta ~ dunif(0, 1)
  for (k in 1:4) {
    beta_lambda0[k] ~ dnorm(0, sd=10)
    beta_rho[k] ~ dnorm(0, sd=10)
  } # k
  sigma ~ dgamma(.01, .01)
  theta ~ dunif(0, 1)

  # Process model
  for (i in 1:nsite) {
    z[i] ~ dbinom(zeta, 1)
    lambda0[i] <- exp(
      beta_lambda0[1] + 
      beta_lambda0[2] * gra[i,1] + 
      beta_lambda0[3] * tem[i,1] + 
      beta_lambda0[4] * tem[i,1] * tem[i,1])
    N[i,1] ~ dpois(lambda0[i] * z[i] + 5e-3 * (1 - z[i]))
    for (t in 2:nyear) {
      rho[i,t-1] <- exp(
        beta_rho[1] + 
        beta_rho[2] * gra[i,t] + 
        beta_rho[3] * tem[i,t] + 
        beta_rho[4] * tem[i,t] * tem[i,t])
      N[i,t] ~ dpois(N[i,t-1] * rho[i,t-1] * z[i] + 5e-3 * (1 - z[i]))
    } # t
  } # i

  # Observation model, IMBCR
  for (d in 1:ndist) {
    intval[d] <- sigma^2 * (exp(-1*(breaks[d]^2)/(2*sigma^2)) - exp(-1*(breaks[d+1]^2)/(2*sigma^2)))
    psi[d] <- 2 * intval[d] / (cutoff ^ 2)
  } # d
  psi_sum <- sum(psi[1:ndist])
  for(d in 1:ndist) {
    psi_prop[d] <- psi[d] / psi_sum
  } # d

  for (r in 1:ntime) {
    phi[r] <- (1 - theta) ^ (r - 1) * theta
  } # r
  phi_sum <- sum(phi[1:ntime])
  for (r in 1:ntime) {
    phi_prop[r] <- phi[r] / phi_sum
  } # r

  for (d in 1:ndist) {
    for (r in 1:ntime) {
      pi[(r-1)*ndist+d] <- psi_prop[d] * phi_prop[r]
    } # r
  } # d

  for(k in 1:nobs_imbcr) {
    imbcr_cnt_sum[k] ~ dbinom(psi_sum * phi_sum, N[sites_imbcr[k], years_imbcr[k]])
    imbcr_cnt[k,1:(ndist*ntime)] ~ dmultinom(pi[1:(ndist*ntime)], imbcr_cnt_sum[k])
  } # k

}) # nimbleCode

#============
# Run Nimble
#============
# Data
constants <- list(
  nsite=nsite, nyear=nyear, 
  ndist=ndist, ntime=ntime, breaks=breaks, cutoff=cutoff, 
  nobs_imbcr=nobs_imbcr, sites_imbcr=sites_imbcr, years_imbcr=years_imbcr, 
  gra=gra, tem=tem
)

data <- list(
  imbcr_cnt_sum=imbcr_cnt_sum, imbcr_cnt=imbcr_cnt
)

# Initial values
Ni <- matrix(5, nsite, nyear)
for (i in 1:nsite) {
  for (t in 1:nyear) {
    if (length(imbcr_cnt_sum[which(sites_imbcr == i & years_imbcr == t)]) > 0) {
      Ni[i,t] <- max(imbcr_cnt_sum[which(sites_imbcr == i & years_imbcr == t)])
    }
  } # t
} # i
zi <- ifelse(rowSums(Ni) == 0, 0, 1)
Ni <- Ni * 2 + 10
inits <- list(
  zeta=0.2, 
  beta_lambda0=c(2,0.5,-1,-1), beta_rho=c(0.1,0,-0.1,-0.1), 
  sigma=60, theta=0.4, 
  N=Ni, z=zi)

model <- nimbleModel(code, constants=constants, data=data, inits=inits)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

mcmcConf$removeSamplers(c('beta_lambda0', 'beta_rho'))
mcmcConf$printSamplers()

mcmcConf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmcConf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmcConf$printSamplers()

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

chain <- 3
nmcmc <- 100000
fit <- runMCMC(compiled$mcmc, nchains = chain, niter = nmcmc, nburnin = 0)

dev.off()


