rm(list=ls())
library(nimble)
library(rstan)

#setwd('')

#==============
# Read in data
#==============
imbcr <- read.csv('data/ready/imbcr.csv')
bbs <- read.csv('data/ready/bbs.csv')
gra <- read.csv('data/ready/grass.csv')
tem <- read.csv('data/ready/temperature.csv')

### Standardize covariates
x1 <- as.matrix(gra[,-c(1,2)])
x2 <- as.matrix(tem[,-c(1,2)])
x1 <- (x1 - 0.5) * 4
x2 <- (x2 - mean(x2)) / sd(x2)
gra[,-c(1,2)] <- x1
tem[,-c(1,2)] <- x2

gra <- gra[which(gra$fish_id %in% sort(unique(c(imbcr$fish_id, bbs$fish_id)))),]
tem <- tem[which(tem$fish_id %in% sort(unique(c(imbcr$fish_id, bbs$fish_id)))),]

info <- data.frame(fish_id = gra$fish_id, site = 1:dim(gra)[1])
gra <- as.matrix(gra[,-c(1:2)])
tem <- as.matrix(tem[,-c(1:2)])
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

bbs_t <- merge(info, bbs)
bbs_t <- bbs_t[order(bbs_t$year, bbs_t$site),]
nobs_bbs <- dim(bbs_t)[1]
years_bbs <- bbs_t$year
sites_bbs <- bbs_t$site
route_bbs <- as.numeric(as.factor(bbs_t$route))
nroute_bbs <- length(unique(route_bbs))
obser_bbs <- as.numeric(as.factor(bbs_t$obs))
nobser_bbs <- length(unique(obser_bbs))
new_bbs <- bbs_t$new
bbs_cnt_ori <- as.matrix(bbs_t[,-c(1:6)])
bbs_cnt <- matrix(0, nobs_bbs, 10)
for (i in 1:10) {
  bbs_cnt[,i] <- rowSums(bbs_cnt_ori[,(i-1)*5+1:5])
} # i
nstop_bbs <- dim(bbs_cnt)[2]
stops_bbs <- seq(-2, 2, length.out=nstop_bbs)

year0 <- min(c(years_imbcr, years_bbs))
years_imbcr <- years_imbcr - year0 + 1
years_bbs <- years_bbs - year0 + 1

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
  chi_mu ~ dnorm(0, sd=10)
  chi_route_sd ~ dgamma(.01, .01)
  chi_obser_sd ~ dgamma(.01, .01)
  chi_new ~ dnorm(0, sd=10)
  chi_stop_beta[1] ~ dnorm(0, sd=10)
  chi_stop_beta[2] ~ T(dnorm(0, sd=10), 0,)

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

  # Observation model, BBS
  for (k in 1:nroute_bbs) {
    chi_route_epsilon[k] ~ dnorm(0, sd=chi_route_sd)
  } # k
  for (k in 1:nobser_bbs) {
    chi_obser_epsilon[k] ~ dnorm(0, sd=chi_obser_sd)
  } # k
  for (j in 1:nstop_bbs) {
    chi_stop_sd[j] <- exp(chi_stop_beta[1] + chi_stop_beta[2] * stops_bbs[j])
  } # j
  for (k in 1:nobs_bbs) {
    log_chi_mu[k] <- 
      chi_mu + 
      chi_route_epsilon[route_bbs[k]] + 
      chi_obser_epsilon[obser_bbs[k]] + 
      chi_new * new_bbs[k]
    for (j in 1:nstop_bbs) {
      log_chi[k,j] ~ dnorm(log_chi_mu[k], chi_stop_sd[j])
      chi[k,j] <- exp(log_chi[k,j])
      bbs_cnt[k,j] ~ dpois(N[sites_bbs[k], years_bbs[k]] * chi[k,j])
    } # j
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
  nobs_bbs=nobs_bbs, sites_bbs=sites_bbs, years_bbs=years_bbs, 
  nroute_bbs=nroute_bbs, route_bbs=route_bbs, nobser_bbs=nobser_bbs, obser_bbs=obser_bbs, 
  new_bbs=new_bbs, nstop_bbs=nstop_bbs, stops_bbs=stops_bbs, 
  gra=gra, tem=tem
)

data <- list(
  imbcr_cnt_sum=imbcr_cnt_sum, imbcr_cnt=imbcr_cnt, 
  bbs_cnt=bbs_cnt
)

# Initial values
Ni <- matrix(5, nsite, nyear)
for (i in 1:nsite) {
  for (t in 1:nyear) {
    if (length(c(imbcr_cnt_sum[which(sites_imbcr == i & years_imbcr == t)], bbs_cnt[which(sites_bbs == i & years_bbs == t),])) > 0) {
      Ni[i,t] <- max(c(imbcr_cnt_sum[which(sites_imbcr == i & years_imbcr == t)], bbs_cnt[which(sites_bbs == i & years_bbs == t),]))
    }
  } # t
} # i
zi <- ifelse(rowSums(Ni) == 0, 0, 1)
Ni <- Ni * 2 + 10
inits <- list(
  zeta=0.2, 
  beta_lambda0=c(2,0.5,-1,-1), beta_rho=c(0.1,0,-0.1,-0.1), 
  sigma=60, theta=0.4, 
  chi_mu=-6, chi_route_sd=1, chi_obser_sd=1, chi_new=0.5, chi_stop_beta=c(-1.5,0.2), 
  N=Ni, z=zi)

model <- nimbleModel(code, constants=constants, data=data, inits=inits)
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

mcmcConf$removeSamplers(c('beta_lambda0', 'beta_rho', 'chi_stop_beta'))
mcmcConf$printSamplers()

mcmcConf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmcConf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmcConf$addSampler(target = c('chi_stop_beta'), type = "RW_block")
mcmcConf$printSamplers()

mcmc <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, mcmc)

chain <- 3
nmcmc <- 100000
fit <- runMCMC(compiled$mcmc, nchains = chain, niter = nmcmc, nburnin = 0)

dev.off()


