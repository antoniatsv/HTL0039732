setwd("C:/Users/Antonia/OneDrive - University of Cambridge, MRC Biostatistics Unit/Desktop/HTL code/my code")

library(tidyverse)
library(rjags)
library(runjags)
library(purrr)
library(dplyr)


#dose levels 
doses <- c(20,40,80,160,320,640,960)
#number of patients per current dose
n.patients.per.dose <- c(0,0,2,2,0,0,0)
#total number of patients across all dose levels
N <- sum(n.patients.per.dose) 
#number of DLTs for each patient/dose level
s <- c(0,0,0,0,0)
#current dose for each patient
d <- c(80, 80, 160, 160)
#reference dose
D.ref <- 960
#binary indicator I (0=montherapy, 1 = combination therapy)
I.combo <- c(0,0,0,0)

#vector of prior means 
mu0 <- -1.75
mu1 <- 0.00
mu2 <- 0.00
prior.mean <- c(mu0, mu1, mu2)

#covariance matrix
sigma0 <- 1.60^2
sigma1 <- 0.2^2
sigma2 <- 1.5^2
corr <- 0
cov.matrix <- matrix(c(sigma0, sigma0*sigma1*corr, sigma0*sigma2*corr,
                       sigma0*sigma1*corr, sigma1, sigma1*sigma2*corr,
                       sigma0*sigma2*corr, sigma1*sigma2*corr, sigma2), ncol = 3, byrow = TRUE)
priorVar<-cov.matrix
#iterations
iter<-4*10^6
#overdose probability threshold
c.overdose<-0.40 
#number of DLTs on the current dose levels
no.DLT = 0

############################################################################################################

blrm.data <- list(
  D.ref = D.ref,
  d = d,
  I.combo = I.combo,
  N = N,
  mu0 = mu0,
  mu1 = mu1,
  mu2 = mu2,
  sigma0 = sigma0,
  sigma1 = sigma1,
  sigma2 = sigma2,
  priorVar = priorVar
)



############################################################################################################
blrm.model.string <- "
model {
  ## Priors ##
  alpha0 ~ dnorm(mu0, 1/sigma0^2)
  alpha1 ~ dnorm(mu1, 1/sigma1^2)
  log.alpha2 ~ dnorm(log(mu2), 1/sigma2^2)
  alpha2 <- exp(log.alpha2)

  Omega[1:3, 1:3] <- inverse(priorVar)

  ## Likelihood ##
  for (j in 1:N) {
    lin.pred.p.1[j] <- alpha0 + alpha1 * log(d[j]/D.ref) + alpha2 * I.combo[j]
    lin.pred.p.2[j] <- ifelse(lin.pred.p.1[j] < -10, -10, ifelse(lin.pred.p.1[j] > 10, 10, lin.pred.p.1[j])) # lin pred after applying constraints to prevent extreme values to make sure the values for p1 are within a certain range
    s[j] ~ dbern(p[j])
    p[j] <- exp(lin.pred.p.2[j]) / (1 + exp(lin.pred.p.2[j])) 
  }
}
"

#create a text connection object from the string
blrm.model.spec <- textConnection(blrm.model.string)
#data to be used by the blrm model 
blrm.data <- list(
  D.ref = D.ref,
  d = d,
  I.combo = I.combo,
  N = N,
  mu0 = mu0,
  mu1 = mu1,
  mu2 = mu2,
  sigma0 = sigma0,
  sigma1 = sigma1,
  sigma2 = sigma2,
  priorVar = priorVar
)
#create a JAGs model object 
jags <- jags.model(file = textConnection(blrm.model.string), data = blrm.data)
#MCMC iterations 'burn-in' to allow the sampler to reach its stationary distribution before collecting samples for inference
update(jags, iter)
#sampling from the posterior distribution of the specified alpha parameters from the JAGS model
tt <-jags.samples(jags,c('alpha0','alpha1','alpha2'),iter,progress.bar="none")
#extracting samples of the alpha parameters from the 'tt' object. it selects the first chain and all iterations for that chain. 
a01<-tt$alpha0[1,,]
a11<-(tt$alpha1[1,,])
a21<-(tt$alpha2[1,,])


########################################
#create empty vectors to store estimates 
LCI<-UCI<-Pi.MTD<-Pi.MTD.store<-Pi.above<-toxicity<-mat.or.vec(length(d),1)
all.Pi <-mat.or.vec(iter,length(d))
stop <-0

for (j in 1:length(d)) {
  LP <- a01 + a11 * log(d[j]/D.ref) + a21 * I.combo[j] 
  #make sure LP stays between -10 and 10
  LP <- pmin(pmax(LP, -10), 10)
  #calculate the odds
  odds <- exp(LP)
  #calculate the predicted probabilities for specific j
  Pi <- odds / (1 + odds)
  #Pi for all j's
  all.Pi[, j] <- Pi
  #calculate Pi for DLT for individual patient j at their current dose level
  Pi.MTD.store[j] <- Pi.MTD[j] <- mean(Pi < 0.30) - mean(Pi < 0.20)
  #probability of toxicity greater than 30%
  Pi.above[j] <- mean(Pi > 0.30)
  #average predicted probability of toxicity across all patients j
  toxicity[j] <- mean(Pi)
  #Lower CI
  LCI[j] <- quantile(Pi, 0.025)
  #Upper CI
  UCI[j] <- quantile(Pi, 0.075)
}

Pi.MTD[Pi.MTD == 0] <- 0.001
Pi.MTD[Pi.above > c.overdose] <- 0

#68

Pi.MTD.mono <- Pi.MTD[1:length(d)]
Pi.MTD.mono.store <- Pi.MTD.store[1:length(d)]
Pi.above.mono <- Pi.above[1:length(d)]
toxicity.mono <- toxicity[1:length(d)]
n.mono <- n.patients.per.dose[1:length(d)]
s.mono <- s[1:length(d)]
LCI.mono <- LCI[1:length(d)]
UCI.mono <- UCI[1:length(d)]

#######
