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
#whether a patient had a DLT or not 
s <- c(0,0,0,0)
#current dose for each patient
d <- c(80, 80, 160, 160)
#reference dose
D.ref <- 960
#binary indicator I (0=montherapy, 1 = combination therapy)
#I.combo <- c(0,0,0,0)
#I.combo <- c(rep(0, length(d)), rep(1, length(d)))

#vector of prior means 
mu0 <- -1.75
mu1 <- 0.00
mu2 <- 0.00
priorMean <- c(mu0, mu1, mu2)

#covariance matrix
sigma0 <- 1.60
sigma1 <- 0.2
sigma2 <- 1.5
priorVar <- matrix(c(sigma0^2, 0, 0,
                     0, sigma1^2, 0,
                     0, 0, sigma2^2), ncol = 3, byrow = TRUE)

priorPrec <- solve(priorVar)

#iterations
iter<-4*10^6
#overdose probability threshold
c.overdose<-0.40 
#number of DLTs on the current dose levels per each patient - no.DLT and s are essentially the same, no? 
no.DLT = c(0,0,0,0,0)

############################################################################################################

blrm.data <- list(
  D.ref = D.ref,
  s = s,
  d = d,
  N = N,
  priorPrec = priorPrec,
  priorMean = priorMean
)



# Model specification for JAGS
blrm.model.string <- "
model {
  ## Priors ##
  theta[1:3] ~ dmnorm(priorMean[], priorPrec[1:3,1:3])

  
  alpha0 <- theta[1]
  alpha1 <- exp(theta[2])
  alpha2 <- exp(theta[3])

  ## Likelihood ##
  for (j in 1:N) {
    lin.pred.p.1[j] <- alpha0 + alpha1 * log(d[j]/D.ref)
    lin.pred.p.2[j] <- ifelse(lin.pred.p.1[j] < -10, -10, ifelse(lin.pred.p.1[j] > 10, 10, lin.pred.p.1[j]))
    p[j] <- exp(lin.pred.p.2[j]) / (1 + exp(lin.pred.p.2[j])) 
    s[j] ~ dbern(p[j])
    
  }
}
"

blrm.model <- jags.model(textConnection(blrm.model.string), data = blrm.data, n.chains = 1)
update(blrm.model, iter)
samples <- jags.samples(blrm.model, variable.names = c("alpha0", "alpha1", "alpha2"), n.iter = iter)

# Extract posterior samples
a01 <- samples$alpha0[1, ,]
a11 <- samples$alpha1[1, , ]
a21 <- samples$alpha2[1, , ]



########################################
#create empty vectors to store estimates 
LCI<-UCI<-Pi.MTD<-Pi.MTD.store<-Pi.above<-toxicity<-mat.or.vec(length(doses),1)
all.Pi <-mat.or.vec(iter,length(doses))
allstop <-0

for (j in 1:length(doses)) {
  LP <- a01 + a11 * log(doses[j]/D.ref)
  #make sure LP stays between -10 and 10
  LP <- pmin(pmax(LP, -10), 10)
  #calculate the odds
  odds <- exp(LP)
  #calculate the predicted probabilities for specific j
  Pi <- odds / (1 + odds)
  
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

#create empty vectors to store new variables
Pi.MTD.mono <- Pi.MTD[1:length(doses)]
Pi.MTD.mono.store <- Pi.MTD.store[1:length(doses)]
Pi.above.mono <- Pi.above[1:length(doses)]
toxicity.mono <- toxicity[1:length(doses)]
n.mono <- n.patients.per.dose[1:length(doses)]
s.mono <- s[1:length(doses)]
LCI.mono <- LCI[1:length(doses)]
UCI.mono <- UCI[1:length(doses)]

###