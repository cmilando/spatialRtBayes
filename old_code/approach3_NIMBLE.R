library(nimble)
library(tidyr)
library(ggplot2)
library(tidyverse)

incidence_data <- read.csv("https://raw.githubusercontent.com/zwzhou-biostat/hrt/main/data/example_data.csv")

# incidence_data <- read.csv("../raw_data/example_data.csv")
# 
# usethis::use_data(incidence_data, overwrite = T)


#incidence_data <- read.csv("example_data.csv")

## the serial interval
w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w / sum(w)
w

## related to how W is used in matrix multiplication
## w = 0 if tau is too large, see the simulation formula
w <- rev(c(w, rep(0, nrow(incidence_data) - length(w))))

w

## The mobility matrix
P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3, byrow = T)

## 
N = length(incidence_data$a)  
K = dim(P)[1]

dim(incidence_data) # t0 and then 214 by 3 regions

L <- 8 ## ? smoothing window?

N ## ? number of days?
K ## ? number of regions?
L ## ? 8 is the smoothing window ?

## expanded ?
incidence_data_expand <- array(NA, dim = c(N-L, K, L))
head(incidence_data_expand)

for(i in 1:(N-L)){
  for(k in 1:K){
    for(j in 1:L){
      incidence_data_expand[i, k, j] <- incidence_data[i + j - 1, k]
    }
  }
}

incidence_data_expand

## if you use 2 dimensions it wont smooth out in NIMBLE
## because NIMBLE doesn't recognize 

###########################################################
# a bunch of implied variables:
#
# N, L, K
# convolution: with N-L-1  and K 
# lambda: with N-L-1, L, K
# R, P, w
# x
# 
# functions
# dpois

## and this runs for NITER times

## convolution

rtCode <- nimbleCode({ 
  
  # Define relationships between nodes
  ## for 2 : Ndays - smoothing window
  for (i in 2:(N-L)) {
    ## across K regions
    
    ## (a) first you do all of convolution
    ## multiplies the serial interval with the incidence from the today and past
    for(k in 1:K){
      convolution[i - 1, k] <- sum(
        (lambda[1:(i - 1), k, 1] * (w[ (N - i + 2):N ]) ) ## << contains lambda, expected incidence each day
      )
    }
    ### (b) then you do all of lambda
    ## inside the sum
    for(k in 1:K){
      for(j in 1:L){
        lambda[i, k, j] <- sum( t(P[1:K, k]) * 
                                  R[i, 1:K] *  ## << this is what you are estimating, R(t)
                                  convolution[i-1, 1:K]
                                )
      }
    }
  }
  
  ## ^ this specifies how these relationship looks like
  
  ## the code below specifies how where you assume the distribution of x comes from
  ## this is the first thing that happens, and then it trickles down into the top part
  
  ## now this is strange, you are passing incidence data in
  ## but are not using it ...
  ## also x is not used in any of these other loops ...
  ## syntax for NIMBLE so this is an INPUT
  ## x follows this distribution with a parameter you want to estimate
  ## so when you specify, this is where Beta follows this distribution
  
  
  ## right so this is m_j(t) ~ Poisson(N_j(t))
  for(i in 1:(N-L)){
    for(k in 1:K){
      for(j in 1:L){
        x[i, k, j] ~ dpois(lambda[i, k, j])
      }
    }
  }
  
  # Set priors
  ## for each region
  for(k in 1:K) {
    ## for each of first L smoothing days
    for(j in 1:L) {
      lambda[1, k, j] <- 10
    }
    beta[k] ~ dnorm(1, 0.5) ## << sigma, because it has just K dimensions
    ## then for all days
    for(i in 1:N){
      alpha[i, k] ~ dnorm(0, 1) # beta, because it has 2 dimensions
      theta[i, k] ~ dnorm(alpha[i, k], beta[k]) # mean = alpha, sigma = beta
      
      ## and this is log R_j(t) ~ N(Beta_j(t), sigma_j)
      R[i, k] <- exp(theta[i, k])
    }
  }
  
})

# takes 
# 30 minutes

## the "~"

## just do GIBBS
## not sure what rejection criteria would make sense, so just try GIBBS

###########################################################
rtConsts <- list(N = length(incidence_data$a),
                 w = w,
                 K = dim(P)[1])

rtData <- list(x = incidence_data_expand,
               P = P)

rtInits <- list(alpha = matrix(1, rtConsts$N, rtConsts$K), 
                beta = rep(1, rtConsts$K),
                theta = matrix(0.1, rtConsts$N, rtConsts$K))

###########################################################
## Really no idea how long this is supposed to take
## and doesn't seem like "convolutions" is defined anywhere"
## or lambda
nimbleMCMC_samples <- nimbleMCMC(code = rtCode, 
                                 constants = rtConsts, 
                                 progressBar = T,
                                 data = rtData, 
                                 inits = rtInits,
                                 setSeed = 1,
                                 nburnin = 1000, niter = 2000,
                                 monitors = c("R", "lambda"),
                                 summary = T, samples = T)
###########################################################

### so where are the random draws occurring
## **  weekly?! R(t)
## ** 

#### RIGHT RIGHT RIGHT
## so the two outputs here are R[i, j]
## so R in space and time

## and lambda -- which I think is the expected incidence (cases)

dat <- data.frame(a = nimbleMCMC_samples$summary[paste0("R[",15:rtConsts$N,", 1]"),1],
                  b = nimbleMCMC_samples$summary[paste0("R[",15:rtConsts$N,", 2]"),1],
                  c = nimbleMCMC_samples$summary[paste0("R[",15:rtConsts$N,", 3]"),1],
                  day = 1:201)

dat.long <- dat %>% pivot_longer(cols = a:c)

p1 <- ggplot(aes(x = day, y = value, color = name), data = dat.long) + geom_line() 

p1

dat <- data.frame(a = nimbleMCMC_samples$summary[paste0("lambda[",15:180,", 1, 1]"),1],
                  b = nimbleMCMC_samples$summary[paste0("lambda[",15:180,", 2, 1]"),1],
                  c = nimbleMCMC_samples$summary[paste0("lambda[",15:180,", 3, 1]"),1],
                  day = 1:166)

dat.long <- dat %>% pivot_longer(cols = a:c)

p2 <- ggplot(aes(x = day, y = value, color = name), data = dat.long) + geom_line()

p2

