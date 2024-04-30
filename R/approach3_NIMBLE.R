library(nimble)
library(tidyr)
library(ggplot2)
library(tidyverse)

incidence_data <- read_csv("https://raw.githubusercontent.com/zwzhou-biostat/hrt/main/data/example_data.csv")

incidence_data <- read.csv("../raw_data/example_data.csv")

usethis::use_data(incidence_data, overwrite = T)


#incidence_data <- read.csv("example_data.csv")
w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)
w <- rev(c(w, rep(0, nrow(incidence_data) - length(w))))

P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3, byrow = T)


N = length(incidence_data$a)
K = dim(P)[1]

L <- 8

N
K
L

incidence_data_expand <- array(NA, dim = c(N-L, K, L))
head(incidence_data_expand)

for(i in 1:(N-L)){
  for(k in 1:K){
    for(j in 1:L){
      incidence_data_expand[i, k, j] <- incidence_data[i+j-1, k]
    }
  }
}

###########################################################
rtCode <- nimbleCode({ 
  
  # Define relationships between nodes
  for (i in 2:(N-L)){
    for(k in 1:K){
      convolution[i-1, k] <- sum((lambda[1:(i-1), k, 1]*(w[(N-i+2):N])))
    }
    for(k in 1:K){
      for(j in 1:L){
        lambda[i, k, j] <- sum(t(P[1:K, k])*R[i, 1:K]*convolution[i-1, 1:K])
      }
    }　
  }
  
  for(i in 1:(N-L)){
    for(k in 1:K){
      for(j in 1:L){
        x[i, k, j] ~ dpois(lambda[i, k, j])
      }
    }
  }
  
  # Set priors
  for(k in 1:K){
    for(j in 1:L){
      lambda[1, k, j] <- 10
    }
    beta[k] ~ dnorm(1, 0.5)
    for(i in 1:N){
      alpha[i, k] ~ dnorm(0, 1)
      theta[i, k] ~ dnorm(alpha[i, k], beta[k])
      R[i, k] <- exp(theta[i, k])
    }
  }
  
})

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
                                 data = rtData, 
                                 inits = rtInits,
                                 setSeed = 1,
                                 nburnin = 5000, niter = 10000,
                                 monitors = c("R", "lambda"),
                                 summary = T, samples = T)
###########################################################

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

