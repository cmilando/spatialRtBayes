source("R/henrik.R")

set.seed(123L)
##
design <- cbind(sapply(1:3, function(x) runif(400)))
dim(design) # observed X
head(design)

## true betas
perfect_position <- c(1,2,3)

## the data
target <- design %*% perfect_position + rnorm(400)
dim(target)

plot(target)

#### test ####

###
# ** TARGET (Y) and DESIGN (X) ARE IMPLIED
log_linear <- function(position, sigma = 10L){
  sum((-2*sigma^(-2))*(target - design %*% position)^2)
}

# ** TARGET (Y) and DESIGN (X) ARE IMPLIED
partial_deriv_lin <- function(position,sigma =10L){
  -sigma^(-2) * t(design) %*% (-target + design %*% position)
}
###

### SO, in order to make this work, I need to define:
## DESIGN: 
## TARGET:
## log_posterior_density
## gradient

## these must be defined by the equations and the data


log_posterior_density = log_linear
gradient = partial_deriv_lin

##
N_ITER <- 20000
N_BURN_IN = 5000  # this avoids the non-convergence at the beginning
N_SKIP = 2        # this makes sure samples are independent
rows = seq(from = N_BURN_IN, to = N_ITER, by = N_SKIP)

linear_sample_l <- sample_NUT(c(5, 1, 10), stepsize = 0.15, iteration = N_ITER)

head(linear_sample_l)
head(target)

hist(linear_sample_l$X1[rows], breaks = 30)
quantile(linear_sample_l$X1[rows], probs = c(0.025, 0.5, 0.975))
abline(v = 1, col = 'red')

hist(linear_sample_l$X2[rows], breaks = 30)
quantile(linear_sample_l$X2[rows], probs = c(0.025, 0.5, 0.975))
abline(v = 2, col = 'red')

hist(linear_sample_l$X3[rows], breaks = 30)
quantile(linear_sample_l$X3[rows], probs = c(0.025, 0.5, 0.975))
abline(v = 3, col = 'red')

