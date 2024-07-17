library(MASS)

## ok lets make a test case


# Beta is J x T
Jx = 2
Tx = 1000 ## the higher this is, the better the estimates are. 
## makes sense. more data!

# 
true_beta  <- matrix(NA, nrow = Tx, ncol = Jx)
for(j in 1:Jx) {
  for(t in 1:Tx) {
    true_beta[t, j] <- t*j + 0*runif(1, -5, 5)
  }
}

plot(true_beta[1:t, 1], type='l', col = 'blue')
if(Jx > 1) points(true_beta[1:t, 2], type='l', col = 'red')
if(Jx > 2) points(true_beta[1:t, 3], type='l', col = 'green')

true_sigma <- matrix(c(2, 4), nrow = 1, ncol = Jx)

# so calculate observed R[t,j]
obs_Y  <- matrix(NA, nrow = Tx, ncol = Jx)
for(j in 1:Jx) {
  for(t in 1:Tx) {
    obs_Y[t, j] <- rnorm(1, mean = true_beta[t, j], true_sigma[1, j])
  }
}

plot(obs_Y[1:t, 1], col = 'blue')
if(Jx > 1) points(obs_Y[1:t, 2], col = 'red')
if(Jx > 2) points(obs_Y[1:t, 3], col = 'green')

# ******
# Ok so now,
# basically staY with uninformed estimates of beta and sigma and iterate
set.seed(1)
N_ITER <- 20000
MAX_SIGMA <- 5

est_beta  <- matrix(rnorm(Tx * Jx), nrow = Tx, ncol = Jx)
est_beta

# *** for testing
est_beta  <- true_beta
# ***

# est_sigma <- matrix(sapply(1:Jx, function(j) 
#   rgamma(n = 1, scale = 1, shape = 1)), nrow = 1, ncol = Jx)

est_sigma_array <- array(runif(N_ITER * Jx, 0, MAX_SIGMA), dim = c(N_ITER, Jx))
dim(est_sigma_array)

# *******
# wait first question -- sigma has to be positive ...
# but its log
stopifnot(all(est_sigma_array > 0))
# ********

## loglikelihood of normal
## https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood
# get_ll <- function(this_sd, mu_vec, y_vec) {
#     ##sigma2 = this_sd^2
#   sum(sapply(1:length(mu_vec), function(t) {
#     this_y = y_vec[t]
#     this_mu = mu_vec[t]
#     dnorm(this_y, this_mu, sd = this_sd, log = TRUE)
#     ##-.5 * log(2 * pi) -.5 * log(sigma2) - (1/(2 * sigma2)) * (this_y - this_mu)^2
#   }))
# }

get_ll_all_j <- function(sigma_vec, mu_vec, y_vec) {
  o <- 0
  for(j in 1:ncol(mu_vec)) {
    this_sd = sigma_vec[j]
    for(t in 1:nrow(mu_vec)) {
      this_y = y_vec[t, j]
      this_mu = mu_vec[t, j]
      o <- o + dnorm(this_y, this_mu, sd = this_sd, log = TRUE)
    }
  }
  return(o)
}

## *******************************************************
## if you could sample directly, you could just do this
##

test_j = 2
mu = true_beta[, test_j]
y = obs_Y[, test_j]

sd_seq = sort(unique(c(seq(from = 0.1, to = MAX_SIGMA, length.out = 1000),
                       seq(from = 0.1, to = MAX_SIGMA, by = 1))))
#hist(sd_seq)
### why is this influenced by the number of points in sd_seq?

xx <- sapply(sd_seq, function(s) get_ll(s, mu_vec = mu, y_vec = y))

plot(x = sd_seq, y = xx, type = 'l')

true_sigma[, test_j]
sd_seq[which(xx == max(xx))]

## *******************************************************
## But you can't, so iterate

## This is called SLICE sampling


for (i in 2:N_ITER) {
  
  if(i %% 1000 == 0) cat(i, "\t")
  
  old_sigma = est_sigma_array[i - 1, 1:Jx]
  
  for (j in 1:Jx) {
    
    # so first get a new estimate of J
    # print(j)
    
    # clone est_sigma
    new_sigma <- old_sigma
    
    # test_sigma
    new_sigma[j] <- runif(1, 0, MAX_SIGMA)
    
    # now calculate lvl 0  M-H
    
    # its, are the updated est_sigma and est_beta more likely to be 
    # the parameters that estimate obs_Y
    # than the previous version of est_sigma and est_beta.
    
    # new_Y is not explicitly calculated, 
    # instead, you use the likelihood function to compare the
    # observations against 
    # expected value given the parameters 
    
    ## at its best, LL_NEW is the negative cloest to zero
    ## LL_OLD: a negative number less close to zero
    LL_NEW <- get_ll_all_j(new_sigma, mu_vec = est_beta, y_vec = obs_Y)
    LL_OLD <- get_ll_all_j(old_sigma, mu_vec = est_beta, y_vec = obs_Y)
    
    # ********
    # ok now do the likelihood ratio test
    ratio = exp(LL_NEW - LL_OLD)
    
    ## this covers you if ratio = Inf
    ratio <- min(1, ratio)
    
    # get a random U(0,1) variable to compare against
    compU <- runif(1, 0, 1)
    
    ## if u < Ratio, details details ...
    # so, as LL_NEW gets 
    decision <- compU < ratio
    
    ## if decision == TRUE, then update params
    if(decision) {
      
      print(LL_NEW)
      
      ## update parameter = NEW PARAMS
      est_sigma_array[i, j] = new_sigma[j]
      
      ### ********
      ### THEN DO THE LEVEL BELOW IF AND ONLY IF YOU ACCEPTED ABOVE
      ### ********
      
      ## since these are all normal, you could try and do as a block
      
      
      
    } else {
      
      ## set parameter = OLD PARAMS
      est_sigma_array[i, j] = old_sigma[j]
      
    }
    
    # set old_sigma for this j only if it needs to be updated
    old_sigma[j] = est_sigma_array[i, j]
    
  }
  
}

N_BURN_IN = 5000  # this avoids the non-convergence at the beginning
N_SKIP = 2        # this makes sure samples are independent
rows = seq(from = N_BURN_IN, to = N_ITER, by = N_SKIP)

quantile(est_sigma_array[rows, 1], probs = c(0.025, .5, 0.975))
quantile(est_sigma_array[rows, 2], probs = c(0.025, .5, 0.975))

plot(est_sigma_array[rows, 1], main = paste("true sigma = ", true_sigma[1]))
abline(a = true_sigma[1], b = 0, col = 'red')
plot(est_sigma_array[rows, 2], main = paste("true sigma = ", true_sigma[2]))
abline(a = true_sigma[2], b = 0, col = 'red')

hist(est_sigma_array[rows, 1], main = paste("true sigma = ", true_sigma[1]))
hist(est_sigma_array[rows, 2], main = paste("true sigma = ", true_sigma[2]))

