library(rstan)
library(parallel)
library(invgamma)
library(tidyverse)
library(matrixStats)
library(lemon)
library(shinystan)

num_cores <- parallel::detectCores()
options(mc.cores = min(4, num_cores))
rstan_options(auto_write = TRUE)

# -------------------------------------------------------------
#### INPUT
window_length <- 1  ## smoothing window for estimaiton
si_shape <- 2       ## shape parameter for serial interval assuming gamma distribution
si_rate  <- 0.5     ## rate parameter for serial interval assuming gamma distribution
si_t_max <- 14      ## maximum number of days with non-zero probability for serial interval
beta_mu  <- 0.      ## prior mean for Rj(t) ~ N
beta_sd  <- 1.      ## prior sd   for Rj(t) ~ N
sigma_shape <- 2.   ## prior for sigma_j assuming inverse gamma
sigma_scale <- 1.   ## prior for sigma_j assuming inverse gamma

source('R/00_generate_data.R')

Y <- as.matrix(N)
for(i in 1:nrow(Y)) {
  for(j in 1:ncol(Y)) {
    Y[i,j] <- as.integer(Y[i,j])
  }
}
NN <- as.integer(nrow(Y))

# make daily
P_list <- lapply(1:NN, function(x) P)
P <- do.call(rbind, P_list)

#
week_vec = as.integer(ceiling((1:NN)/7))
n_weeks = as.integer(length(unique(week_vec)))

# Data list for Stan
stan_data <- list(
  N = NN,
  NW = n_weeks,
  week_vec = week_vec,
  J = J,
  Y = Y,
  P = P,
  S = S,
  W = w,
  beta_mu = beta_mu,
  beta_sd = beta_sd,
  sigma_shape = sigma_shape, 
  sigma_scale = sigma_scale,
  Z = as.integer(window_length),
  init_cases = c(10, 10, 10),
  fixed_sigma = c(0.1, 0.1, 0.1)
)

# -------------------------------------------------------------
# run the STAN code, takes ~ 10min
N_ITER = 1000
N_WARMUP = 1000
N_CHAINS = 4
m_hier <- stan(file="stan_sliding_v4.stan",
               data = stan_data, 
               chains = N_CHAINS,
               warmup = N_WARMUP,
               iter = N_ITER + N_WARMUP)

# extract output
# saveRDS(m_hier, 'stan_out_weekly.RDS')
# 
# shinystan::launch_shinystan(m_hier)




