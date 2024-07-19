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
si_shape <- 2       ## shape parameter for serial interval assuming gamma distribution
si_rate  <- 0.5     ## rate parameter for serial interval assuming gamma distribution
si_t_max <- 14      ## maximum number of days with non-zero probability for serial interval

source('00_generate_data.R')

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
  init_cases = c(10, 10, 10)
)

# -------------------------------------------------------------
# run the STAN code, takes ~ 10min
N_ITER = 4000
N_WARMUP = 2000
N_CHAINS = 4

# compile first
model1 <- stan_model(file="../src/stan_sliding_v4nc1.stan",
                     model_name = "spatialRt")

# then run
m_hier <- sampling(object = model1,
               data = stan_data, 
               chains = N_CHAINS,
               warmup = N_WARMUP,
               iter = N_ITER + N_WARMUP,
               control = list(max_treedepth = 20,
                              adapt_delta = 0.99))

# extract output
# saveRDS(m_hier, 'stan_out_weekly_nc.RDS')
# 
# shinystan::launch_shinystan(m_hier)




