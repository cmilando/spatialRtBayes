
library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(gridExtra)
library(parallel)
m <- 3 # number of regions
# pop <- rnorm(m, 50000, 10000) %>% round()
# initial_inf_rate <- rbeta(m, 1, 500)
# initial_inf_num <- round(pop*initial_inf_rate)

# initial_inf_num <- c(10, 10, 10)
# N_day <- 200


w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)


P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

N_simulated_data <- read.csv("example_data.csv")


stan_data<-list()
stan_data$N  <- nrow(N_simulated_data)
stan_data$casesa <- N_simulated_data$a
# stan_data$casesa <- rep(2, 182)
stan_data$casesb <- N_simulated_data$b
# stan_data$casesb <- rep(2, 182)
stan_data$casesc <- N_simulated_data$c
# stan_data$N2 <- 60

stan_data$cases <- matrix(c(N_simulated_data$a, N_simulated_data$b, N_simulated_data$c), 
                          nrow = length(N_simulated_data$a))

head(stan_data$cases)

stan_data$K <- dim(N_simulated_data)[2]

stan_data$P <- P

# 
# stan_data$logra <- log(R_est_sim_data[, 1])
# stan_data$logrb <- log(R_est_sim_data[, 2])
# stan_data$logrc <- log(R_est_sim_data[, 3])



w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)
w <- c(w, rep(0, nrow(N_simulated_data) - length(w)))

stan_data$SI <- w

# saveRDS(N_simulated, "/rprojectnb2/zwzhou-thesis/rstan/simulated_data_exo.rds")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## no imported infections

stan_model1 <- stan_model(model_code = 
                            "// Stan model for simple linear regression

data {
  int<lower=1> N; // days of observed data
  int<lower=1> K;
  int cases[N, K];
  matrix[K, K] P;
  real SI[N]; // fixed SI using empirical data
}
transformed data {
  vector[N] SI_rev; // SI in reverse order
  for(i in 1:N){
    SI_rev[i] = SI[N-i+1];
  }
}
parameters {
  matrix[N, K] weekly_effect;
  matrix<lower=0, upper=1>[N, K] weekly_rho;
  vector<lower=0>[K] weekly_sd;
}

transformed parameters {
  matrix[N, K] prediction=rep_matrix(1e-5, N, K);
  matrix[N, K] convolution;
  matrix<lower=0>[N, K] Rt;
  {
    Rt = exp(weekly_effect);
    for (i in 2:N) {
      convolution[i, ] = tail(SI_rev, i-1)'*prediction[1:(i-1),];
      prediction[i, ] = prediction[i, ] + (convolution[i, ] .* Rt[i, ])*P;
    }
  }
}
model {
  weekly_sd ~ normal(0,1);
  for (j in 1:K){
    weekly_rho[,j] ~ normal(0, 0.5);
    weekly_effect[1:N, j] ~ normal(weekly_rho[1:N, j] , weekly_sd[j]);
  }
  
  for (i in 1:N) {
    for (j in 1:K){
      cases[i:min(i+8, N), j] ~ poisson(prediction[i, j]);
    }
  }
  
}"
)

fit <- sampling(stan_model1,data=stan_data)
