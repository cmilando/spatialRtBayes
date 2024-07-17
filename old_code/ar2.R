library(rstan)
options(mc.cores = parallel::detectCores())

stan_string = "
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  for (n in 2:N) {
    y[n] ~ normal(alpha + beta * y[n-1], sigma);
  }
}
"

set.seed(123)
N <- 100
alpha_true <- 0.5
beta_true <- 0.8
sigma_true <- 1.0
y <- numeric(N)
y[1] <- rnorm(1, 0, sigma_true)
for (n in 2:N) {
  y[n] <- rnorm(1, alpha_true + beta_true * y[n-1], sigma_true)
}
plot(y, type = 'l')

# source("R/00_generate_data.R")

# Prepare data list for STAN
stan_data <- list(N = N, y = y)

# Fit the model using STAN
fit <- stan(model_code = stan_string, data = stan_data, chains = 4, seed = 123)

# Print the summary of the model fit

out <- rstan::extract(fit)


y_pred <- numeric(N)
y_pred[1] <- y[1]
for (n in 2:N) {
  y_pred[n] <- mean(out$alpha) + mean(out$beta) * y_pred[n-1]
}

lines(y_pred, type = 'l', col = 'blue')



Rmatrix

stan_data = list(
  N = as.integer(nrow(R_rev) - 1),
  y = R_rev[-1, 1]
)

R_rev[-1, 1]
y = R_rev[-1, 1]
out <- rstan::stan(model_code = stan_string, data = stan_data)

out_ex <- rstan::extract(out)

dim(out_ex$alpha)
dim(out_ex$beta)

mean_alpha = mean(out_ex$alpha[2001:4000])
mean_beta = mean(out_ex$beta[2001:4000])
ypred = 
for()