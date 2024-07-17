# Load necessary library
library(rstan)

# Define the STAN model as a string
stan_model_code <- "
data {
  int<lower=0> Q;       // num previous noise terms
  int<lower=3> N;       // num observations
  vector[N] y;          // observation at time t
}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  vector[Q] theta;      // error coeff, lag -t
}
transformed parameters {
  vector[N] epsilon;    // error term at time t
  for (n in 1:N) {
    epsilon[n] = y[n] - mu;
    for (q in 1:min(n - 1, Q)) {
      epsilon[n] = epsilon[n] - theta[q] * epsilon[n - q];
    }
  }
}
model {
  vector[N] eta;
  mu ~ cauchy(0, 2.5);
  theta ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 2.5);
  for (n in 1:N) {
    eta[n] = mu;
    for (q in 1:min(n - 1, Q)) {
      eta[n] = eta[n] + theta[q] * epsilon[n - q];
    }
  }
  y ~ normal(eta, sigma);
}
generated quantities {
  vector[N] y_pred;
  for (n in 1:N) {
    real eta_n = mu;
    for (q in 1:min(n - 1, Q)) {
      eta_n = eta_n + theta[q] * epsilon[n - q];
    }
    y_pred[n] = normal_rng(eta_n, sigma);
  }
}
"

# Generate sine wave data with mild noise
set.seed(123)
N <- 200
x <- rep(seq(0, 2 * pi, length.out = N/4), 4)
y <- sin(x) + rnorm(N, 0,1) # sine wave with mild noise
Q <- 10
plot(y, type = "l")

# Prepare data for STAN
stan_data <- list(Q = Q, N = N, y = y)

# Fit the model
fit <- stan(model_code = stan_model_code, data = stan_data, 
            chains = 4, iter = 20000)

library(shinystan)
shinystan::launch_shinystan(fit)
# Print the results
#print(fit)

# Extract the generated quantities (predictions)
y_pred <- extract(fit)$y_pred

# Plot the original data and the predictions
plot(1:N, y, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "y", main = "Observed vs Predicted y")
lines(1:N, apply(y_pred, 2, mean), col = "red", lwd = 2)
lines(1:N, apply(y_pred, 2, quantile, probs = 0.025), col = "red", lwd = 1)
lines(1:N, apply(y_pred, 2, quantile, probs = 0.975), col = "red", lwd = 1)
legend("topleft", legend = c("Observed", "Predicted"), col = c("blue", "red"), lwd = 2)
