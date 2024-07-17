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
  vector[Q] theta_raw;  // raw error coefficients
  real<lower=0> sigma_theta; // scale for theta
}
transformed parameters {
  vector[Q] theta;      // error coeff, lag -t
  theta = theta_raw * sigma_theta; // scale the raw theta

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
  mu ~ normal(0, 2); // informative prior
  theta_raw ~ normal(0, 1); // standard normal prior for raw coefficients
  sigma_theta ~ cauchy(0, 2.5); // prior for the scale
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
  vector[N] eta;
  vector[N] epsilon_pred;
  for (n in 1:N) {
    epsilon_pred[n] = y[n] - mu;
    for (q in 1:min(n - 1, Q)) {
      epsilon_pred[n] = epsilon_pred[n] - theta[q] * epsilon_pred[n - q];
    }
  }
  for (n in 1:N) {
    eta[n] = mu;
    for (q in 1:min(n - 1, Q)) {
      eta[n] = eta[n] + theta[q] * epsilon_pred[n - q];
    }
    y_pred[n] = normal_rng(eta[n], sigma);
  }
}
"

# Generate some example data
set.seed(123)
N <- 100
Q <- 2
mu <- 0.5
theta <- c(0.3, 0.2)
sigma <- 1
epsilon <- numeric(N)
y <- numeric(N)
epsilon[1:Q] <- rnorm(Q, 0, sigma)
for (n in 1:Q) {
  y[n] <- rnorm(1, mu, sigma)
}
for (n in (Q + 1):N) {
  epsilon[n] <- rnorm(1, 0, sigma)
  y[n] <- rnorm(1, mu + sum(theta * epsilon[(n - Q):(n - 1)]), sigma)
}

# Prepare data for STAN
stan_data <- list(Q = Q, N = N, y = y)

# Fit the model
fit <- stan(model_code = stan_model_code, data = stan_data, 
            iter = 2000, chains = 4)

# Print the results
print(fit)

# Extract the generated quantities (predictions)
y_pred <- extract(fit)$y_pred

# Plot the original data and the predictions
plot(1:N, y, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "y", main = "Observed vs Predicted y")
lines(1:N, apply(y_pred, 2, mean), col = "red", lwd = 2)
legend("topleft", legend = c("Observed", "Predicted"), col = c("blue", "red"), lwd = 2)
