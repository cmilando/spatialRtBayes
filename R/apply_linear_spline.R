# this gets your x and y
KNOTS <- 5
L_OUT <- 20
source("R/linear_spline.R")
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
head(ll_df, 21)

dim(ll_df)

stan_data = list(
  N_k = as.integer(KNOTS),
  N = as.integer(nrow(ll_df)),
  Y = ll_df$new_y,
  X = ll_df$new_x,
  Kx = as.integer(x),
  Xmax = max(ll_df$new_x),
  Xmin = min(ll_df$new_x)
)

m1 <- stan_model("R/linearInterSpline.stan")

m1 <- stan_model("linear_spline_backup.stan")

fit1 <- sampling(m1, data = stan_data, iter = 5000,
                 control = list(max_treedepth = 30, 
                                adapt_delta = 0.9999))

x
## so its not a question of going up or down
## it works for 3 KNOTS but starts to crap out at 4
## which I think just means it needs more time !!!!!
## NICE JOB!!!!
## this could also be a function of the number of X points that you have
## huh, doesn't seem to work well for L_OUT = 3
## but does work for L_OUT = 100
### OOOF NO IT WAS A PROBLEM WITH DUPLICATES IN L, which won't happen for me

## No the number of iterations does scale with the number of knots
## and the max_treedepth


# sets that work:
# 4 knots, L = 50, iter = 5000, max_depth = 30, adapt_delta = 0.9999
## but ... not every time it seems .... huh

# 5 knots, L = 50, iter = 5000, max_depth = 50, adapt_delta = 0.9999 ## not quite but close

# doesn't work:
# 5 knots, L = 50, iter = 5000, max_depth = 30, adapt_delta = 0.9999

oo <- rstan::extract(fit1)

#dim(oo$betas)

Y_out <- t(apply(oo$Y_out, 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}))

dim(Y_out)
dim(ll_df)
head(ll_df$new_y)
head(Y_out)

lines(x = ll_df$new_x, y = Y_out[, 2], col = 'blue')
lines(x = ll_df$new_x, y = Y_out[, 1], col = 'green')
lines(x = ll_df$new_x, y = Y_out[, 3], col = 'green')

t(apply(oo$Ky, 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}))


t(apply(oo$Kx, 2, function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}))

