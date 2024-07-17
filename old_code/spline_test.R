# ttps://github.com/milkha/Splines_in_Stan/blob/master/splines_in_stan.Rmd

library(rstan)
library(splines)
set.seed(1234)
num_knots <- 10 # true number of knots
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
X <- seq(from=-10, to=10, by=.1)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
num_data <- as.integer(length(X))
a0 <- 0.2
a <- rnorm(num_basis, 0, 1)
B_true <- t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
Y_true <- as.vector(a0*X + a%*%B_true)
Y <- Y_true + rnorm(length(X), 0, 0.2)

num_knots <- 10; # number of knots for fitting
num_basis <- num_knots + spline_degree - 1
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))
rstan_options(auto_write = TRUE);

options(mc.cores = 4);

#spline_model<-stan_model("b_spline.stan")
spline_penalized_model<-stan_model("b_spline_penalized.stan")

#fit_spline<-sampling(spline_model,iter=500,control = list(adapt_delta=0.95))
fit_spline_penalized<-sampling(spline_penalized_model,
                               data = list(
                                 num_data = num_data,             #// number of data points
                                 num_knots = num_knots,            #// num of knots
                                 knots = knots,  #// the sequence of knots
                                 spline_degree = spline_degree,  #// the degree of spline (is equal to order - 1)
                                 Y = Y,
                                 X = X
                               ),
                               #iter=500,
                               control = list(adapt_delta=0.95))

# oo <- rstan::extract(fit_spline_penalized)
# dim(oo$a_raw)
# dim(oo$a0)
# 
# dim(oo$a_raw
#     )
# dim(oo$Y_hat)
# dim(oo$a_raw)



ff<-extract(fit_spline_penalized)
Y_hat_med_pen <- array(NA, length(Y))
Y_hat_ub_pen <- array(NA, length(Y))
Y_hat_lb_pen <- array(NA, length(Y))
for (i in 1:length(Y)) {
  Y_hat_med_pen[i] <- median(ff$Y_hat[,i]);
  Y_hat_lb_pen[i] <- quantile(ff$Y_hat[,i],probs = 0.25)
  Y_hat_ub_pen[i] <- quantile(ff$Y_hat[,i],probs = 0.75)
}
plot(X,Y, col="azure4")
#lines(X, Y_hat_med, col="Red", lw=2, lty=1)
lines(X, Y_hat_med_pen, col="blue",lw=2, lty=1)
legend(-1,-2,0,legend=c("fit without smoothing prior", "fit with smoothing prior"), col=c("red", "blue"), lw=4, box.lty=0, lty=1, border="white")

# soo ... solve for the exact 


## just do linear piecewise
## wait no you do know the knots


library(splines)

## right, splines are just a way of breaking up X
## they don't 

## so no, you construct a spline from points
## can do this with matrices


