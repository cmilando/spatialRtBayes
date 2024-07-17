# https://hedibert.org/wp-content/uploads/2021/02/stan-rstan-examples.html#example-5-ar1-plus-noise-model
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# simulating some data
set.seed(123)
n   = 200
mu  = 0.05
phi = 0.95
sig = 1.0
tau = 0.7
d   = 0.0001
mean = mu/(1-phi) 
sd   = tau/sqrt(1-phi^2)
x = numeric(n)
x[1]   = rnorm(1 ,mu/(1-phi),tau/sqrt(1-phi^2))
for(t in 2:n)
  x[t] = rnorm(1,mu + phi*x[t-1] , tau)
z = x
z[abs(z) < d] = 0
y  = rnorm(n, z, sig)

plot(x, type = 'l')
lines(y, type = 'l', col = 'red')


# trying to predict Y

# Running stan code
data <- list(n=n,x=x,y=y)

fit <- stan(file='ex5-ar1plusnoise.stan',data=data)

la <- extract(fit, permuted = TRUE) # return a list of arrays 

dim(la$z)
zz <- colMeans(la$z)
length(zz)
plot(zz)
lines(y,col='red')

acfs = matrix(0,51,n)
for (t in 1:n)
  acfs[,t] = acf(la$u[,t],lag.max=50,plot=FALSE)$acf
qx = t(apply(la$u,2,quantile,c(0.05,0.5,0.95)))
qz = t(apply(la$z,2,quantile,c(0.05,0.5,0.95)))


par(mfrow=c(1,1))
ts.plot(qx,ylab="states")
lines(qz[,1],col=2)
lines(qz[,2],col=2)
lines(qz[,3],col=2)
lines(x,col=3)
points(z,col=4)
legend("topright",legend=c("Posterior quantiles",expression(x[t])),col=1:2,lty=1)
