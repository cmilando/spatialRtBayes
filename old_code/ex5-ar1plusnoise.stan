data {
  int <lower=0> n;
  real x[n];
  real y[n];
}

parameters {
  real u_err[n];
  real <lower=0> sig;
  real <lower=0> tau;
  real<lower=0, upper=1> phi;
  real mu;
  real d;
}

transformed parameters{
  real u[n]; //Level
  real s[n];
  real z[n];
  u[1] = u_err[1];
  s[1] = 1;
  z[1] = u[1];
  for (t in 2:n) {
    u[t] = mu + phi*u[t-1] + tau*u_err[t];
    if (fabs(u[t])>=d){
      s[t] = 1;
    }else{
      s[t]=0;
    }
    z[t] = u[t]*s[t];
  }
}

model{
  u_err ~ normal(0,1);
  y     ~ normal(z,sig);
  mu    ~ normal(0,1);
  phi   ~ normal(0,1);
  tau   ~ gamma(1,1);
  sig   ~ gamma(1,1);
  d ~ uniform(0,fabs(mu/(1-phi))+3*tau/sqrt(1-phi^2));
}