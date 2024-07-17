data {
  int<lower=1> N_k;      // number of knots
  int<lower=1> N;        // number of points
  real Y[N];             // Y values
  real<lower=0.01> X[N]; // X values
  //real Kx[N_k];
  real Xmax;
  real Xmin;
}

parameters {
  real<lower=Xmin,upper=Xmax> Kx_raw[N_k - 2]; // x values
  real Ky[N_k]; // y values,
  real<lower=0.1> ysigma;
  real<lower=0.01> tau;
}

transformed parameters {
  real Y_calc[N];      //
  real y_intercepts[N_k - 1];
  real y_slopes[N_k - 1];
  matrix[2, 2] Xmat;
  vector[2] Ymat;
  vector[2] calc_beta;
  
  // spline X adjustment factor
  real<lower=Xmin,upper=Xmax> Kx[N_k]; // x values
  Kx[1] = Xmin;
  Kx[N_k] = Xmax;
  for (i in 2:(N_k - 1))
    Kx[i] = Kx[i-1] + Kx_raw[i - 1] * tau;
  //print("Kx = ", Kx);
  //print("Ky = ", Ky);
  
  // First, solve for Betas
  // spline_design = 1
  // you'd have to modify more if you wanted it to be different
  for(k in 1:(N_k - 1)) {
      
      //
      Xmat[1, 1] = 1;          // a
      Xmat[2, 1] = 1;          // c
      Xmat[1, 2] = Kx[k];      // b
      Xmat[2, 2] = Kx[k + 1];  // d 
      
      //
      Ymat[1] = Ky[k];
      Ymat[2] = Ky[k + 1];
      
      // Compute calc_beta = (Xmat' * Xmat)^-1 * Xmat' * Ymat
      calc_beta = mdivide_left_spd(Xmat' * Xmat, Xmat' * Ymat);
      
      //
      //betas[1, k] = calc_beta[1];
      y_intercepts[k] = calc_beta[1];
      y_slopes[k] = calc_beta[2];
      //betas[2, k] = calc_beta[2];
      //////////////////////

  }
  
  //print("y_intercepts = ", y_intercepts);
  //print("y_slopes = ", y_slopes);
  
  // Then plug in across all X to get Y_calc
  for(i in 1:N) {
    
    // for the bottom bracket is <=
    if(X[i] >= Kx[1] && X[i] <= Kx[2]) {
      //Y_calc[i] = betas[1, 1] + X[i] * betas[1, 2];
      Y_calc[i] = y_intercepts[1] + X[i] * y_slopes[1];
    }
    
    for(k in 2:(N_k - 1)) {
    
      // just < than now
      if(X[i] > Kx[k] && X[i] <= Kx[k+1]) {
        //Y_calc[i] = betas[1, k] + X[i] * betas[2, k];
        Y_calc[i] = y_intercepts[k] + X[i] * y_slopes[k];
      }
      
    }
    
  }
  
  //print("Y_calc = ", Y_calc);
      
}


model {
  
  // priors with truncation
  // https://mc-stan.org/docs/stan-users-guide/truncation-censoring.html#truncation.section
  Kx_raw ~ normal(0, 1) T[Xmin, Xmax];;
  tau ~ normal(0, 1) T[0,];
  Ky ~ normal(0, 1);
  ysigma ~ inv_gamma(1, 1);
  
  // likelihood
  Y ~ normal(Y_calc, ysigma);
  
}

generated quantities {
  
   real Y_out[N];
   
   // then this is same as in transformed parameters
  for(i in 1:N) {
    
    // for the bottom bracket is >=
    if(X[i] >= Kx[1] && X[i] <= Kx[2]) {
      //Y_out[i] = betas[1, 1] + X[i] * betas[1, 2];
      Y_out[i] = y_intercepts[1] + X[i] * y_slopes[1];
    }
    
    for(k in 2:(N_k - 1)) {
    
      // just > than now
      if(X[i] > Kx[k] && X[i] <= Kx[k+1]) {
        //Y_out[i] = betas[1, k] + X[i] * betas[2, k];
        Y_out[i] = y_intercepts[k] + X[i] * y_slopes[k];
      }
      
    }
    
  }
  
}
