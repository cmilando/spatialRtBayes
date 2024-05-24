data {
  int<lower=1> N_k;      // number of knots
  int<lower=1> N;        // number of points
  real Y[N];             // Y values
  real<lower=0.01> X[N]; // X values
  real Kx[N_k];
  real Xmax;
  real Xmin;
}

parameters {
  //real<lower=Xmin,upper=Xmax> Kx[N_k]; // x values
  real Ky[N_k]; // y values,
  real<lower=0.1> ysigma;
}

transformed parameters {
  real Y_calc[N];      //
  real y_intercepts[N_k - 1];
  real y_slopes[N_k - 1];
  //real betas[2, N_k - 1];  // betas. 2 because this is a linear 
  matrix[2, 2] Xmat;
  //matrix[2, 2] XtX;
  vector[2] Ymat;
  vector[2] calc_beta;
  
  // First, solve for Betas
  // spline_design = 1
  // you'd have to modify more if you wanted it to be different
  for(k in 1:(N_k - 1)) {
      
      //print("kx = ", Kx);
      
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
      
}


model {
  
  // priors
  //Kx ~ normal(10, 5);
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
