data {
  int<lower=1> N;          // the number of observations. 
  int<lower=1> J;          // the number of states
  int<lower=0> Y[N,J];     // observed cases. 
  matrix[N*J,J] P;         // transfer matrix, changes with time
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
  real beta_mu;            // prior mu for beta
  real beta_sd;            // prior sd for beta
  real sigma_shape;        // prior shape for sigma
  real sigma_scale;        // prior scale for sigma
  int<lower=0,upper=14> Z; // sliding window size, 0 to 14
}

parameters {
  // sd, beta, and logR
  vector<lower=0.01>[J] xsigma; // state-specific st-dev
  matrix[N,J] xbeta;            // time-state specific beta
  matrix[N,J] logR;             // time-state specific R values, in Log space
}

transformed parameters {
  // M and R go here
  matrix[N,J] M;           // expected value of cases
  matrix[N,J] R;           // time and state specific R values, expressed normally
  //vector[J] inner_vec;       // gets the vectors of the m_w sums
  //vector[J] outer_vec; 
  matrix[J, J] RR;

  // get R in exp() space
  for(j in 1:J) {
    R[1:(N), j] = exp(logR[1:(N), j]);
  }
  
  // ------ EQUATION (11b) -------------
  // initialize M at t = 0
  for(j in 1:J) {
    M[1, j] = 1e-5;
    R[1, j] = 1e-5;
  }
  
  // calculate m[t,j] for t >= 1
  for(n in 2:N) {
    
    // calculate the inner part first
    int tau_end = min(S, n - 1);
    
    // Create diagonal matrix RR
    RR = diag_matrix(to_vector(R[n, ]));
    
    // Create MM based on t value
    // MM is m(t - 1), ..., m(1)
    // where rows are regions
    // and columns are time
    matrix[J, tau_end] MM;
    for(tt in 1:tau_end) {
      MM[, tt] = to_vector(M[n - tt, ]);
    }

    // Create WW vector
    matrix[tau_end, 1] WW;
    WW[1:tau_end, 1] = to_vector(W[1:tau_end]);

    // Calculate result
    int start_P = J*(n - 1) + 1;
    int end_P   = J*(n - 1) + J;
    M[n, ] = to_row_vector(P[start_P:end_P, ]' * RR * MM * WW);
    
    /*
    for (jx in 1:J) {
      inner_vec[jx] = 0;
      for(tau in 1:tau_end) {
        inner_vec[jx] = inner_vec[jx] + R[n,jx] * M[n - tau,jx] * W[tau];
      }
    }
    
    // then the outer loop
    for(j in 1:J) {
      outer_vec[j] = 0;
      for(jx in 1:J) {
        outer_vec[j] = outer_vec[j] + P[jx, j] * inner_vec[jx];
      }
      M[n, j] = outer_vec[j];
    }
    */
    
  }

}

model {
   
   
  // ------ EQUATION (11a) -------
  // priors and sample 
  xsigma ~ inv_gamma(sigma_shape, sigma_scale); 

  for(j in 1:J) {
      
      // weak prior on beta 
      xbeta[, j] ~ normal(beta_mu, beta_sd); 
      
      // sample logR
      logR[1:N, j] ~ normal(xbeta[1:N, j], xsigma[j]);
      
  }
  
  // ------ EQUATION (11c) -------
  for(j in 1:J) {
      // -- FORWARD WINDOW
      for(n in 1:N) {
        int max_row = min(N, n + Z - 1);
        target += poisson_lpmf(Y[n:max_row, j] | M[n, j]);
      }
  }
    

}

