data {
  int<lower=1> N;          // the number of observations. 
  int<lower=1> NW;          // the number of weeks. 
  int week_vec[N];          // the week vector 
  int<lower=1> J;          // the number of states
  int<lower=0> Y[N,J];     // observed cases. 
  matrix[N*J,J] P;         // transfer matrix, changes with time
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
  //real beta_mu;            // prior mu for beta
  //real beta_sd;            // prior sd for beta
  real sigma_shape;        // prior shape for sigma
  real sigma_scale;        // prior scale for sigma
  int<lower=0> Z;          // smoothing window size
  vector[J] init_cases;    // initial
  vector[J] fixed_sigma;   // fixed sigma
}

parameters {
  // set up for non-centered parameterization of BETA
  matrix[NW,J] beta_mu;              // time-state specific beta
  matrix<lower=0.01> [NW,J] beta_sd; // time-state specific beta
  matrix[NW,J] beta_raw;              // time-state specific beta
  //
  vector<lower=0.01>[J] xsigma; // state-specific st-dev
  //
  matrix[NW,J] logR_raw;             // time-state specific R values, in Log space
}

transformed parameters {
  matrix[NW,J] xbeta;      // time-state specific beta
  matrix[NW,J] logR;       // time-state specific R values, in Log space
  matrix[N,J] M;           // expected value of cases
  matrix[N,J] R;           // time and state specific R values, expressed normally
  matrix[J, J] RR;         // Diag R matrix
  
  //
  for(j in 1:J) {
    for(w_i in 1:NW) {
      
      xbeta[w_i, j] = beta_mu[w_i, j] + beta_sd[w_i, j]  * beta_raw[w_i, j];
      
      logR[w_i, j] = xbeta[w_i, j]  + xsigma[j]  * logR_raw[w_i, j];
      
    }
  }
  
  
  // get R in exp() space
  for(j in 1:J) {
    for(n in 1:N) {
      R[n, j] = exp(logR[week_vec[n], j]);
    }
  }
  
  // ------ EQUATION (11b) -------------
  // initialize M at t = 0
  for(j in 1:J) {
    M[1, j] = init_cases[j];
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

    // Create w_i vector
    matrix[tau_end, 1] w_i;
    w_i[1:tau_end, 1] = to_vector(W[1:tau_end]);

    // Calculate result
    int start_P = J*(n - 1) + 1;
    int end_P   = J*(n - 1) + J;
    M[n, ] = to_row_vector(P[start_P:end_P, ]' * RR * MM * w_i);
    
  }

}

model {
   
   
  // ------ EQUATION (11a) -------
  // priors and sample 
  xsigma ~ inv_gamma(sigma_shape, sigma_scale); 
  
  for(j in 1:J) {
      //  non-centered parms
       beta_mu[, j] ~ std_normal();
       beta_sd[, j] ~ std_normal(); // this should be something else
       beta_raw[, j] ~ std_normal();
       logR_raw[, j] ~ std_normal();
  }
  
  // ------ EQUATION (11c) -------
  for(j in 1:J) {
      // -- REAR WINDOW
      for(n in 1:N) {
        int min_row = max(1, n - Z + 1);
        Y[min_row:n, j] ~ poisson(M[n, j]);

      }
  }
    

}

