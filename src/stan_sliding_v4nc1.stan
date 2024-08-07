data {
  int<lower=1> N;          // the number of observations. 
  int<lower=1> NW;         // the number of weeks. 
  int week_vec[N];         // the week vector 
  int<lower=1> J;          // the number of regions
  int<lower=0> Y[N,J];     // observed cases. 
  matrix[N*J,J] P;         // transfer matrix, changes with time
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
  vector[J] init_cases;    // initial cases
}

parameters {
  
  // CENTRAL RT
  real<lower=0.00> sigma_logRt_central;  // how do we scale it
  vector[NW] logRt_central_error;       // draws from the st normal
  real logRt_central_intercept;

  // DEVIATION
  vector<lower=0.00>[J] sigma_logRt;    // central, could be region specific
  matrix[NW, J] logRt_error; // for reach region, whats the avg deviation
  
}

transformed parameters {
  //
  vector[NW] logR_central; 
  matrix[NW, J] logR; 
  
  //
  matrix[N,J] R;           // time and region specific R values, expressed normally
  //
  matrix[N,J] M;           // expected value of cases
  matrix[J, J] RR;         // Diag R matrix
  
  // ------ CALCULATE R(t) -------------
  
  // CENTRAL RT -- TIME VARYING
  logR_central[1] = logRt_central_intercept;
  for(w_i in 2:NW) {
    logR_central[w_i] = logR_central[w_i-1] + sigma_logRt_central * logRt_central_error[w_i];
  }
  
  // DEVIATION -- SPACE VARYING
  for(j in 1:J) {
    for(w_i in 1:NW) {
      logR[w_i, j] = logR_central[w_i] + sigma_logRt[j] * logRt_error[w_i, j];
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
   
  // ------ EQUATION (11a - modified to be weekly) -------
  // Central R(t)
  sigma_logRt_central ~ normal(0, 0.1);   // week to week variation in Rt
  logRt_central_error ~ std_normal();     // non-centered parameterization
  logRt_central_intercept ~ normal(1, 1); // starting point

  sigma_logRt ~ normal(0, 0.1);   // week to week variation in Rt
  to_vector(logRt_error) ~ std_normal();     // non-centered parameterization
  
  // ------ EQUATION (11c) -------
  for(j in 1:J) {
      for(n in 1:N) {
        Y[n, j] ~ poisson(M[n, j]);
      }
  }
    

}

