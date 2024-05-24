data {
  int<lower=1> N;       // the number of observations. 
  int<lower=1> J;       // the number of states
  int<lower=0> Y[N,J];  // observed cases. first row padded 
  matrix[N*J,J] P;      // transfer matrix, in this example static in time
  int<lower=1> S;       // length of serial interval
  vector[S] W;          // serial interval
}

parameters {
  // sd, beta, and logR
  //vector<lower=0.01>[J] xsigma; // state-specific st-dev
  matrix[N,J] xbeta;            // time-state specific beta
  matrix[N,J] logR;             // time-state specific R values, in Log space
}

transformed parameters {
  // M and R actually go here
  matrix[N,J] M;           // expected value of cases
  matrix[N,J] R;           // time and state specific R values, expressed normally
  vector[J] sum_m_w;       // gets the vectors of the m_w sums
  real prod_R_p_sum;       // another inner variable
  row_vector[J] r_times_p; // gets the R*p product

  
  // get R in exp() space
  for(j in 1:J) {
    R[1:(N), j] = exp(logR[1:(N), j]);
  }
  
  // ------ EQUATION (11b) -------------
  // initialize M at t = 0
  for(j in 1:J) {
      M[1:(N),j] = rep_vector(1e-5,N);
  }
  
  // calculate m[t,j] for t >= 1
  for(n in 2:N) {
    
    // for first get the sum of the m * w vector
    // this can be a dot_prodcut of the last t-Tau timepoints of w and m
    // for m[(n - S):(n - 1)]
    // but of course, n-tau can never be lower than 1
    for (jx in 1:J) {
      int start_index = max(1, n - S + 1);
      int end_index = n - 1;
      sum_m_w[jx] = dot_product(M[start_index:end_index, jx], 
                                W[1:(end_index - start_index + 1)]); 
    }
    
    // so this is the outer loop, and there is actually
    // and inner loop as well
    for(j in 1:J) {
      
      // you need yet another inner loop, because Pj is dependent on this
      // r_times_p is a row_vector
      int start_index = J*(n - 1) + 1;
      int end_index   = J*(n - 1) + J;
      r_times_p = R[n, ] .* P[start_index:end_index, j]';
      prod_R_p_sum = dot_product(r_times_p, sum_m_w);
      //for(jx in 1:J) {
      //  prod_R_p += R[n, jx] * P[jx, j] * sum_m_w[jx];
      //}
      
      M[n,j] = prod_R_p_sum;
  
    }
  }
  
  
}

model {
   
  // ------ EQUATION (11a) -------
  // priors and sample R[n, j]
  // xsigma ~ inv_gamma(1, 1); // weak prior on st-dev
  // also tried gamma? and nothing ..
  for(j in 1:J) {
    // weak prior on beta 
    xbeta[1:(N), j] ~ normal(0, 1); 
    // sample logR
    //logR[1:(N), j] ~ normal(xbeta[1:(N), j], xsigma[j]);
    logR[1:(N), j] ~ normal(xbeta[1:(N), j], 1);
  }
  
  // ------ EQUATION (11c) -------
  // likelihood comparison
  //Y ~ poisson(M);
  for(n in 1:N) {
    for(j in 1:J) {
      target += poisson_lpmf(Y[n, j] | M[n, j]);
      //Y[n, j] ~ poisson(M[n, j])
      //;
      
      // n = 2
      // window is 1, 2, 3
      
      // calculate the mean of M[c(1, 2, 3), j]
      
      //Y[c(1, 2, 3), j] ~ poi
      
    }
  }

  
}

