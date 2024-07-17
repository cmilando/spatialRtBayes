data {
  int<lower=1> N; // days of observed data

  int casesa[N]; // reported cases
  int casesb[N]; // reported cases
  int casesc[N]; // reported cases

  real SI[N]; // fixed SI using empirical data
  
}

transformed data {
  vector[N] SI_rev; // SI in reverse order
  for(i in 1:N)
  SI_rev[i] = SI[N-i+1];
  
}

parameters {

  vector[N] weekly_effecta;
  vector[N] weekly_effectb;
  vector[N] weekly_effectc;
  
  vector<lower=0, upper=1>[N] weekly_rhoa;
  vector<lower=0, upper=1>[N] weekly_rhob;
  vector<lower=0, upper=1>[N] weekly_rhoc;

  real<lower=0> weekly_sda;
  real<lower=0> weekly_sdb;
  real<lower=0> weekly_sdc;

}

transformed parameters {

  vector[N] predictiona=rep_vector(1e-5,N);
  vector[N] predictionb=rep_vector(1e-5,N);
  vector[N] predictionc=rep_vector(1e-5,N);
  
  vector<lower=0>[N] Rta;
  vector<lower=0>[N] Rtb;
  vector<lower=0>[N] Rtc;
  
  {
  
    Rta[1:N] = exp(weekly_effecta[1:N]);
    Rtb[1:N] = exp(weekly_effectb[1:N]);
    Rtc[1:N] = exp(weekly_effectc[1:N]);

    for (i in 2:N) {
    
      real convolutiona = dot_product(predictiona[1:(i-1)], tail(SI_rev, i-1));
      real convolutionb = dot_product(predictionb[1:(i-1)], tail(SI_rev, i-1));
      real convolutionc = dot_product(predictionc[1:(i-1)], tail(SI_rev, i-1));
      
      predictiona[i] = predictiona[i] + Rta[i] * convolutiona * 0.8 + Rtb[i] * convolutionb * 0.2 + Rtc[i] * convolutionc * 0.1;
      predictionb[i] = predictionb[i] + Rta[i] * convolutiona * 0.15 + Rtb[i] * convolutionb * 0.6 + Rtc[i] * convolutionc * 0.3;
      predictionc[i] = predictionc[i] + Rta[i] * convolutiona * 0.05 + Rtb[i] * convolutionb * 0.2 + Rtc[i] * convolutionc * 0.6;
      
    }
    
  }
}

model {

  weekly_sda ~ normal(0,1);
  weekly_sdb ~ normal(0,1);
  weekly_sdc ~ normal(0,1);
  
  weekly_rhoa ~ normal(0.8, 0.5);
  weekly_rhob ~ normal(0.8, 0.5);
  weekly_rhoc ~ normal(0.8, 0.5);
  
  weekly_effecta[1:(N)] ~ normal(weekly_rhoa[1:(N)] , weekly_sda);
  weekly_effectb[1:(N)] ~ normal(weekly_rhob[1:(N)] , weekly_sdb);
  weekly_effectc[1:(N)] ~ normal(weekly_rhoc[1:(N)] , weekly_sdc);
  
  casesa ~ poisson(predictiona);
  casesb ~ poisson(predictionb);
  casesc ~ poisson(predictionc);
}