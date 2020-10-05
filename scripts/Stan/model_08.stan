data {
  int<lower=1> N;                 // Number of measurements
  int<lower=1> N_group;           // Number of groupings
  int<lower=1, upper=N_group> group[N]; // group ID: 1: control, 2: patient, 3: patient+aza)
  int<lower=0> n[N];              // Binomial n: number peptides tested
  int<lower=0> y[N];              // Binomial utcome: number peptides positive
}

parameters {
  vector<lower=0,upper=1>[N_group] theta;    // Proportion positive for each group
}

model{
  // prior
  theta ~ beta(1,1);
  
  // likelihood
  for (i in 1:N) {
    y[i] ~ binomial(n[i], theta[group[i]]);
  }
}

generated quantities {
  real d_pt_con = theta[2] - theta[1];
  real d_aza_con = theta[3] - theta[1];
  real d_aza_pt = theta[3] - theta[2];
  real rr_pt_con = theta[2] / theta[1];
  real rr_aza_con = theta[3] / theta[1];
  real rr_aza_pt = theta[3] / theta[2];
}
