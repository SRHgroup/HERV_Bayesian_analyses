data {
  int<lower=1> N;                 // Number of individuals (controls + patients)
  int<lower=1> N_group;           // Number of groupings 
  int<lower=1> N_hla;             // Number of HLA alleles
  int<lower=1, upper=N_group> group[N]; // group ID: 1: control, 2: patient, 3: patient+aza)
  int<lower=1, upper=N_hla> hla[N];     // Allele: 1=A0101, 2=A0201, 3=B0702, 4=B0801
  int<lower=0> n[N];              // Binomial n: number of tested peptides
  int<lower=0> y[N];              // Binomial utcome: number of positive peptides
}

parameters {
  real b0;                           // Intercept: Background proportion of positive peptides
  vector[N_group - 1] b_group_raw;   // Slope for effect of being a patient
  vector[N_hla - 1] b_hla_raw;       // Slopes for effects of alleles
}

transformed parameters {
  // Work around for sum to zero constraints on deflection parameters
  vector[N_group] b_group = append_row( b_group_raw, -sum( b_group_raw ));
  vector[N_hla] b_hla = append_row( b_hla_raw, -sum( b_hla_raw ));
}

model {
  // Priors
  b0 ~ normal(-4, 3);
  b_group ~ normal(0, inv(sqrt(1 - inv(N_group)))); // Ensures standard normal prior
  b_hla ~ normal(0, inv(sqrt(1 - inv(N_hla))));

  // Likelihood
  for (i in 1:N) {
      y[i] ~ binomial_logit(n[i], b0 + b_group[group[i]] + b_hla[hla[i]]);
  }
}

generated quantities {
  real<lower=0, upper=1> p_control = inv_logit(b0 + b_group[1]);
  real<lower=0, upper=1> p_patient = inv_logit(b0 + b_group[2]);
  real<lower=0, upper=1> p_aza = inv_logit(b0 + b_group[3]);
  real d_pt_con = p_patient - p_control;
  real d_aza_con = p_aza - p_control;
  real d_aza_pt = p_aza - p_patient;
  real b_pt_con = b_group[2] - b_group[1];
  real b_aza_con = b_group[3] - b_group[1];
  real b_aza_pt = b_group[3] - b_group[2];
  real rr_pt_con = p_patient / p_control;
  real rr_aza_con = p_aza / p_control;
  real rr_aza_pt = p_aza / p_patient;
}

