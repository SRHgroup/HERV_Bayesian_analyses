data {
  int<lower=0> N;                      // Number of data points
  int<lower=0, upper=1> y[N];          // Binary outcoume: 1 = responder
  matrix[N, 4] hla;                    // Indicator variables for 4 HLA alleles
  matrix[N, 4] pepcount;               // Number of peptides recognized (standardized to z-score)
}

parameters {
  real b0;
  real mu_hla;
  real mu_pep;
  real<lower=0> sigma_hla;
  real<lower=0> sigma_pep;
  vector[4] b_hla_raw;
  vector[4] b_pep_raw;
}

transformed parameters {
  // Non-centered parameterization
  vector[4] b_hla = b_hla_raw * sigma_hla + mu_hla;
  vector[4] b_pep = b_pep_raw * sigma_pep + mu_pep;
}

model {
  // Priors
  b0 ~ student_t(3,0,2.5);
  mu_hla ~ student_t(3,0,2.5);
  mu_pep ~ student_t(3,0,2.5);
  sigma_hla ~ student_t(4,0,1);
  sigma_pep ~ student_t(4,0,1);
  b_hla_raw ~ normal(0, 1);
  b_pep_raw ~ normal(0, 1);

  // Likelihood
  y ~ bernoulli_logit(b0 + hla * b_hla + pepcount * b_pep);
}
