data {
  int<lower=0> N;                      // Number of data points
  int<lower=0, upper=1> y[N];          // Binary outcoume: 1 = responder
  matrix[N, 4] hla;                    // Indicator variables for 4 HLA alleles
  vector<lower=0, upper=1>[N] herv;       // Indicator for whether there was any HERV response
  vector<lower=0, upper=1>[N] vir;        // Indicator for whether there was any viral antigen response
}

transformed data {
  matrix[N, 3] response;
  vector<lower=0, upper=1>[N] interaction = herv .* vir;
  response = append_col( append_col(herv, vir), interaction);
}

parameters {
  real b0;
  real mu_hla;
  real mu_response;
  real<lower=0> sigma_hla;
  real<lower=0> sigma_response;
  vector[4] b_hla_raw;
  vector[3] b_response_raw;                 // regression coeffs for: herv, vir, interaction. Hierarchical model
}

transformed parameters {
  vector[4] b_hla = b_hla_raw * sigma_hla + mu_hla;
  vector[3] b_response = b_response_raw * sigma_response + mu_response;
}

model {
  // Priors
  b0 ~ student_t(3,0,2.5);
  mu_hla ~ student_t(3,0,2.5);
  mu_response ~ student_t(3,0,2.5);
  sigma_hla ~ student_t(4,0,1);
  sigma_response ~ student_t(4,0,1);
  b_hla_raw ~ normal(0, 1);
  b_response_raw ~ normal(0, 1);

  // Likelihood
  y ~ bernoulli_logit(b0 + hla * b_hla + response * b_response);
}
