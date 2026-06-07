// =============================================================================
// MODEL 2 — No Pooling / Separate (non-hierarchical, sensitivity version)
// =============================================================================
data {
  int<lower=0> N;
  int<lower=1> K;                             
  vector[N]    bill_length;
  vector[N]    bill_depth;                     
  vector[N]    sex;                           
  array[N] int<lower=1, upper=K> species;

  real          prior_alpha_mu;
  real<lower=0> prior_alpha_sd;
  real          prior_beta_depth_mu;
  real<lower=0> prior_beta_depth_sd;
  real          prior_beta_sex_mu;
  real<lower=0> prior_beta_sex_sd;
  real<lower=0> prior_sigma_sd;                
}
parameters {
  vector[K]     alpha;        
  vector[K]     beta_depth;   
  real          beta_sex;     
  real<lower=0> sigma;        
}
model {
  alpha      ~ normal(prior_alpha_mu,      prior_alpha_sd);
  beta_depth ~ normal(prior_beta_depth_mu, prior_beta_depth_sd);
  beta_sex   ~ normal(prior_beta_sex_mu,   prior_beta_sex_sd);
  sigma      ~ normal(0, prior_sigma_sd);
  
  for (n in 1:N) {
    bill_length[n] ~ normal(
      alpha[species[n]]
      + beta_depth[species[n]] * bill_depth[n]
      + beta_sex * sex[n],
      sigma
    );
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    real mu_n = alpha[species[n]]
                + beta_depth[species[n]] * bill_depth[n]
                + beta_sex * sex[n];
    log_lik[n] = normal_lpdf(bill_length[n] | mu_n, sigma);
    y_rep[n]   = normal_rng(mu_n, sigma);
  }
}
