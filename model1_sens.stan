// =============================================================================
// MODEL 1 — Complete Pooling (sensitivity version)
// =============================================================================
data {
  int<lower=0> N;
  vector[N] bill_length;
  vector[N] bill_depth_std;

  real          prior_alpha_mu;     
  real<lower=0> prior_alpha_sd;     
  real          prior_beta_mu;      
  real<lower=0> prior_beta_sd;      
  real<lower=0> prior_sigma_sd;     
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(prior_alpha_mu, prior_alpha_sd);
  beta  ~ normal(prior_beta_mu,  prior_beta_sd);
  sigma ~ normal(0, prior_sigma_sd);

  bill_length ~ normal(alpha + beta * bill_depth_std, sigma);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  for (n in 1:N) {
    real mu_n  = alpha + beta * bill_depth_std[n];
    log_lik[n] = normal_lpdf(bill_length[n] | mu_n, sigma);
    y_rep[n]   = normal_rng(mu_n, sigma);
  }
}
