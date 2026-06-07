// MODÈLE 1 : Complete Pooling (non-hiérarchique)

data {
  int<lower=0> N;
  vector[N] bill_length;
  vector[N] bill_depth_std;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(45, 10);
  beta  ~ normal(0, 5);
  sigma ~ normal(0, 5);
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
