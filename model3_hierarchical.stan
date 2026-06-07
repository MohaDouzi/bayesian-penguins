# MODÈLE 3 : Hiérarchique (Partial Pooling)


data {
  int<lower=0> N;
  int<lower=1> K;              
  vector[N] bill_length;
  vector[N] bill_depth;        
  vector[N] sex;               
  array[N] int<lower=1, upper=K> species;
}

parameters {
  real mu_alpha;              
  real<lower=0> sigma_alpha;   

  vector[K] alpha;             
  real beta_depth;            
  real beta_sex;               
  real<lower=0> sigma;        
}

model {
  mu_alpha    ~ normal(45, 10);    
  sigma_alpha ~ normal(0, 5);

  alpha ~ normal(mu_alpha, sigma_alpha);

  beta_depth ~ normal(0, 5);
  beta_sex   ~ normal(0, 3);
  sigma       ~ normal(0, 5);

  for (n in 1:N) {
    bill_length[n] ~ normal(
      alpha[species[n]]
      + beta_depth * bill_depth[n]
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
                + beta_depth * bill_depth[n]
                + beta_sex * sex[n];
    log_lik[n] = normal_lpdf(bill_length[n] | mu_n, sigma);
    y_rep[n]   = normal_rng(mu_n, sigma);
  }
}
