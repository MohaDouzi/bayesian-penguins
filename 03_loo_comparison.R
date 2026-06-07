# =============================================================================
# Bayesian Statistics — Antarctic Penguins
# 03_loo_comparison.R : Model comparison via PSIS-LOO
# =============================================================================

library(loo)
library(rstan)

# =============================================================================
# 1. COMPUTE LOO FOR EACH MODEL
# =============================================================================
log_lik_m1 <- extract_log_lik(fit1, merge_chains = FALSE)
log_lik_m2 <- extract_log_lik(fit2, merge_chains = FALSE)
log_lik_m3 <- extract_log_lik(fit3, merge_chains = FALSE)

loo_m1 <- loo(log_lik_m1, r_eff = relative_eff(exp(log_lik_m1)))
loo_m2 <- loo(log_lik_m2, r_eff = relative_eff(exp(log_lik_m2)))
loo_m3 <- loo(log_lik_m3, r_eff = relative_eff(exp(log_lik_m3)))

print(loo_m1); print(loo_m2); print(loo_m3)

# =============================================================================
# 2. PARETO-k DIAGNOSTICS
# =============================================================================
par(mfrow = c(1, 3))
plot(loo_m1, main = "Pareto-k — M1")
plot(loo_m2, main = "Pareto-k — M2")
plot(loo_m3, main = "Pareto-k — M3")
par(mfrow = c(1, 1))

high_k_m2 <- which(loo_m2$diagnostics$pareto_k > 0.7)
high_k_m3 <- which(loo_m3$diagnostics$pareto_k > 0.7)
if (length(high_k_m2) > 0) print(penguins[high_k_m2, ])
if (length(high_k_m3) > 0) print(penguins[high_k_m3, ])

# =============================================================================
# 3. MODEL COMPARISON
# =============================================================================
print(loo_compare(loo_m1, loo_m2, loo_m3))

# =============================================================================
# 4. PRACTICAL PREDICTIVE ACCURACY (approximate LOO-RMSE, in mm)
# =============================================================================
loo_rmse <- function(loo_obj, sigma) {
  m <- mean(loo_obj$pointwise[, "elpd_loo"])
  inside <- -2 * m - log(2 * pi * sigma^2)
  sigma * sqrt(max(inside, 0))                
}
sigma_m2 <- mean(rstan::extract(fit2)$sigma)
sigma_m3 <- mean(rstan::extract(fit3)$sigma)

cat(sprintf("Approx. LOO-RMSE  M2: %.2f mm | M3: %.2f mm\n",
            loo_rmse(loo_m2, sigma_m2), loo_rmse(loo_m3, sigma_m3)))

# In-sample RMSE (optimistic) — if RMSE ~ sigma the model is not overfitting.
mu_hat_m2 <- colMeans(rstan::extract(fit2)$y_rep)
mu_hat_m3 <- colMeans(rstan::extract(fit3)$y_rep)
cat(sprintf("In-sample RMSE    M2: %.2f mm (sigma %.2f) | M3: %.2f mm (sigma %.2f)\n",
            sqrt(mean((penguins$bill_length - mu_hat_m2)^2)), sigma_m2,
            sqrt(mean((penguins$bill_length - mu_hat_m3)^2)), sigma_m3))

