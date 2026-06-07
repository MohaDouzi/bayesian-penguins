# =============================================================================
# Bayesian Statistics — Antarctic Penguins
# 01_modelling.R : Complete Pooling (M1), No Pooling (M2), Hierarchical (M3)
# =============================================================================

library(rstan)
library(ggplot2)
library(dplyr)
library(bayesplot)
library(tidyr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# =============================================================================
# DATA PREPARATION
# =============================================================================

penguins <- readRDS("penguins.RDS")

# Summary statistics for interpretation only. Priors are NOT derived from them
depth_mean  <- mean(penguins$bill_depth)
depth_sd    <- sd(penguins$bill_depth)

# Standardise bill_depth (mean = 0, sd = 1):
penguins$bill_depth_std <- (penguins$bill_depth - depth_mean) / depth_sd

# Species index: Adelie = 1, Chinstrap = 2, Gentoo = 3
penguins$species_num <- as.integer(factor(penguins$species,
                                          levels = c("Adelie", "Chinstrap", "Gentoo")))

# Sex: female = 0, male = 1
penguins$sex_num <- ifelse(penguins$sex == "male", 1, 0)

species_summary <- penguins %>%
  group_by(species) %>%
  summarise(n = n(), mean = round(mean(bill_length), 2),
            sd = round(sd(bill_length), 2), .groups = "drop")
print(species_summary)

species_colors <- c("Adelie" = "darkorange", "Chinstrap" = "purple", "Gentoo" = "steelblue")
species_labels <- c("Adelie", "Chinstrap", "Gentoo")

set.seed(42)
N_sim <- 2000

# =============================================================================
# MODEL 1 — COMPLETE POOLING
#   bill_length_i ~ N(alpha + beta * depth_std_i, sigma)
# =============================================================================

# --- Prior predictive check ---
alpha_ppc <- rnorm(N_sim, 45, 10)
beta_ppc  <- rnorm(N_sim,  0,  5)
sigma_ppc <- abs(rnorm(N_sim, 0, 5))            # Half-Normal(0, 5)
y_ppc_m1  <- rnorm(N_sim, alpha_ppc + beta_ppc * 0, sigma_ppc)

round(quantile(y_ppc_m1, c(0.025, 0.975)), 1)  
round(mean(y_ppc_m1 < 10 | y_ppc_m1 > 100) * 100, 1)

# --- Fit ---
stan_data_m1 <- list(
  N              = nrow(penguins),
  bill_length    = penguins$bill_length,
  bill_depth_std = penguins$bill_depth_std
)
fit1 <- stan(file = "model1_complete_pooling.stan", data = stan_data_m1,
             iter = 2000, warmup = 1000, chains = 4, seed = 42)

# --- Diagnostics (target: Rhat < 1.01, n_eff > 400, 0 divergences) ---
print(fit1, pars = c("alpha", "beta", "sigma"), digits_summary = 3)
check_hmc_diagnostics(fit1)

# --- Posterior ---
s1 <- rstan::extract(fit1)
mean(s1$alpha);  quantile(s1$alpha, c(0.025, 0.975))   # intercept at average depth
mean(s1$beta);   quantile(s1$beta,  c(0.025, 0.975))   # NEGATIVE: Simpson's paradox artefact
mean(s1$beta) / depth_sd                               # slope in mm per mm of depth
mean(s1$sigma);  quantile(s1$sigma, c(0.025, 0.975))   # large (~global SD): explains little

# =============================================================================
# MODEL 2 — NO POOLING (separate per species)
#   bill_length_i ~ N(alpha[k] + beta_depth[k]*depth_std_i + beta_sex*sex_i, sigma)
# =============================================================================

# --- Prior predictive check ---
alpha_ppc2   <- rnorm(N_sim, 45, 10)
beta_k_ppc   <- rnorm(N_sim,  0,  5)
beta_sex_ppc <- rnorm(N_sim,  0,  3)
sigma_ppc2   <- abs(rnorm(N_sim, 0, 5))
y_ppc_female <- rnorm(N_sim, alpha_ppc2 + beta_sex_ppc * 0, sigma_ppc2)
y_ppc_male   <- rnorm(N_sim, alpha_ppc2 + beta_sex_ppc * 1, sigma_ppc2)
round(mean(y_ppc_female < 10 | y_ppc_female > 100) * 100, 1)
round(mean(y_ppc_male   < 10 | y_ppc_male   > 100) * 100, 1)

# --- Fit ---
stan_data_m2 <- list(
  N           = nrow(penguins), K = 3,
  bill_length = penguins$bill_length,
  bill_depth  = penguins$bill_depth_std,
  sex         = penguins$sex_num,
  species     = penguins$species_num
)
fit2 <- stan(file = "model2_no_pooling.stan", data = stan_data_m2,
             iter = 2000, warmup = 1000, chains = 4, seed = 42)

# --- Diagnostics ---
print(fit2, pars = c("alpha", "beta_depth", "beta_sex", "sigma"), digits_summary = 3)
check_hmc_diagnostics(fit2)
diag_m2 <- summary(fit2)$summary
diag_m2 <- diag_m2[!is.nan(diag_m2[, "Rhat"]), ]
max(diag_m2[, "Rhat"]); sum(diag_m2[, "Rhat"] > 1.01); min(diag_m2[, "n_eff"])

# --- Posterior ---
s2 <- rstan::extract(fit2)

# Intercepts = predicted bill_length for a FEMALE of species k at average depth.
sapply(1:3, function(k) round(mean(s2$alpha[, k]), 2))  
sapply(1:3, function(k) round(quantile(s2$alpha[, k], c(.025,.975)), 2))

sapply(1:3, function(k) round(mean(s2$beta_depth[, k]), 2))
sapply(1:3, function(k) round(mean(s2$beta_depth[, k]) / depth_sd, 2))
sapply(1:3, function(k) round(quantile(s2$beta_depth[, k], c(.025,.975)), 2))  # Adelie CI overlaps 0

mean(s2$beta_sex);  quantile(s2$beta_sex, c(0.025, 0.975))   
mean(s2$sigma);     quantile(s2$sigma,    c(0.025, 0.975)) 
# --- Plots ---
df_alpha_m2 <- data.frame(value = as.vector(s2$alpha),
                          species = rep(species_labels, each = nrow(s2$alpha)))
ggplot(df_alpha_m2, aes(value, fill = species)) +
  geom_density(alpha = 0.7, color = "white") +
  geom_vline(data = species_summary, aes(xintercept = mean, color = species),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) + guides(color = "none") +
  labs(title = "M2 — Posterior distributions of species intercepts",
       subtitle = "Dashed line = observed species mean",
       x = "Intercept alpha[k] (mm)", y = "Density", fill = "Species") +
  theme_bw(base_size = 13)

df_beta_m2 <- data.frame(value = as.vector(s2$beta_depth),
                         species = rep(species_labels, each = nrow(s2$beta_depth)))
ggplot(df_beta_m2, aes(value, fill = species)) +
  geom_density(alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = species_colors) +
  labs(title = "M2 — Posterior distributions of bill_depth slopes per species",
       subtitle = "All slopes positive within species (vs M1's negative slope)",
       x = "beta_depth[k] (mm/SD)", y = "Density", fill = "Species") +
  theme_bw(base_size = 13)

x_grid    <- seq(min(penguins$bill_depth_std), max(penguins$bill_depth_std), length.out = 100)
x_grid_mm <- x_grid * depth_sd + depth_mean
lines_m2 <- do.call(rbind, lapply(1:3, function(k) {
  a <- mean(s2$alpha[, k]); b <- mean(s2$beta_depth[, k]); bsx <- mean(s2$beta_sex)
  rbind(data.frame(x = x_grid_mm, y = a + b * x_grid,       species = species_labels[k], sex = "female"),
        data.frame(x = x_grid_mm, y = a + b * x_grid + bsx, species = species_labels[k], sex = "male"))
}))
ggplot() +
  geom_point(data = penguins, aes(bill_depth, bill_length, color = species, shape = sex),
             alpha = 0.5, size = 2) +
  geom_line(data = lines_m2, aes(x, y, color = species, linetype = sex), linewidth = 1) +
  scale_color_manual(values = species_colors) +
  scale_linetype_manual(values = c(female = "dashed", male = "solid")) +
  labs(title = "M2 — Posterior regression lines (mean estimate)",
       subtitle = "One line per species x sex combination",
       x = "bill_depth (mm)", y = "bill_length (mm)", color = "Species", linetype = "Sex") +
  theme_bw(base_size = 13)

# =============================================================================
# MODEL 3 — HIERARCHICAL (partial pooling)
#   bill_length_i ~ N(alpha[k] + beta_depth*depth_std_i + beta_sex*sex_i, sigma)
#   alpha[k] ~ N(mu_alpha, sigma_alpha)
# =============================================================================

# --- Prior predictive check  ---
mu_alpha_ppc    <- rnorm(N_sim, 45, 10)
sigma_alpha_ppc <- abs(rnorm(N_sim, 0, 5))
alpha_ppc3      <- rnorm(N_sim, mu_alpha_ppc, sigma_alpha_ppc)
beta_depth_ppc  <- rnorm(N_sim, 0, 5)
beta_sex_ppc3   <- rnorm(N_sim, 0, 3)
sigma_ppc3      <- abs(rnorm(N_sim, 0, 5))
y_ppc3_female   <- rnorm(N_sim, alpha_ppc3 + beta_sex_ppc3 * 0, sigma_ppc3)
y_ppc3_male     <- rnorm(N_sim, alpha_ppc3 + beta_sex_ppc3 * 1, sigma_ppc3)
round(mean(y_ppc3_female < 10 | y_ppc3_female > 100) * 100, 1)
round(mean(y_ppc3_male   < 10 | y_ppc3_male   > 100) * 100, 1)
round(mean(sigma_alpha_ppc < 10) * 100, 1)      # >95% of mass below 10 mm

# --- Fit ---
stan_data_m3 <- stan_data_m2                     # same data structure as M2
fit3 <- stan(file = "model3_hierarchical.stan", data = stan_data_m3,
             iter = 2000, warmup = 1000, chains = 4, seed = 42)

# --- Diagnostics ---
print(fit3, pars = c("mu_alpha", "sigma_alpha"), digits_summary = 3)
print(fit3, pars = c("alpha", "beta_depth", "beta_sex", "sigma"), digits_summary = 3)
check_hmc_diagnostics(fit3)
diag_m3 <- summary(fit3)$summary
diag_m3 <- diag_m3[!is.nan(diag_m3[, "Rhat"]), ]
max(diag_m3[, "Rhat"]); sum(diag_m3[, "Rhat"] > 1.01); min(diag_m3[, "n_eff"])

# --- Posterior ---
s3 <- rstan::extract(fit3)
mean(s3$mu_alpha);    quantile(s3$mu_alpha, c(0.025, 0.975))      
mean(s3$sigma_alpha); quantile(s3$sigma_alpha, c(0.025, 0.975))   
sapply(1:3, function(k) round(mean(s3$alpha[, k]), 2))            
mean(s3$beta_depth);  mean(s3$beta_depth) / depth_sd; quantile(s3$beta_depth, c(0.025, 0.975))
mean(s3$beta_sex);    quantile(s3$beta_sex, c(0.025, 0.975))
mean(s3$sigma);       quantile(s3$sigma, c(0.025, 0.975))         

# --- Plots ---
df_alpha_m3 <- data.frame(value = as.vector(s3$alpha),
                          species = rep(species_labels, each = nrow(s3$alpha)))
ggplot(df_alpha_m3, aes(value, fill = species)) +
  geom_density(alpha = 0.7, color = "white") +
  geom_vline(xintercept = mean(s3$mu_alpha), linetype = "dashed",
             color = "black", linewidth = 0.8) +
  annotate("text", x = mean(s3$mu_alpha) + 0.4, y = 0.02, label = "mu_alpha",
           hjust = 0, size = 3.5) +
  scale_fill_manual(values = species_colors) +
  labs(title = "M3 — Posterior distributions of species intercepts",
       subtitle = "Dashed line = posterior mean of mu_alpha (shrinkage target)",
       x = "Intercept alpha[k] (mm)", y = "Density", fill = "Species") +
  theme_bw(base_size = 13)

ggplot(data.frame(value = s3$sigma_alpha), aes(value)) +
  geom_density(fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = mean(s3$sigma_alpha), linetype = "dashed",
             color = "darkblue", linewidth = 0.8) +
  labs(title = "M3 — Posterior of sigma_alpha (between-species variability)",
       subtitle = paste0("Posterior mean = ", round(mean(s3$sigma_alpha), 2),
                          " mm — above 0, hierarchical model is justified"),
       x = "sigma_alpha (mm)", y = "Density") +
  theme_bw(base_size = 13)

lines_m3 <- do.call(rbind, lapply(1:3, function(k) {
  a <- mean(s3$alpha[, k]); b <- mean(s3$beta_depth); bsx <- mean(s3$beta_sex)
  rbind(data.frame(x = x_grid_mm, y = a + b * x_grid,       species = species_labels[k], sex = "female"),
        data.frame(x = x_grid_mm, y = a + b * x_grid + bsx, species = species_labels[k], sex = "male"))
}))
ggplot() +
  geom_point(data = penguins, aes(bill_depth, bill_length, color = species, shape = sex),
             alpha = 0.5, size = 2) +
  geom_line(data = lines_m3, aes(x, y, color = species, linetype = sex), linewidth = 1) +
  scale_color_manual(values = species_colors) +
  scale_linetype_manual(values = c(female = "dashed", male = "solid")) +
  labs(title = "M3 — Posterior regression lines (mean estimate)",
       subtitle = "Shared slope across species; separate intercepts",
       x = "bill_depth (mm)", y = "bill_length (mm)", color = "Species", linetype = "Sex") +
  theme_bw(base_size = 13)

# =============================================================================
# SHRINKAGE: M2 vs M3 INTERCEPTS
# =============================================================================
m2_intercepts <- colMeans(s2$alpha)
m3_intercepts <- colMeans(s3$alpha)
global_mean   <- mean(s3$mu_alpha)
shrinkage_df <- data.frame(
  species   = species_labels,
  n_obs     = c(146, 68, 119),
  M2        = round(m2_intercepts, 2),
  M3        = round(m3_intercepts, 2),
  shrinkage = round(abs(m2_intercepts - m3_intercepts), 2)
)
print(shrinkage_df)   # Chinstrap (n = 68) is pulled most toward the global mean

ggplot(shrinkage_df) +
  geom_segment(aes(x = M2, xend = M3, y = species, yend = species, color = species),
               arrow = arrow(length = unit(0.25, "cm")), linewidth = 1.2) +
  geom_vline(xintercept = global_mean, linetype = "dashed", color = "gray40") +
  annotate("text", x = global_mean + 0.3, y = 0.6, label = "global mean",
           size = 3.5, color = "gray40") +
  scale_color_manual(values = species_colors) +
  labs(title = "Shrinkage effect: M2 (no pooling) -> M3 (partial pooling)",
       subtitle = "Arrows point from M2 estimate to M3 estimate",
       x = "Intercept (mm)", y = NULL, color = "Species") +
  theme_bw(base_size = 13)

# =============================================================================
# CROSS-MODEL SIGMA COMPARISON
# =============================================================================
sigma_comparison <- data.frame(
  model      = c("M1 Complete Pooling", "M2 No Pooling", "M3 Hierarchical"),
  sigma_mean = c(mean(s1$sigma), mean(s2$sigma), mean(s3$sigma)),
  lower95    = c(quantile(s1$sigma, .025), quantile(s2$sigma, .025), quantile(s3$sigma, .025)),
  upper95    = c(quantile(s1$sigma, .975), quantile(s2$sigma, .975), quantile(s3$sigma, .975))
)
ggplot(sigma_comparison, aes(model, sigma_mean, color = model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2, linewidth = 1) +
  scale_color_manual(values = c("red", "steelblue", "darkgreen")) +
  labs(title = "Residual noise (sigma) across models",
       subtitle = "M2 and M3 absorb between-species variance; M1 does not",
       x = NULL, y = "sigma (mm)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")

