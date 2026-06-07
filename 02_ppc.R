# =============================================================================
# Bayesian Statistics — Antarctic Penguins
# 02_ppc.R : Posterior Predictive Checks
# Run after 01_modelling.R (uses fit1, fit2, fit3, penguins).
# =============================================================================

library(bayesplot)
library(ggplot2)
library(gridExtra)
library(dplyr)

theme_set(theme_bw(base_size = 13))

# =============================================================================
# 0. EXTRACT y_rep AND OBSERVED DATA
# =============================================================================
y_obs    <- penguins$bill_length
y_rep_m1 <- rstan::extract(fit1)$y_rep
y_rep_m2 <- rstan::extract(fit2)$y_rep
y_rep_m3 <- rstan::extract(fit3)$y_rep

species_group <- penguins$species
sex_group     <- penguins$sex

# 100 of the 4000 draws for readable density overlays
set.seed(42)
idx_100 <- sample(seq_len(nrow(y_rep_m1)), 100)

# =============================================================================
# 1. DENSITY OVERLAY — overall shape (the data are bimodal)
# =============================================================================
p_dens_m1 <- ppc_dens_overlay(y_obs, y_rep_m1[idx_100, ]) +
  labs(title = "M1 — Complete Pooling",
       subtitle = "Dark = real data | Light blue = 100 simulations",
       x = "bill_length (mm)") + xlim(25, 65)
p_dens_m2 <- ppc_dens_overlay(y_obs, y_rep_m2[idx_100, ]) +
  labs(title = "M2 — No Pooling",
       subtitle = "Dark = real data | Light blue = 100 simulations",
       x = "bill_length (mm)") + xlim(25, 65)
p_dens_m3 <- ppc_dens_overlay(y_obs, y_rep_m3[idx_100, ]) +
  labs(title = "M3 — Hierarchical",
       subtitle = "Dark = real data | Light blue = 100 simulations",
       x = "bill_length (mm)") + xlim(25, 65)
grid.arrange(p_dens_m1, p_dens_m2, p_dens_m3, ncol = 3,
             top = "PPC Figure 1 — Density overlay")

# =============================================================================
# 2. GLOBAL MEAN — central tendency
# =============================================================================
p_mean_m1 <- ppc_stat(y_obs, y_rep_m1, stat = "mean") + labs(title = "M1 — global mean", x = "Simulated mean (mm)")
p_mean_m2 <- ppc_stat(y_obs, y_rep_m2, stat = "mean") + labs(title = "M2 — global mean", x = "Simulated mean (mm)")
p_mean_m3 <- ppc_stat(y_obs, y_rep_m3, stat = "mean") + labs(title = "M3 — global mean", x = "Simulated mean (mm)")
grid.arrange(p_mean_m1, p_mean_m2, p_mean_m3, ncol = 3,
             top = "PPC Figure 2 — Global mean")

# =============================================================================
# 3. GROUPED MEAN BY SPECIES — the KEY diagnostic
# =============================================================================
p_grp_m1 <- ppc_stat_grouped(y_obs, y_rep_m1, group = species_group, stat = "mean") +
  labs(title = "M1 — by species", subtitle = "Vertical line = real species mean", x = "Simulated mean (mm)")
p_grp_m2 <- ppc_stat_grouped(y_obs, y_rep_m2, group = species_group, stat = "mean") +
  labs(title = "M2 — by species", subtitle = "Vertical line = real species mean", x = "Simulated mean (mm)")
p_grp_m3 <- ppc_stat_grouped(y_obs, y_rep_m3, group = species_group, stat = "mean") +
  labs(title = "M3 — by species", subtitle = "Vertical line = real species mean", x = "Simulated mean (mm)")
grid.arrange(p_grp_m1, p_grp_m2, p_grp_m3, ncol = 1,
             top = "PPC Figure 3 — Mean by species (KEY diagnostic)")

# =============================================================================
# 4. GROUPED MEAN BY SEX
# =============================================================================
p_sex_m1 <- ppc_stat_grouped(y_obs, y_rep_m1, group = sex_group, stat = "mean") + labs(title = "M1 — by sex", x = "Simulated mean (mm)")
p_sex_m2 <- ppc_stat_grouped(y_obs, y_rep_m2, group = sex_group, stat = "mean") + labs(title = "M2 — by sex", x = "Simulated mean (mm)")
p_sex_m3 <- ppc_stat_grouped(y_obs, y_rep_m3, group = sex_group, stat = "mean") + labs(title = "M3 — by sex", x = "Simulated mean (mm)")
grid.arrange(p_sex_m1, p_sex_m2, p_sex_m3, ncol = 3,
             top = "PPC Figure 4 — Mean by sex")

# =============================================================================
# 5. STANDARD DEVIATION (global, then by species)
# =============================================================================
p_sd_m1 <- ppc_stat(y_obs, y_rep_m1, stat = "sd") + labs(title = "M1 — SD", x = "Simulated SD (mm)")
p_sd_m2 <- ppc_stat(y_obs, y_rep_m2, stat = "sd") + labs(title = "M2 — SD", x = "Simulated SD (mm)")
p_sd_m3 <- ppc_stat(y_obs, y_rep_m3, stat = "sd") + labs(title = "M3 — SD", x = "Simulated SD (mm)")
grid.arrange(p_sd_m1, p_sd_m2, p_sd_m3, ncol = 3,
             top = "PPC Figure 5 — Standard deviation")

p_sd_sp_m1 <- ppc_stat_grouped(y_obs, y_rep_m1, group = species_group, stat = "sd") +
  labs(title = "M1 — SD by species", subtitle = "Vertical line = observed species SD", x = "Simulated SD (mm)")
p_sd_sp_m2 <- ppc_stat_grouped(y_obs, y_rep_m2, group = species_group, stat = "sd") +
  labs(title = "M2 — SD by species", subtitle = "Vertical line = observed species SD", x = "Simulated SD (mm)")
p_sd_sp_m3 <- ppc_stat_grouped(y_obs, y_rep_m3, group = species_group, stat = "sd") +
  labs(title = "M3 — SD by species", subtitle = "Vertical line = observed species SD", x = "Simulated SD (mm)")
grid.arrange(p_sd_sp_m1, p_sd_sp_m2, p_sd_sp_m3, ncol = 1,
             top = "PPC Figure — SD by species")

# =============================================================================
# 6. MINIMUM vs MAXIMUM
# =============================================================================
p_2d_m1 <- ppc_stat_2d(y_obs, y_rep_m1, stat = c("min", "max")) +
  ggtitle("M1 — Min vs Max") + xlab("min (mm)") + ylab("max (mm)") + theme(legend.position = "bottom")
p_2d_m2 <- ppc_stat_2d(y_obs, y_rep_m2, stat = c("min", "max")) +
  ggtitle("M2 — Min vs Max") + xlab("min (mm)") + ylab("max (mm)") + theme(legend.position = "bottom")
p_2d_m3 <- ppc_stat_2d(y_obs, y_rep_m3, stat = c("min", "max")) +
  ggtitle("M3 — Min vs Max") + xlab("min (mm)") + ylab("max (mm)") + theme(legend.position = "bottom")
grid.arrange(p_2d_m1, p_2d_m2, p_2d_m3, ncol = 3,
             top = "PPC Figure 6 — Min vs Max (joint)")

# =============================================================================
# 7. BAYESIAN P-VALUES
# =============================================================================
ppc_pvalue <- function(y, y_rep, stat_fn) {
  round(mean(apply(y_rep, 1, stat_fn) > stat_fn(y)), 3)
}

cat("Global statistics (p-value):\n")
for (st in c("mean", "sd", "min", "max")) {
  f <- match.fun(st)
  cat(sprintf("  %-4s  M1: %.3f | M2: %.3f | M3: %.3f\n", st,
              ppc_pvalue(y_obs, y_rep_m1, f),
              ppc_pvalue(y_obs, y_rep_m2, f),
              ppc_pvalue(y_obs, y_rep_m3, f)))
}

cat("\nMean by species (p-value — most discriminating):\n")
for (sp in c("Adelie", "Chinstrap", "Gentoo")) {
  idx <- which(penguins$species == sp)
  cat(sprintf("  %-9s  M1: %.3f | M2: %.3f | M3: %.3f\n", sp,
              ppc_pvalue(y_obs[idx], y_rep_m1[, idx], mean),
              ppc_pvalue(y_obs[idx], y_rep_m2[, idx], mean),
              ppc_pvalue(y_obs[idx], y_rep_m3[, idx], mean)))
}

