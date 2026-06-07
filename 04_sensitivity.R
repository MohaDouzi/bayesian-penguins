# =============================================================================
# Sensitivity analysis — Normal models (report Section 6)
# 3 prior regimes x 3 models -> PPC plots one by one
# =============================================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ── Data ──────────────────────────────────────────────────────────────────────
penguins <- readRDS("penguins.RDS")
penguins$bill_depth_std <- (penguins$bill_depth - mean(penguins$bill_depth)) / sd(penguins$bill_depth)
penguins$species_num    <- as.integer(factor(penguins$species,
                                             levels = c("Adelie", "Chinstrap", "Gentoo")))
penguins$sex_num <- ifelse(penguins$sex == "male", 1, 0)
y_obs <- penguins$bill_length
N <- nrow(penguins)
K <- 3

# ── Prior regimes ─────────────────────────────────────────────────────────────
# base  = main models (reference)
# large = sd x5, quasi non-informative
# serre = sd /3, more informative
priors <- list(
  base = list(
    p_alpha_mu = 45, p_alpha_sd = 10,
    p_beta_mu = 0,   p_beta_sd = 5,
    p_beta_sex_mu = 0, p_beta_sex_sd = 3,
    p_sigma_sd = 5,
    p_mu_alpha_mu = 45, p_mu_alpha_sd = 10,
    p_sigma_alpha_sd = 5
  ),
  diffuse = list(
    p_alpha_mu = 45, p_alpha_sd = 50,
    p_beta_mu = 0,   p_beta_sd = 25,
    p_beta_sex_mu = 0, p_beta_sex_sd = 15,
    p_sigma_sd = 25,
    p_mu_alpha_mu = 45, p_mu_alpha_sd = 50,
    p_sigma_alpha_sd = 25
  ),
  informative = list(
    p_alpha_mu = 45, p_alpha_sd = 3,
    p_beta_mu = 0,   p_beta_sd = 2,
    p_beta_sex_mu = 0, p_beta_sex_sd = 1,
    p_sigma_sd = 1.67,
    p_mu_alpha_mu = 45, p_mu_alpha_sd = 3,
    p_sigma_alpha_sd = 1.67
  )
)

# ── Stan models ───────────────────────────────────────────────────────────────
sm1 <- stan_model("model1_sens.stan")
sm2 <- stan_model("model2_sens.stan")
sm3 <- stan_model("model3_sens.stan")

# ── Fitting helper ────────────────────────────────────────────────────────────
run_fit <- function(sm, data_list) {
  sampling(sm, data = data_list,
           iter = 2000, warmup = 1000, chains = 4,
           seed = 42, refresh = 0,
           control = list(adapt_delta = 0.95))
}

# ── Fit all 9 combinations ────────────────────────────────────────────────────
fits <- list(M1 = list(), M2 = list(), M3 = list())

for (pn in names(priors)) {
  p <- priors[[pn]]
  
  fits$M1[[pn]] <- run_fit(sm1, list(
    N = N, bill_length = y_obs,
    bill_depth_std   = penguins$bill_depth_std,
    prior_alpha_mu   = p$p_alpha_mu,
    prior_alpha_sd   = p$p_alpha_sd,
    prior_beta_mu    = p$p_beta_mu,
    prior_beta_sd    = p$p_beta_sd,
    prior_sigma_sd   = p$p_sigma_sd
  ))
  
  fits$M2[[pn]] <- run_fit(sm2, list(
    N = N, K = K,
    bill_length         = y_obs,
    bill_depth          = penguins$bill_depth_std,
    sex                 = penguins$sex_num,
    species             = penguins$species_num,
    prior_alpha_mu      = p$p_alpha_mu,
    prior_alpha_sd      = p$p_alpha_sd,
    prior_beta_depth_mu = p$p_beta_mu,
    prior_beta_depth_sd = p$p_beta_sd,
    prior_beta_sex_mu   = p$p_beta_sex_mu,
    prior_beta_sex_sd   = p$p_beta_sex_sd,
    prior_sigma_sd      = p$p_sigma_sd
  ))
  
  fits$M3[[pn]] <- run_fit(sm3, list(
    N = N, K = K,
    bill_length          = y_obs,
    bill_depth           = penguins$bill_depth_std,
    sex                  = penguins$sex_num,
    species              = penguins$species_num,
    prior_mu_alpha_mu    = p$p_mu_alpha_mu,
    prior_mu_alpha_sd    = p$p_mu_alpha_sd,
    prior_sigma_alpha_sd = p$p_sigma_alpha_sd,
    prior_beta_depth_mu  = p$p_beta_mu,
    prior_beta_depth_sd  = p$p_beta_sd,
    prior_beta_sex_mu    = p$p_beta_sex_mu,
    prior_beta_sex_sd    = p$p_beta_sex_sd,
    prior_sigma_sd       = p$p_sigma_sd
  ))
}
print(fits)
# ── Colour palette ────────────────────────────────────────────────────────────
col_base  <- "steelblue"
col_large <- "tomato"
col_serre <- "seagreen"
prior_colors <- c(base = col_base, diffuse = col_large, informative = col_serre)

# ── Robustness table ──────────────────────────────────────────────────────────
sigma_alpha_results <- data.frame(
  prior = names(priors),
  posterior_mean_mm = sapply(names(priors), function(pn)
    round(mean(rstan::extract(fits$M3[[pn]])$sigma_alpha), 2))
)
print(sigma_alpha_results)

# ── density overlay ───────────────────────────────────────────────────
plot_stat_overlay <- function(model_label, stat_fn, stat_name, x_label) {
  stat_obs <- stat_fn(y_obs)
  
  df_list <- lapply(names(priors), function(pn) {
    yr   <- rstan::extract(fits[[model_label]][[pn]])$y_rep
    vals <- apply(yr, 1, stat_fn)
    data.frame(value = vals, prior = pn)
  })
  df <- do.call(rbind, df_list)
  df$prior <- factor(df$prior, levels = c("base", "diffuse", "informative"))
  
  ggplot(df, aes(x = value, fill = prior, colour = prior)) +
    geom_density(alpha = 0.5, linewidth = 1) +         
    geom_vline(xintercept = stat_obs, colour = "black",
               linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = stat_obs, y = Inf,
             label = paste0("T(y) = ", round(stat_obs, 1), " mm"),
             hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = prior_colors) +
    scale_colour_manual(values = prior_colors) +
    guides(
      fill   = guide_legend(override.aes = list(alpha = 0.8, linewidth = 2)),
      colour = guide_legend(override.aes = list(alpha = 1,   linewidth = 2))
    ) +
    labs(
      title    = paste(model_label, "—", stat_name),
      subtitle = "Dashed line = observed value | Colours = prior regime",
      x        = x_label, y = "Density",
      fill = "Prior", colour = "Prior"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.key.size = unit(1.2, "lines"))
}

# ── min-max 2D scatter ────────────────────────────────────────────────
plot_minmax_overlay <- function(model_label) {
  obs_min <- min(y_obs)
  obs_max <- max(y_obs)
  
  df_list <- lapply(names(priors), function(pn) {
    yr <- rstan::extract(fits[[model_label]][[pn]])$y_rep
    data.frame(
      min_val = apply(yr, 1, min),
      max_val = apply(yr, 1, max),
      prior   = pn
    )
  })
  df <- do.call(rbind, df_list)
  df$prior <- factor(df$prior, levels = c("base", "diffuse", "informative"))
  
  ggplot() +
    geom_point(data = df,
               aes(x = min_val, y = max_val, colour = prior),
               alpha = 0.25, size = 1.2) +              
    geom_point(aes(x = obs_min, y = obs_max),
               colour = "black", size = 5, shape = 16) +
    annotate("text", x = obs_min, y = obs_max,
             label = "T(y)", hjust = -0.3, vjust = -0.5,
             size = 3.5, fontface = "bold") +
    scale_colour_manual(values = prior_colors) +
    guides(
      colour = guide_legend(
        override.aes = list(alpha = 1, size = 4)
      )
    ) +
    labs(
      title    = paste(model_label, "— Min vs Max"),
      subtitle = "Black dot = observed (min, max) | Coloured cloud = T(y_rep)",
      x = "min (mm)", y = "max (mm)", colour = "Prior"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.key.size = unit(1.2, "lines"))
}

# ── Affichage ─────────────────────────────────────────────────────────────────
for (m in c("M1", "M2", "M3")) {
  print(plot_stat_overlay(m, mean, "Mean", "Simulated mean (mm)"))
  print(plot_stat_overlay(m, min,  "Min",  "Simulated min (mm)"))
  print(plot_stat_overlay(m, max,  "Max",  "Simulated max (mm)"))
  print(plot_minmax_overlay(m))
}

# ── Sauvegarde ────────────────────────────────────────────────────────────────
for (stat_name in c("Mean", "Min", "Max")) {
  stat_fn  <- match.fun(tolower(stat_name))
  x_lab    <- paste0("Simulated ", tolower(stat_name), " (mm)")
  plot_row <- lapply(c("M1", "M2", "M3"), function(m)
    plot_stat_overlay(m, stat_fn, stat_name, x_lab))
  g <- arrangeGrob(grobs = plot_row, ncol = 3,
                   top = paste0("Sensitivity PPC — ", stat_name,
                                "\n3 prior regimes overlaid | base / diffuse / informative"))
  ggsave(paste0("fig_sensitivity_", tolower(stat_name), ".png"),
         plot = g, width = 15, height = 5, dpi = 150)
}

plot_row_2d <- lapply(c("M1", "M2", "M3"), plot_minmax_overlay)
g_2d <- arrangeGrob(grobs = plot_row_2d, ncol = 3,
                    top = "Sensitivity PPC — Min vs Max\n3 prior regimes overlaid")
ggsave("fig_sensitivity_minmax.png", plot = g_2d, width = 15, height = 5, dpi = 150)

message("Done — 4 figures saved.")

# =============================================================================
# NUMERICAL ROBUSTNESS SUMMARY
# =============================================================================

library(HDInterval) 

# Parameters to monitor per model
params_m1 <- c("alpha", "beta", "sigma")
params_m2 <- c("alpha[1]", "alpha[2]", "alpha[3]",
               "beta_depth[1]", "beta_depth[2]", "beta_depth[3]",
               "beta_sex", "sigma")
params_m3 <- c("mu_alpha", "sigma_alpha",
               "alpha[1]", "alpha[2]", "alpha[3]",
               "beta_depth", "beta_sex", "sigma")

extract_summary <- function(fit, params) {
  draws <- rstan::extract(fit)
  do.call(rbind, lapply(params, function(par) {
    # handle indexed params like "alpha[1]"
    match <- regmatches(par, regexec("^(\\w+)(?:\\[(\\d+)\\])?$", par))[[1]]
    var_name <- match[2]
    idx      <- if (nchar(match[3]) > 0) as.integer(match[3]) else NA
    vals <- if (is.na(idx)) draws[[var_name]] else draws[[var_name]][, idx]
    ci   <- hdi(vals, credMass = 0.90)
    data.frame(
      parameter = par,
      mean      = round(mean(vals), 3),
      sd        = round(sd(vals),   3),
      hpd_lo    = round(ci[1],      3),
      hpd_hi    = round(ci[2],      3)
    )
  }))
}

build_table <- function(model_label, params) {
  rows <- lapply(names(priors), function(pn) {
    s        <- extract_summary(fits[[model_label]][[pn]], params)
    s$prior  <- pn
    s$model  <- model_label
    s
  })
  df <- do.call(rbind, rows)
  
  # relative % change of the mean vs base
  base_means <- df[df$prior == "base", c("parameter", "mean")]
  names(base_means)[2] <- "base_mean"
  df <- merge(df, base_means, by = "parameter")
  df$pct_change <- round(100 * abs(df$mean - df$base_mean) /
                           (abs(df$base_mean) + 1e-9), 1)
  df$flag <- ifelse(df$pct_change > 10 & df$prior != "base", "***", "")
  df[order(df$parameter, df$prior), ]
}

# Build tables for all 3 models
tab_m1 <- build_table("M1", params_m1)
tab_m2 <- build_table("M2", params_m2)
tab_m3 <- build_table("M3", params_m3)

print_robustness <- function(tab, model_label) {
  message("\n===== ", model_label, " — Sensitivity summary =====")
  out <- tab[, c("model", "parameter", "prior",
                 "mean", "sd", "hpd_lo", "hpd_hi",
                 "pct_change", "flag")]
  print(out, row.names = FALSE)
}

print_robustness(tab_m1, "M1")
print_robustness(tab_m2, "M2")
print_robustness(tab_m3, "M3")

message("\n===== Parameters with largest prior sensitivity (non-base only) =====")
all_tabs <- rbind(tab_m1, tab_m2, tab_m3)
sensitive <- all_tabs[all_tabs$prior != "base", ]
sensitive  <- sensitive[order(-sensitive$pct_change), ]
print(sensitive[, c("model", "parameter", "prior",
                    "mean", "base_mean", "pct_change", "flag")],
      row.names = FALSE)

for (pn in names(priors)) {
  vals <- rstan::extract(fits$M3[[pn]])$sigma_alpha
  ci   <- hdi(vals, credMass = 0.90)
  excludes_zero <- ci[1] > 0
  message(sprintf("  %-12s  90%% HPD = [%.2f, %.2f]  excludes 0: %s",
                  pn, ci[1], ci[2],
                  ifelse(excludes_zero, "YES", "NO")))
}

for (pn in names(priors)) {
  # M2: 3 species slopes
  for (k in 1:3) {
    vals <- rstan::extract(fits$M2[[pn]])$beta_depth[, k]
    sp   <- c("Adelie", "Chinstrap", "Gentoo")[k]
    message(sprintf("  M2 %-12s %-10s  P(beta>0) = %.3f",
                    pn, sp, mean(vals > 0)))
  }
  # M3: shared slope
  vals <- rstan::extract(fits$M3[[pn]])$beta_depth
  message(sprintf("  M3 %-12s shared     P(beta>0) = %.3f",
                  pn, mean(vals > 0)))
}
