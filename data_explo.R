data <- readRDS("penguins.RDS")
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(rstan)
library(loo)
library(bayesplot) 
# ============================================================
# 1. DATA EXPLORATION
# ============================================================
##General case with all the penguins
summary(data)

ggplot(data, aes(x = bill_depth, y = bill_length)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  labs(
    title  = "Bill Depth vs Bill Length",
    x      = "Bill depth (mm)",
    y      = "Bill Length (mm)",
    color  = "Species",
    shape  = "Sex"
  )

modele_global <- lm(bill_length ~ bill_depth, data = data)
summary(modele_global)

##Statistics for each species
couleurs <- c("Adelie" = "#E07B54", "Chinstrap" = "#9B59B6", "Gentoo" = "#3A9AD9")

stats_species <- data %>%
  group_by(species) %>%
  summarise(
    n            = n(),
    mean_length  = round(mean(bill_length), 2),
    sd_length    = round(sd(bill_length), 2),
    mean_depth   = round(mean(bill_depth), 2),
    sd_depth     = round(sd(bill_depth), 2)
  )
print(stats_species)

#Scatter plot Bill Length vs Bill depth
ggplot(data, aes(x = bill_depth, y = bill_length, color = species)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = species), method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = couleurs) +
  labs(
    title  = "Bill Depth vs Bill Length (by species)",
    x      = "Bill depth (mm)",
    y      = "Bill Length (mm)",
    color  = "Species",
    shape  = "Sex"
  )

##Statistics by species and by sex
stats_sex <- data %>%
  group_by(species, sex) %>%
  summarise(
    n           = n(),
    mean_length = round(mean(bill_length), 2),
    sd_length   = round(sd(bill_length), 2),
    .groups = "drop"
  )
print(stats_sex)

#Scatter plot Bill Length vs Bill Depth for each sex
ggplot(data, aes(x = bill_depth, y = bill_length, color = species)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = couleurs) +
  facet_wrap(~ sex) +
  labs(
    title = "Relation between Bill Depth / Length depending on the sex",
    x     = "Bill Depth (mm)",
    y     = "Bill Length (mm)",
    color = "Species"
  )

#Boxplot Bill Length ~ species + sex
ggplot(data, aes(x = species, y = bill_length, fill = sex)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Bill length depending on species and sex",
       x = "Species", y = "Bill length (mm)")

#Boxplot Bill Depth ~ species + sex
ggplot(data, aes(x = species, y = bill_depth, fill = sex)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Bill depth depending on species and sex",
       x = "Species", y = "Bill depth (mm)")



# ── Figure : histogram per species ──────────────────
p_hist <- ggplot(penguins, aes(x = bill_length, fill = species)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, colour = "white", alpha = 0.8) +
  geom_density(colour = "black", linewidth = 0.8, fill = NA) +
  facet_wrap(~species, ncol = 3) +
  scale_fill_manual(values = species_colors) +
  labs(
    title    = "Distribution of bill length by species",
    subtitle = "Histogram + kernel density — bell-shaped profile within each group",
    x        = "Bill length (mm)",
    y        = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

print(p_hist)
ggsave("fig_hist_by_species.png", p_hist, width = 12, height = 4, dpi = 150)


# ── Figure : per species and per sexe ─────────────────────────────────────────


p_hist_sex <- ggplot(penguins, aes(x = bill_length, fill = species)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1.5, colour = "white", alpha = 0.8) +
  geom_density(colour = "black", linewidth = 0.7, fill = NA) +
  facet_grid(sex ~ species) +
  scale_fill_manual(values = species_colors) +
  labs(
    title    = "Distribution of bill length by species and sex",
    subtitle = "Bell-shaped profile within each species-sex combination",
    x        = "Bill length (mm)",
    y        = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

print(p_hist_sex)
ggsave("fig_hist_by_species_sex.png", p_hist_sex, width = 12, height = 7, dpi = 150)
