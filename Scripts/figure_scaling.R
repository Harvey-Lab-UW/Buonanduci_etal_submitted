# Scaling analysis
# Figure illustrating quantile regression scaling relationships

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(cowplot)


# Load data --------------

# Landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv")%>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Pull out NW and NWC fires
fire_metrics_NW <- filter(fire_metrics, NW_Cascadia == 0) 
fire_metrics_NWC <- filter(fire_metrics, NW_Cascadia == 1) 

# Pull out unique fire regimes
fire_metrics_high <- filter(fire_metrics_NW, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics_NW, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics_NW, Fire_Regime == "Low")

# Create facet labels
FRG_names <- c(
  `High` = "High-severity regime",
  `Mixed` = "Mixed-severity regime",
  `Low` = "Low-severity regime"
)

# Color for plotting
mycol <- "#0d3a67"


# Area-wtd mean patch size ----
mod_high <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_patch_area_AW_mean) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))


# Plotting
p1 <- ggplot(pred) +
  facet_grid( ~ Fire_Regime, labeller = as_labeller(FRG_names)) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = " \nPatch size:\nArea-wtd mean (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`), fill = mycol, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`), fill = mycol, alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.25`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.75`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.95`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.5`), color = mycol) +
  geom_point(aes(x = log_fire_area, y = log_patch_area_AW_mean), color = "gray30", size = 0.6, alpha = 0.1) +
  geom_point(data = select(fire_metrics_NWC, -Fire_Regime), 
             mapping = aes(x = log_fire_area, y = log_patch_area_AW_mean), shape = 1, size = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 9),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Patch distribution: beta ----
mod_high <- gcrq(beta ~ ps(log_fire_area),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions
pred_high <- fire_metrics_high %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_high$fitted.values) ) 
pred_mix <- fire_metrics_mix %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_mix$fitted.values) ) 
pred_low <- fire_metrics_low %>% drop_na(beta) %>% 
  bind_cols( as_tibble(mod_low$fitted.values) ) 

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p2 <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(name = " \nPatch size:\n\u03B2 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`), fill = mycol, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`), fill = mycol, alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.25`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.75`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.95`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.5`), color = mycol) +
  geom_point(aes(x = log_fire_area, y = beta), color = "gray30", size = 0.6, alpha = 0.1) +
  geom_point(data = select(fire_metrics_NWC, -Fire_Regime), 
             mapping = aes(x = log_fire_area, y = beta), shape = 1, size = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Patch distribution: psi ----
mod_high <- gcrq(psi ~ ps(log_fire_area),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(psi ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(psi ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions
pred_high <- fire_metrics_high %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_high$fitted.values) ) 
pred_mix <- fire_metrics_mix %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_mix$fitted.values) ) 
pred_low <- fire_metrics_low %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_low$fitted.values) ) 

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p3 <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(limits = c(-0.21, 1),
                     name = " \nPatch size:\n\u03C8 parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`), fill = mycol, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`), fill = mycol, alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.25`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.75`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.95`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.5`), color = mycol) +
  geom_point(aes(x = log_fire_area, y = psi), color = "gray30", size = 0.6, alpha = 0.1) +
  geom_point(data = select(fire_metrics_NWC, -Fire_Regime), 
             mapping = aes(x = log_fire_area, y = psi), shape = 1, size = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



# Core area ----
mod_high <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p4 <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = " \nPatch structure:\nTotal core area (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`), fill = mycol, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`), fill = mycol, alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.25`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.75`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.95`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.5`), color = mycol) +
  geom_point(aes(x = log_fire_area, y = log_total_core), color = "gray30", size = 0.6, alpha = 0.1) +
  geom_point(data = select(fire_metrics_NWC, -Fire_Regime), 
             mapping = aes(x = log_fire_area, y = log_total_core), shape = 1, size = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



# Seed decay coefficient ----
mod_high <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Predictions
pred_high <- fire_metrics_high %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_high$fitted.values) )
pred_mix <- fire_metrics_mix %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_mix$fitted.values) )
pred_low <- fire_metrics_low %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_low$fitted.values) )

pred <- bind_rows(pred_high, pred_mix, pred_low) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High")))

# Plotting
p5 <- ggplot(pred) +
  facet_grid( ~ Fire_Regime) +
  scale_y_continuous(breaks = c(-3, -2.5, -2, -1.5),
                     labels = c("0.001", "0.003", "0.01", "0.03"),
                     name = " \nPatch structure:\nSDC parameter") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1,000", "10,000", "100,000"),
                     name = "Fire size (ha)") +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.05`, ymax = `0.95`), fill = mycol, alpha = 0.15) +
  geom_ribbon(aes(x = log_fire_area, ymin = `0.25`, ymax = `0.75`), fill = mycol, alpha = 0.2) +
  geom_line(aes(x = log_fire_area, y = `0.05`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.25`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.75`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.95`), lwd = 0.1, color = mycol) +
  geom_line(aes(x = log_fire_area, y = `0.5`), color = mycol) +
  geom_point(aes(x = log_fire_area, y = log_SDC), color = "gray30", size = 0.6, alpha = 0.1) +
  geom_point(data = select(fire_metrics_NWC, -Fire_Regime), 
             mapping = aes(x = log_fire_area, y = log_SDC), shape = 1, size = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Combine plots ---------

plot_grid(p1, p2, p3, p4, p5, labels = c("(a)", "(b)", "(c)", "(d)", "(e)"), 
          label_size = 10,
          ncol = 1, align = "v", rel_heights = c(1.2, 1, 1, 1, 1.3))

ggsave(paste0("Figures/figure3_scaling.png"),
       width = 6.5, height = 7, units = "in", dpi = 320)






