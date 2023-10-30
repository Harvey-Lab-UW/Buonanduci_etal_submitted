# Scaling analysis
# Evaluate prediction error of regime-specific quantile regression models
# to evaluate fit of models for Northwestern Cascadia

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(cowplot)


# Load data ----------------

# Landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv") %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out NW and WC fires
fire_metrics_NW <- filter(fire_metrics, NW_Cascadia == 0) 
fire_metrics_NWC <- filter(fire_metrics, NW_Cascadia == 1) 

# Pull out unique regimes
fire_metrics_high <- filter(fire_metrics_NW, Fire_Regime == "High")
fire_metrics_mix <- filter(fire_metrics_NW, Fire_Regime == "Mixed")
fire_metrics_low <- filter(fire_metrics_NW, Fire_Regime == "Low")


# Quantile loss functions for evaluating NW Cascadia fit ----
loss_fun <- function(tau, obs, pred){
  mean( pmax(tau * (obs - pred), (tau - 1) * (obs - pred)) )
}

model_loss_fun <- function(model, df, var){
  # Drop any NAs for response variable
  df <- df[ !is.na(df[[var]]), ]
  # Make predictions from model
  pred <- as_tibble(predict.gcrq(model, newdata = as.data.frame(df)))
  # If only 1 value of tau, assign name
  if (length(pred) == 1){ names(pred) <- as.character(model$taus) }
  # Calculate loss for each quantile (tau)
  loss <- c()
  for (tau in names(pred)){
    loss <- c(loss, loss_fun(tau = as.numeric(tau), obs = df[[var]],
                             pred = pred[[tau]]) )
  }
  # Return values of tau and calculated loss
  return(tibble(tau = as.numeric(names(pred)),
                loss = loss))
}



# Area-wtd mean patch size ----
mod_high <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_patch_area_AW_mean ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Quantile loss calculations
loss_high <- model_loss_fun(mod_high, df = fire_metrics_NWC, var = "log_patch_area_AW_mean") %>%
  mutate(FRG_max4 = "High")
loss_mix <- model_loss_fun(mod_mix, df = fire_metrics_NWC, var = "log_patch_area_AW_mean") %>%
  mutate(FRG_max4 = "Mixed")
loss_low <- model_loss_fun(mod_low, df = fire_metrics_NWC, var = "log_patch_area_AW_mean") %>%
  mutate(FRG_max4 = "Low")

loss <- bind_rows(loss_high, loss_mix, loss_low) %>%
  mutate(FRG_max4 = factor(FRG_max4, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau))

p1 = ggplot(loss, aes(y = loss, x = tau, group = FRG_max4, color = FRG_max4)) +
  geom_line() + geom_point(shape = 1, size = 2) +
  scale_color_viridis_d(name = "Fire regime", direction = -1, option = "inferno", end = 0.8) +
  scale_y_continuous(name = "Quantile loss", limits = c(0, max(loss$loss))) +
  scale_x_discrete(name = "Quantile", labels = c('0.05', '0.25', '0.50', '0.75', '0.95')) +
  ggtitle("Area-weighted mean patch size") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.position = "none")




# Patch distribution: beta ----
mod_high <- gcrq(beta ~ ps(log_fire_area),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(beta ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Quantile loss calculations
loss_high <- model_loss_fun(mod_high, df = fire_metrics_NWC, var = "beta") %>%
  mutate(FRG_max4 = "High")
loss_mix <- model_loss_fun(mod_mix, df = fire_metrics_NWC, var = "beta") %>%
  mutate(FRG_max4 = "Mixed")
loss_low <- model_loss_fun(mod_low, df = fire_metrics_NWC, var = "beta") %>%
  mutate(FRG_max4 = "Low")

loss <- bind_rows(loss_high, loss_mix, loss_low) %>%
  mutate(FRG_max4 = factor(FRG_max4, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau))

p2 = ggplot(loss, aes(y = loss, x = tau, group = FRG_max4, color = FRG_max4)) +
  geom_line() + geom_point(shape = 1, size = 2) +
  scale_color_viridis_d(name = "Fire regime", direction = -1, option = "inferno", end = 0.8) +
  scale_y_continuous(name = "Quantile loss", limits = c(0, max(loss$loss))) +
  scale_x_discrete(name = "Quantile", labels = c('0.05', '0.25', '0.50', '0.75', '0.95')) +
  ggtitle("\u03B2 parameter") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.position = "none")



# Patch distribution: psi ----
mod_high <- gcrq(psi ~ ps(log_fire_area),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(psi ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(psi ~ ps(log_fire_area),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Quantile loss calculations
loss_high <- model_loss_fun(mod_high, df = fire_metrics_NWC, var = "psi") %>%
  mutate(FRG_max4 = "High")
loss_mix <- model_loss_fun(mod_mix, df = fire_metrics_NWC, var = "psi") %>%
  mutate(FRG_max4 = "Mixed")
loss_low <- model_loss_fun(mod_low, df = fire_metrics_NWC, var = "psi") %>%
  mutate(FRG_max4 = "Low")

loss <- bind_rows(loss_high, loss_mix, loss_low) %>%
  mutate(FRG_max4 = factor(FRG_max4, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau))

p3 =  ggplot(loss, aes(y = loss, x = tau, group = FRG_max4, color = FRG_max4)) +
  geom_line() + geom_point(shape = 1, size = 2) +
  scale_color_viridis_d(name = "Fire regime", direction = -1, option = "inferno", end = 0.8) +
  scale_y_continuous(name = "Quantile loss", limits = c(0, max(loss$loss))) +
  scale_x_discrete(name = "Quantile", labels = c('0.05', '0.25', '0.50', '0.75', '0.95')) +
  ggtitle("\u03C8 parameter") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.position = "none")


# Core area ----
mod_high <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Quantile loss calculations
loss_high <- model_loss_fun(mod_high, df = fire_metrics_NWC, var = "log_total_core") %>%
  mutate(FRG_max4 = "High")
loss_mix <- model_loss_fun(mod_mix, df = fire_metrics_NWC, var = "log_total_core") %>%
  mutate(FRG_max4 = "Mixed")
loss_low <- model_loss_fun(mod_low, df = fire_metrics_NWC, var = "log_total_core") %>%
  mutate(FRG_max4 = "Low")

loss <- bind_rows(loss_high, loss_mix, loss_low) %>%
  mutate(FRG_max4 = factor(FRG_max4, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau))

p4 = ggplot(loss, aes(y = loss, x = tau, group = FRG_max4, color = FRG_max4)) +
  geom_line() + geom_point(shape = 1, size = 2) +
  scale_color_viridis_d(name = "Fire regime", direction = -1, option = "inferno", end = 0.8) +
  scale_y_continuous(name = "Quantile loss", limits = c(0, max(loss$loss))) +
  scale_x_discrete(name = "Quantile", labels = c('0.05', '0.25', '0.50', '0.75', '0.95')) +
  ggtitle("Total core area") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.position = "none")


# Distance-to-seed: SDC parameter ----
mod_high <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                 lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_high)
mod_mix <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_mix)
mod_low <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = c(0.05,0.25,0.5,0.75,0.95), data = fire_metrics_low)

# Quantile loss calculations
loss_high <- model_loss_fun(mod_high, df = fire_metrics_NWC, var = "log_SDC") %>%
  mutate(FRG_max4 = "High")
loss_mix <- model_loss_fun(mod_mix, df = fire_metrics_NWC, var = "log_SDC") %>%
  mutate(FRG_max4 = "Mixed")
loss_low <- model_loss_fun(mod_low, df = fire_metrics_NWC, var = "log_SDC") %>%
  mutate(FRG_max4 = "Low")

loss <- bind_rows(loss_high, loss_mix, loss_low) %>%
  mutate(FRG_max4 = factor(FRG_max4, levels = c("Low", "Mixed", "High"))) %>%
  mutate(tau = factor(tau))

p5 = ggplot(loss, aes(y = loss, x = tau, group = FRG_max4, color = FRG_max4)) +
  geom_line() + geom_point(shape = 1, size = 2) +
  scale_color_viridis_d(name = "Fire regime", direction = -1, option = "inferno", end = 0.8) +
  scale_y_continuous(name = "Quantile loss", limits = c(0, max(loss$loss))) +
  scale_x_discrete(name = "Quantile", labels = c('0.05', '0.25', '0.50', '0.75', '0.95')) +
  ggtitle("SDC parameter") +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.position = "bottom")

# Combine plots ---------
plot_grid(p1, p2, p3, p4, p5, 
          ncol = 1, align = "v", label_size = 10,
          rel_heights = c( 1, 1, 1, 1, 1.3),
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)"))

ggsave("Figures/figureS1_NWC_quantile_loss.png", bg = "white",
       width = 5, height = 8, units = "in", dpi = 320)
