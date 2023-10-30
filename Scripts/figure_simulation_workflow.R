# Scaling analysis
# Figure illustrating simulation approach

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(mgcv)
require(tidymv)
require(scales) # for numeric axis labels with commas
source("Scripts/truncated_lnorm_functions.R")


# Load data --------------

# Landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv") %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out NW and WC fires
fire_metrics_NW <- filter(fire_metrics, NW_Cascadia == 0) 
fire_metrics_NWC <- filter(fire_metrics, NW_Cascadia == 1) 

# Pull out high severity FRGs
fire_metrics_high <- filter(fire_metrics_NW, Fire_Regime == "High")




# Simulation approach for total core area -----

# Combine NW high-severity data with NWC data
fire_metrics_NWC_high <- bind_rows(fire_metrics_NWC, fire_metrics_high)

# Fit quantile regression model to core area
mod_core <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)

# Extract core area fitted values
fit_core <- fire_metrics_NWC_high %>% 
  drop_na(log_total_core) %>% bind_cols( as_tibble(mod_core$fitted.values) )

# Pivot longer to plot all quantiles
fit_core_long <- fit_core %>%
  pivot_longer(cols = (length(fit_core)-98):length(fit_core),
               names_to = "tau",
               values_to = "fit")

# Plot all quantiles
ggplot() +
  scale_y_continuous(breaks = c(0, 2, 4),
                     labels = c("1", "100", "10,000"),
                     name = "Total core area (ha)") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_line(data = fit_core_long, mapping =
              aes(x = log_fire_area, y = fit, group = tau), color = "gray30", alpha = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("Figures/figure2_core.png",
       width = 2.7, height = 2.5, units = "in", dpi = 320)




# Simulation approach for distance-to-seed distributions -------------

# Combine NW high-severity data with NWC data
fire_metrics_NWC_high <- bind_rows(fire_metrics_NWC, fire_metrics_high) %>%
  drop_na(log_SDC)

# Step 1: Fit quantile regression model to SDC parameter
mod_SDC <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)

# Extract SDC fitted values
fit_SDC <- fire_metrics_NWC_high %>% 
  drop_na(log_SDC) %>% bind_cols( as_tibble(mod_SDC$fitted.values) )

# Pivot longer to plot all percentiles
fit_SDC_long <- fit_SDC %>%
  pivot_longer(cols = (length(fit_SDC)-98):length(fit_SDC),
               names_to = "tau",
               values_to = "fit")

# Plot all quantiles
ggplot() +
  scale_y_continuous(limits = c(-3.1, -1.8),
                     breaks = c(-3, -2.5, -2),
                     labels = c("0.001", "0.003", "0.01"),
                     name = "SDC") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_line(data = fit_SDC_long, mapping =
              aes(x = log_fire_area, y = fit, group = tau), color = "gray30", alpha = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("Figures/figure2_SDC.png",
       width = 2.7, height = 2.5, units = "in", dpi = 320)


# Step 2: Fit gam model to high-severity & forested proportion (as function of fire size & SDC)
mod_HS_forest_prop <- gam(HS_forest_prop ~ s(log_fire_area) + s(log_SDC), 
                          family = betar(link="logit"), data = fire_metrics_NWC_high)

inv_logit <- function(x){1/(1+exp(-x))}

mod_HS_forest_prop_p <- predict_gam(mod_HS_forest_prop)

mod_HS_forest_prop_p %>%
  ggplot(aes(log_fire_area, log_SDC, z = inv_logit(fit))) +
  geom_raster(aes(fill = inv_logit(fit))) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(option = "magma", name = expression(P[FHS])) +
  scale_y_continuous(breaks = c(-3, -2.5, -2),
                     labels = c("0.001", "0.003", "0.01"),
                     name = "SDC") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

ggsave("Figures/figure2_SDC_FHSprop.png",
       width = 3.3, height = 2.5, units = "in", dpi = 320)


# Step 3: Calculate distance-to-seed distribution

# An example with specified parameters

# Distance bins
dist <- seq(0, 1000, by=10)
# Calculate forested & high-severity area
# as inverse cumulative proportions, using modified logistic function
sim_SDC <- 0.003
P <- 1 / (10^(sim_SDC*dist))
# Convert inverse cum. proportions to inverse cumulative areas
sim_HS_forest_area = 20000 * 0.5
A <- round(P * sim_HS_forest_area)

DTS_dist <- tibble(dist = dist, P = P, A = A)

ggplot(DTS_dist) +
  geom_line(mapping = aes(x = dist, y = A / 1000)) +
  scale_y_continuous(breaks = c(0, 5, 10),
                     name = "Inverse cumulative area\n(1,000 ha)") +
  scale_x_continuous(breaks = c(0, 500, 1000),
                     labels = c("0", "500", "1,000"),
                     name = "Distance to seed (m)") +
  annotate("text", x = 300, y = 8, hjust = 0, size = 4, color = "gray20",
           label = 'Fire size = 2e4 ha') +
  annotate("text", x = 300, y = 7, hjust = 0, size = 4, color = "gray20",
           label = 'SDC = 0.003') +
  annotate("text", x = 300, y = 6, hjust = 0, size = 4, color = "gray20",
           label = expression(P[FHS] ~ '= 0.5')) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5, unit = "pt"))

ggsave("Figures/figure2_DTS.png",
       width = 3, height = 2.5, units = "in", dpi = 320)



# Simulation approach for patch size distributions  ------

# Combine NW high-severity data with NWC data
fire_metrics_NWC_high <- bind_rows(fire_metrics_NWC, fire_metrics_high) %>%
  drop_na(psi)

# Step 1: Fit quantile regression model to psi parameter

# Psi parameter converges with increasing fire size.
# However, we see unrealistic non-monotonic fitted curves for highest upper quantiles.
# Thus, need to impose monotonicity constraints for high upper quantiles.

# All quantiles, no monotonicity constraints imposed
mod_psi_aq <- gcrq(psi ~ ps(log_fire_area),
                   lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)

# Upper quantiles with monotonicity constraints imposed
mod_psi_uq <- gcrq(psi ~ ps(log_fire_area, monotone = -1),
                   lambda0 = 5, tau = seq(0.5,0.99,0.01), data = fire_metrics_NWC_high)

# Extract psi fitted values without any monotonicity contstraints
fit_psi <- fire_metrics_NWC_high %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_psi_aq$fitted.values) ) 

# Pivot longer to plot all quantiles
fit_psi_long <- fit_psi %>%
  pivot_longer(cols = (length(fit_psi)-98):length(fit_psi),
               names_to = "tau",
               values_to = "fit")

# Plot all quantiles
# This plot demonstrates the unrealistically large values of psi 
# at upper quantiles & large fire sizes
ggplot() +
  scale_y_continuous(name = "\u03C8") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_line(data = fit_psi_long, mapping =
              aes(x = log_fire_area, y = fit, group = tau), color = "gray30", alpha = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# So instead, do this:
# Extract psi fitted values with monotonicity constraints for percentiles >= 96th ptile
# and combined with fitted values without monotonicity constraints for percentiles < 96th ptile
fit_psi <- fire_metrics_NWC_high %>% drop_na(psi) %>% 
  bind_cols( as_tibble(mod_psi_aq$fitted.values) ) %>%
  select(-`0.96`, -`0.97`, -`0.98`, -`0.99`) %>%
  bind_cols( as_tibble(mod_psi_uq$fitted.values) %>%
               select(`0.96`, `0.97`, `0.98`, `0.99`) )

# Pivot longer to plot all quantiles
fit_psi_long <- fit_psi %>%
  pivot_longer(cols = (length(fit_psi)-98):length(fit_psi),
               names_to = "tau",
               values_to = "fit")

# Plot all quantiles
ggplot() +
  scale_y_continuous(name = "\u03C8") +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("1", "10", "100"),
                     name = "Fire size (1,000 ha)") +
  geom_line(data = fit_psi_long, mapping =
              aes(x = log_fire_area, y = fit, group = tau), color = "gray30", alpha = 0.2) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("Figures/figure2_psi.png",
       width = 2.7, height = 2.5, units = "in", dpi = 320)



# Step 2: Fit gam model to beta parameter (as function of psi)
mod_beta <- gam(beta ~ s(psi), data = fire_metrics_NWC_high, gamma = 5)

mod_beta_p <- predict_gam(mod_beta)

mod_beta_p %>%
  ggplot(aes(psi, fit)) +
  scale_y_continuous(name = "\u03B2") + # beta
  scale_x_continuous(name = "\u03C8", limits = c(-0.19, 0.9)) + # psi
  geom_smooth_ci(ci_alpha = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("Figures/figure2_psi_beta.png",
       width = 2.7, height = 2.5, units = "in", dpi = 320)



# Step 3: Fit gam model to HS proportion (as function of psi, beta, fire size)
mod_HS_prop <- gam(HS_prop ~ s(psi) + s(beta) + s(log_fire_area), 
                   family = betar(link="logit"), data = fire_metrics_NWC_high)

mod_HS_prop_p <- predict_gam(mod_HS_prop)

mod_HS_prop_p %>%
  ggplot(aes(psi, beta, z = inv_logit(fit))) +
  geom_raster(aes(fill = inv_logit(fit))) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(option = "magma",
                       breaks = c(0.2, 0.4, 0.6, 0.8),
                       name = expression(P[HS])) +
  scale_y_continuous(breaks = c(-1,0,1,2),
                     name = "\u03B2") + # beta
  scale_x_continuous(breaks = c(0, 0.3, 0.6),
                     name = "\u03C8") + # psi
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))

ggsave("Figures/figure2_psi_beta_HSprop.png",
       width = 3.3, height = 2.5, units = "in", dpi = 320)


# Step 4: Randomly sample from patch size distribution

# Function to calculate inverse empirical CDF
inv_CDF_fun <- function(patch_dist){
  x_n <- length(patch_dist)
  x_ecdf <- c()
  for(x in 1:x_n){
    x_ecdf <- c(x_ecdf,
                sum(patch_dist >= patch_dist[x]) / x_n)
  }
  return(x_ecdf)
}

# Function to simulate patch size distribution, 
# given the following parameters: fire size, psi, beta, HS proportion
sim_patch_fun <- function(log_fire_size, psi, beta, HS_prop){
  
  # Calculate high-severity area
  sim_HS_area <- HS_prop * 10^log_fire_size
  
  # Sample from truncated lnorm distribution until total equals HS area
  sim_dist <- c()
  
  while(sum(sim_dist) < sim_HS_area){
    
    # Draw patch
    patch <- r_truncated_ln(n = 1, xmin = 1, xmax = 10^log_fire_size,
                            beta = beta, psi = psi)
    
    # Check that patch won't put distribution over HS area limit;
    # if not, add patch to distribution;
    # if so, adjust the final patch size to be equal to the 
    # difference between current HS area and target HS area
    if(sum(sim_dist) + patch <= sim_HS_area){
      sim_dist <- c(sim_dist, patch)
    } else {
      sim_dist <- c(sim_dist, sim_HS_area - sum(sim_dist))
    }
    
  }
  
  return(sim_dist)
}

# An example with specified parameters
set.seed(321)
sim_patch_dist <- tibble( patch = sim_patch_fun(log_fire_size = log10(20000),
                                                psi = 0.05, 
                                                beta = 1.4, 
                                                HS_prop = 0.5)) %>%
  mutate( inv_ecdf = inv_CDF_fun(patch))


ggplot(sim_patch_dist, aes(x = patch, y = inv_ecdf)) + 
  geom_point(shape = 1, alpha = 0.6) +
  labs(x = "Patch size (ha)", y = "Inv. cumulative probability") +
  scale_x_log10(limits = c(1, 20000),
                breaks = c(1, 100, 10000),
                labels = comma) + 
  scale_y_log10(limits = c(1e-3, 1.2)) +
  annotate("text", x = 30, y = 10^0, hjust = 0, size = 4, color = "gray20",
           label = "Fire size = 2e4 ha") +
  annotate("text", x = 30, y = 10^-0.35, hjust = 0, size = 4, color = "gray20",
           label = expression(paste(psi, " = 0.05, ", beta, " = 1.4"))) +
  annotate("text", x = 30, y = 10^-0.7, hjust = 0, size = 4, color = "gray20",
           label = expression(P[HS] ~ '= 0.5')) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave("Figures/figure2_patch_dist.png",
       width = 3, height = 2.5, units = "in", dpi = 320)
