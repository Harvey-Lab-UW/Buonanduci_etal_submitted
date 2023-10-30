# Scaling analysis
# Simulate patch size distributions

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(mgcv)
require(gratia) # for simulating from gam models
require(scales) # for numeric axis labels with commas
require(cowplot)
source("Scripts/truncated_lnorm_functions.R")


# Load data --------------

# Landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv") %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out NW and WC fires
fire_metrics_NW <- filter(fire_metrics, NW_Cascadia == 0) 
fire_metrics_NWC <- filter(fire_metrics, NW_Cascadia == 1) 

# Pull out high severity regime
fire_metrics_high <- filter(fire_metrics_NW, Fire_Regime == "High")

# Combine NW high-severity regime fires with NW Cascadia fires
fire_metrics_NWC_high <- bind_rows(fire_metrics_NWC, fire_metrics_high) %>%
  drop_na(beta)



# Patch size distribution models --------

# Psi parameter converges with increasing fire size.
# However, we see unrealistic non-monotonic fitted curves for highest upper quantiles.
# Thus, need to impose monotonicity constraints for high upper quantiles.

# All quantiles
mod_psi_aq <- gcrq(psi ~ ps(log_fire_area),
                   lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)
# Upper quantiles
mod_psi_uq <- gcrq(psi ~ ps(log_fire_area, monotone = -1),
                   lambda0 = 5, tau = seq(0.5,0.99,0.01), data = fire_metrics_NWC_high)

# Predictions at fire sizes 1000, 10000, 100000
pred_psi_3 <- c( predict.gcrq(mod_psi_aq, newdata = data.frame("log_fire_area" = 3))[1:95],
                 predict.gcrq(mod_psi_uq, newdata = data.frame("log_fire_area" = 3))[47:50]) 
pred_psi_4 <- c( predict.gcrq(mod_psi_aq, newdata = data.frame("log_fire_area" = 4))[1:95],
                 predict.gcrq(mod_psi_uq, newdata = data.frame("log_fire_area" = 4))[47:50]) 
pred_psi_5 <- c( predict.gcrq(mod_psi_aq, newdata = data.frame("log_fire_area" = 5))[1:95],
                 predict.gcrq(mod_psi_uq, newdata = data.frame("log_fire_area" = 5))[47:50]) 

# Fit gam model beta parameter (as function of psi)
mod_beta <- gam(beta ~ s(psi), data = fire_metrics_NWC_high)

# Fit gam model to HS prop (as function of psi, beta, fire size)
mod_HS_prop <- gam(HS_prop ~ s(psi) + s(beta) + s(log_fire_area), 
                   family = betar(link="logit"), data = fire_metrics_NWC_high)


# Patch size distribution simulation ---------

# !! Warning !!: patch size distribution simulation can be very slow
# Skip down to line of code that loads previously simulated distributions if desired

# Simulate patch size distributions under three scenarios:
# 10 fires, size 1e5 ha
# 100 fires, size 1e4 ha
# 1000 fires, size 1e3 ha

# Step 1: Simulate psi, randomly generated from quantile regression
# Step 2: Simulate beta from psi
# Step 3: Simulate HS_prop from psi, beta, and fire size
# Step 3: Use psi, beta, and HS_prop to generate patch size distribution
#         from each hypothetical fire, ensuring that total area of patches 
#         equals HS_area

# Function to simulate patch size distribution
sim_patch_dist <- function(log_fire_size, pred_psi, mod_beta, mod_HS_prop){
  
  # Simulate psi parameter
  sim_psi <- pred_psi[round(runif(1, 1, 99), 0)]
  
  # Simulate beta parameter
  sim_beta <- simulate(mod_beta,
                       newdata = data.frame(log_fire_area = log_fire_size,
                                            psi = sim_psi)) %>% as.numeric()
  
  # Simulate high-severity proportion
  sim_HS_prop <- simulate(mod_HS_prop, nsim = 1,
                          newdata = data.frame(log_fire_area = log_fire_size,
                                               psi = sim_psi, beta = sim_beta)) %>% as.numeric()
  
  # Calculate high-severity area
  sim_HS_area <- sim_HS_prop * 10^log_fire_size
  
  # Sample from truncated lnorm distribution until total equals HS area
  sim_dist <- c()
  
  while(sum(sim_dist) < sim_HS_area){
    
    # Draw patch
    patch <- r_truncated_ln(n = 1, xmin = 1, xmax = 10^log_fire_size,
                            beta = sim_beta, psi = sim_psi)
    
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


# 10 fires, size 1e5 ha
sim_patch5 <- tibble()

for (i in 1:100){ # iterations
  for (f in 1:10){ # fires
    sim_patch5 <- bind_rows(sim_patch5,
                            tibble( patch = sim_patch_dist(log_fire_size = 5, 
                                                           pred_psi = pred_psi_5, 
                                                           mod_beta, mod_HS_prop),
                                    f = f, i = i,
                                    fire_n_size = "10 x 100,000 ha"))
  }
}

# 100 fires, size 1e4 ha
sim_patch4 <- tibble()

for (i in 1:100){ # iterations
  for (f in 1:100){ # fires
    sim_patch4 <- bind_rows(sim_patch4,
                            tibble( patch = sim_patch_dist(log_fire_size = 4, 
                                                           pred_psi = pred_psi_4, 
                                                           mod_beta, mod_HS_prop),
                                    f = f, i = i,
                                    fire_n_size = "100 x 10,000 ha"))
  }
}

# 1000 fires, size 1e3 ha
sim_patch3 <- tibble()

for (i in 1:100){ # iterations
  for (f in 1:1000){ # fires
    sim_patch3 <- bind_rows(sim_patch3,
                            tibble( patch = sim_patch_dist(log_fire_size = 3, 
                                                           pred_psi = pred_psi_3, 
                                                           mod_beta, mod_HS_prop),
                                    f = f, i = i,
                                    fire_n_size = "1,000 x 1,000 ha"))
  }
}

# Combined tibble
sim_patch <- bind_rows(sim_patch5, sim_patch4, sim_patch3) %>%
  mutate(fire_n_size = factor(fire_n_size, levels = c("10 x 100,000 ha",
                                                      "100 x 10,000 ha",
                                                      "1,000 x 1,000 ha")))

# Calculate inverse empirical cumulative densities and areas
sim_patch_ecdf <- tibble()

for(fns in levels(sim_patch$fire_n_size)){
  for(d in unique(sim_patch$i)){
    
    patch_dist <- filter(sim_patch, fire_n_size == fns) %>%
      filter(i == d) %>%
      arrange(desc(patch))
    patch_dist <- patch_dist$patch
    
    # Inverse cumulative areas
    x_cumsum <- cumsum(patch_dist)
    
    # Inverse cumulative densities
    x_n <- length(patch_dist)
    x_ecdf <- c()
    for(x in 1:x_n){
      x_ecdf <- c(x_ecdf,
                  sum(patch_dist >= patch_dist[x]) / x_n)
    }
    
    sim_patch_ecdf <- bind_rows(sim_patch_ecdf,
                                tibble(patch = patch_dist,
                                       inv_ecdf = x_ecdf,
                                       inv_ca = x_cumsum / sum(patch_dist),
                                       fire_n_size = fns,
                                       i = d))
  }
}

# Write simulated patches to file
write_csv(sim_patch_ecdf, "Data/simulated_patches.csv")

# Use this line to read in previously simulated patches
#sim_patch_ecdf <- read_csv("Data/simulated_patches.csv")

sim_patch_ecdf <- sim_patch_ecdf %>% 
  mutate(fire_n_size = factor(fire_n_size, levels = c("10 x 100,000 ha",
                                                      "100 x 10,000 ha",
                                                      "1,000 x 1,000 ha")))

# Plotting ---------

# Inverse cumulative probability
p1 = ggplot(sim_patch_ecdf, aes(x = patch, y = inv_ecdf,
                                color = fire_n_size, group = paste(i, fire_n_size))) + 
  geom_line(alpha = 0.5) + 
  scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.9, direction = -1,
                        name = "Number and size of fires") +
  labs(x = "Patch size (ha)", y = "Inverse cumulative probability (P[X ≥ x])") +
  scale_x_log10(limits = c(1, 120000), labels = comma) + 
  scale_y_log10(limits = c(3e-5, 1.2)) +
  guides(color = guide_legend(override.aes = list(lwd = 1.5, alpha = 1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position = c(0.7, 0.83),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


# Inverse cumulative proportion of high-severity area
p2 = ggplot(sim_patch_ecdf, aes(x = patch, y = inv_ca,
                                color = fire_n_size, group = paste(i, fire_n_size))) + 
  geom_line(alpha = 0.5) + 
  scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.9, direction = -1,
                        name = "Number and size\nof fire events") +
  labs(x = "Patch size (ha)", y = "Inverse cumulative proportion\nof high-severity burned area") +
  scale_x_log10(limits = c(1, 120000), labels = comma) + 
  guides(color = guide_legend(override.aes = list(lwd = 1.5, alpha = 1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


plot_grid(p1, p2, nrow = 1,
          labels = c("(a)", "(b)"), label_size = 10)

ggsave("Figures/figure5_sim_patch_dist.png",
       width = 6.5, height = 3.5, units = "in", dpi = 320)






# Summary stats for table and text -------

# Create function to calculate area-weighted means
AW_mean <- function(x){sum( x * (x / sum(x)))}

patch_stats <- sim_patch_ecdf %>% 
  group_by(fire_n_size, i) %>%
  summarise(patch_mean = mean(patch),
            patch_AW_mean = AW_mean(patch),
            patch_median = median(patch),
            patch_max = max(patch)) %>%
  group_by(fire_n_size) %>%
  summarise(patch_mean_mean = mean(patch_mean),
            patch_mean_sd = sd(patch_mean),
            patch_mean_min = min(patch_mean),
            patch_mean_max = max(patch_mean),
            patch_AW_mean_mean = mean(patch_AW_mean),
            patch_AW_mean_sd = sd(patch_AW_mean),
            patch_AW_mean_min = min(patch_AW_mean),
            patch_AW_mean_max = max(patch_AW_mean),
            patch_median_mean = mean(patch_median),
            patch_median_sd = sd(patch_median),
            patch_median_min = min(patch_median),
            patch_median_max = max(patch_median),
            patch_max_mean = mean(patch_max),
            patch_max_sd = sd(patch_max),
            patch_max_min = min(patch_max),
            patch_max_max = max(patch_max))

write_csv(patch_stats, "Tables/sumstats_sim_patch.csv")


# Determine proportion of high-severity area composed of patches >= 1000 ha
sim_patch_ecdf %>%
  filter(patch >= 1000) %>%
  group_by(fire_n_size, i) %>%
  slice(which.min(patch)) %>%
  ungroup() %>% group_by(fire_n_size) %>%
  summarize(prop_mean = mean(inv_ca),
            prop_min = min(inv_ca),
            prop_max = max(inv_ca))


