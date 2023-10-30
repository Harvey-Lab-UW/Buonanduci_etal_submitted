# Scaling analysis
# Simulate core area and distance-to-seed distributions

# Load packages
require(tidyverse)
require(quantreg)
require(quantregGrowth)
require(mgcv)
require(gratia)
require(scales) # for numeric axis labels with commas
require(cowplot)
require(ggimage) # for embedding images in plots


# Load data --------------

# Landscape metrics data file
fire_metrics <- read_csv("Data/fire_metrics.csv") %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("Low", "Mixed", "High"))) 

# Pull out NW and WC fires
fire_metrics_NW <- filter(fire_metrics, NW_Cascadia == 0) 
fire_metrics_NWC <- filter(fire_metrics, NW_Cascadia == 1) 

# Pull out high severity FRGs
fire_metrics_high <- filter(fire_metrics_NW, Fire_Regime == "High")

# Combine NW high-severity regime fires with NW Cascadia fires
fire_metrics_NWC_high <- bind_rows(fire_metrics_NWC, fire_metrics_high)



# Total core area -----

# Fit quantile regression model to core area
mod_core <- gcrq(log_total_core ~ ps(log_fire_area, monotone = 1),
                 lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)

# Predict at fire sizes 1000, 10000, 100000
pred_core <- predict.gcrq(mod_core, 
                          newdata = data.frame("log_fire_area" = c(3, 4, 5)))


# Simulate total core area under three scenarios:
# 10 fires, size 1e5 ha
# 100 fires, size 1e4 ha
# 1000 fires, size 1e3 ha
set.seed(321)

sim_core5 <- rep(NA, 100)
for(i in 1:100){
  fires <- sample(pred_core[3,], 10, replace = TRUE)
  fires <- 10^fires
  sim_core5[i] <- log10(sum(fires))
}
sim_core4 <- rep(NA, 100)
for(i in 1:100){
  fires <- sample(pred_core[2,], 100, replace = TRUE)
  fires <- 10^fires
  sim_core4[i] <- log10(sum(fires))
}
sim_core3 <- rep(NA, 100)
for(i in 1:100){
  fires <- sample(pred_core[1,], 1000, replace = TRUE)
  fires <- 10^fires
  sim_core3[i] <- log10(sum(fires))
}

sim_cores <- tibble(core = c(sim_core3, sim_core4, sim_core5),
                    fire_n_size = c(rep("1,000 x 1,000 ha", 100),
                                    rep("100 x 10,000 ha", 100),
                                    rep("10 x 100,000 ha", 100))) %>%
  mutate(core_ha = 10^core) %>%
  mutate(fire_n_size = factor(fire_n_size, levels = c("1,000 x 1,000 ha", 
                                                      "100 x 10,000 ha", 
                                                      "10 x 100,000 ha"))) 


p1 = ggplot(sim_cores, aes(y = core_ha, x = fire_n_size, color = fire_n_size)) +
  geom_jitter(alpha = 0.3, size = 0.8) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.9) +
  scale_y_continuous(name = "Total core area (ha)\n(Area exceeding DTS of 150 m)",
                     labels = comma, limits = c(50000, 200000)) +
  scale_x_discrete(name = "Number and size of fires",
                   labels = c("1,000 x\n1,000 ha", 
                              "100 x\n10,000 ha", 
                              "10 x\n100,000 ha")) +
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



# Distance-to-seed distributions -----

fire_metrics_NWC_high <- fire_metrics_NWC_high %>% drop_na(log_SDC)

# Fit quantile regression model to SDC parameter
mod_SDC <- gcrq(log_SDC ~ ps(log_fire_area, monotone = -1),
                lambda0 = 5, tau = seq(0.01,0.99,0.01), data = fire_metrics_NWC_high)

# Predict at fire sizes 1000, 10000, 100000
pred_SDC_3 <- predict.gcrq(mod_SDC, newdata = data.frame("log_fire_area" = 3))
pred_SDC_4 <- predict.gcrq(mod_SDC, newdata = data.frame("log_fire_area" = 4))
pred_SDC_5 <- predict.gcrq(mod_SDC, newdata = data.frame("log_fire_area" = 5))

# Fit gam model to high-severity & forested proportion (as function of fire size & SDC)
mod_HS_forest_prop <- gam(HS_forest_prop ~ s(log_fire_area) + s(log_SDC), 
                          family = betar(link="logit"), data = fire_metrics_NWC_high)

# Simulate distance-to-seed distributions under three scenarios:
# 10 fires, size 1e5 ha
# 100 fires, size 1e4 ha
# 1000 fires, size 1e3 ha

# Step 1: Simulate log_SDC, randomly selected from quantile regression
# Step 2: Simulate high-severity forest proportion from log_SDC and fire size
# Step 3: Use log_SDC and HS_forest_prop to generate distance-to-seed distribution
#         from each hypothetical fire

# Distances at which to calculate areas
# FWIW, maximum distance-to-seed observation in dataset is 2700 m
dist <- seq(0, 1500, by=30)

# Function to simulate distance-to-seed distributions
sim_DTS_dist <- function(log_fire_size = 3, pred_SDC, mod_HS_forest_prop){
  
  # Simulate SDC
  sim_log_SDC <- pred_SDC[round(runif(1, 1, 99), 0)]
  sim_SDC <- 10^sim_log_SDC
  
  # Simulate high-severity & forested proportion
  sim_HS_forest_prop <- simulate(mod_HS_forest_prop,
                                 newdata = data.frame(log_fire_area = log_fire_size,
                                                      log_SDC = sim_log_SDC)) %>% as.numeric()
  
  # Calculate high-severity & forested area
  sim_HS_forest_area <- sim_HS_forest_prop * 10^log_fire_size
  
  # Calculate forested & high-severity area
  # as inverse cumulative proportions, using modified logistic function
  P <- 1 / (10^(sim_SDC*dist))
  
  # Convert inverse cumulative proportions to inverse cumulative areas
  A <- round(P * sim_HS_forest_area)
  
  return(A)
}

set.seed(123)

# 10 fires, size 1e5 ha
sim_DTS5 <- tibble(dist = dist)

for (i in 1:100){
  A <- numeric(length(dist))
  for (f in 1:10){
    A <- A + sim_DTS_dist(log_fire_size = 5, 
                          pred_SDC = pred_SDC_5, mod_HS_forest_prop = mod_HS_forest_prop)
  }
  sim_DTS5[[paste0("i", i)]] <- A
}

sim_DTS5 <- pivot_longer(sim_DTS5, cols = starts_with("i"),
                         names_to = "i", values_to = "area") %>%
  mutate(i = paste0("5", i)) %>%
  mutate(fire_n_size = "10 x 100,000 ha")


# 100 fires, size 1e4 ha
sim_DTS4 <- tibble(dist = dist)

for (i in 1:100){
  A <- numeric(length(dist))
  for (f in 1:100){
    A <- A + sim_DTS_dist(log_fire_size = 4, 
                          pred_SDC = pred_SDC_4, mod_HS_forest_prop = mod_HS_forest_prop)
  }
  sim_DTS4[[paste0("i", i)]] <- A
}

sim_DTS4 <- pivot_longer(sim_DTS4, cols = starts_with("i"),
                         names_to = "i", values_to = "area") %>%
  mutate(i = paste0("4", i)) %>%
  mutate(fire_n_size = "100 x 10,000 ha")


# 1000 fires, size 1e3 ha
sim_DTS3 <- tibble(dist = dist)

for (i in 1:100){
  A <- numeric(length(dist))
  for (f in 1:1000){
    A <- A + sim_DTS_dist(log_fire_size = 3, 
                          pred_SDC = pred_SDC_3, mod_HS_forest_prop = mod_HS_forest_prop)
  }
  sim_DTS3[[paste0("i", i)]] <- A
}

sim_DTS3 <- pivot_longer(sim_DTS3, cols = starts_with("i"),
                         names_to = "i", values_to = "area") %>%
  mutate(i = paste0("3", i)) %>%
  mutate(fire_n_size = "1,000 x 1,000 ha")



# Combined figure
sim_DTS <- bind_rows(sim_DTS5, sim_DTS4, sim_DTS3) %>%
  mutate(fire_n_size = factor(fire_n_size, levels = c("1,000 x 1,000 ha", 
                                                      "100 x 10,000 ha", 
                                                      "10 x 100,000 ha"))) 


p2 = ggplot(sim_DTS, aes(x = dist, y = area, group = i, color = fire_n_size)) +
  geom_line(alpha = 0.1) +
  scale_y_continuous(name = "Area exceeding DTS (ha)",
                     labels = comma) +
  scale_x_continuous(name = "Distance to seed (DTS; m)",
                     limits = c(0,1010)) +
  scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.9,
                        name = "Number and size of fires") +
  guides(color = guide_legend(reverse = TRUE, override.aes = list(lwd = 1.5, alpha = 1))) +
  theme_bw() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position = c(0.6, 0.7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


plot_grid(p1, p2, nrow = 1, rel_widths = c(0.4, 0.6),
          labels = c("(a)", "(b)"), label_size = 10)


ggsave("Figures/figure4_sim_core_DTS.png",
       width = 6.5, height = 3, units = "in", dpi = 320)





# Summary stats for table and text --------------

core_stats <- sim_cores %>% 
  group_by(fire_n_size) %>%
  summarize(core_mean = mean(core_ha), 
            core_sd = sd(core_ha),
            core_min = min(core_ha),
            core_max = max(core_ha))

write_csv(core_stats, "Tables/sumstats_sim_core.csv")

summary(aov(core_ha ~ fire_n_size, data = sim_cores))
TukeyHSD(aov(core_ha ~ fire_n_size, data = sim_cores))


DTS_stats <- sim_DTS %>%
  group_by(dist, fire_n_size) %>%
  summarize(core_mean = mean(area), 
            core_sd = sd(area),
            core_min = min(area),
            core_max = max(area)) %>%
  filter(dist == 150 | dist == 450 | dist == 750)

write_csv(DTS_stats, "Tables/sumstats_sim_DTS.csv")

