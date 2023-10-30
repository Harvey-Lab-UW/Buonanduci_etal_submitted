# Study area figure

require(tidyverse)
require(sf)
require(raster)
require(ggsn)
require(cowplot)

# Load data layers -------

# All NW fires
fire_perims <- st_read("Data/Spatial/fire_perims.gpkg") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 500) %>%
  mutate(Fire_Regime = factor(Fire_Regime, levels = c("High", "Mixed", "Low")))

# Separate into NWC and non-NWC fires
NWC_perims <- fire_perims %>% filter(NW_Cascadia == 1)
NW_perims <- fire_perims %>% filter(NW_Cascadia == 0)

# US states
states <- st_read("Data/Spatial/states.gpkg")
wa_or_line <- st_read("Data/Spatial/WA_OR_line.gpkg")

# NW region 
NW_region <- st_read("Data/Spatial/region_NW.gpkg") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 1000)

# NWC region
NWC_region <- st_read("Data/Spatial/region_NW_Cascadia.gpkg") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 1000)

# Cascade Crest (WA and OR)
CCrest <- st_read("Data/Spatial/CCrest.gpkg") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 1000)



# Plotting entire study region -----------
plt_extent1 <- extent(fire_perims)

p1 = ggplot() +
  xlim( c(plt_extent1[1] - 100000, plt_extent1[2]) ) +
  ylim( c(plt_extent1[3], plt_extent1[4]) ) +
  labs( x = NULL, y = NULL ) +
  geom_sf(data = states, fill = "gray90", color = NA) +
  geom_sf(data = NW_region, fill = "gray", color = NA) +
  geom_sf(data = NWC_region, fill = "blue", color = NA, alpha = 0.2) +
  geom_sf(data = states, fill = NA, color = "white") +
  geom_sf(data = NW_perims, mapping = aes(fill = Fire_Regime), color = NA) +
  geom_sf(data = NWC_perims, fill = "blue", color = NA) +
  scale_fill_viridis_d(option = "inferno",
                       name = "Historical fire regime",
                       labels = c("High-severity",
                                  "Mixed-severity",
                                  "Low-severity")) +
  annotate("text", 
           x = plt_extent1[1] - 100000, 
           y = plt_extent1[3] + 0.8*(plt_extent1[4] - plt_extent1[3]), 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Northwestern\nCascadia") +
  ggsn::scalebar(data = NW_perims, dist = 100, dist_unit = "km", 
                 transform = FALSE, location = "bottomright", st.size = 2,
                 st.bottom = FALSE, st.dist = 0.03, border.size = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.13),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, 'mm'),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



# Plotting NW Cascadia study region -----------
plt_extent2 <- extent(NWC_region)

p2 = ggplot() +
  xlim( c(plt_extent2[1], plt_extent2[2] + 5000) ) +
  ylim( c(plt_extent2[3] - 20000, plt_extent2[4]) ) +
  labs( x = NULL, y = NULL ) +
  geom_sf(data = states, fill = "gray90", color = NA) +
  geom_sf(data = NW_region, fill = "gray", color = NA) +
  geom_sf(data = NWC_region, fill = "blue", color = NA, alpha = 0.2) +
  geom_sf(data = wa_or_line, fill = NA, color = "white") +
  geom_sf(data = NWC_perims, fill = "blue", color = NA) +
  geom_sf(data = CCrest, color = "black") +
  annotate("text", 
           x = plt_extent2[1] + 0.73*(plt_extent2[2] - plt_extent2[1]), 
           y = plt_extent2[3] + 0.4*(plt_extent2[4] - plt_extent2[3]), 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Cascade\nCrest") +
  ggsn::scalebar(data = NWC_perims, dist = 50, dist_unit = "km", 
                 transform = FALSE, location = "bottomright", st.size = 2,
                 st.bottom = FALSE, st.dist = 0.03, border.size = 0.5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# Calculate relative widths
width1 <- plt_extent1[2] - plt_extent1[1] + 100000
height1 <- plt_extent1[4] - plt_extent1[3]
width2 <- (plt_extent2[2] + 5000) - plt_extent2[1]
height2 <- plt_extent2[4] - (plt_extent2[3] - 20000)

plot_grid(p1, p2, nrow = 1, rel_widths = c(width1/height1, width2/height2),
          labels = c("(a)", "(b)"), label_size = 10,
          vjust = 3, hjust = c(-2, -1.5))


ggsave(paste0("Figures/figure1_study_region.png"), bg = "white",
       width = 6.5, height = 3, units = "in", dpi = 320)

