library(here)
library(terra)
library(tidyverse)
library(maps)
library(mapdata)
library(akgfmaps)
library(viridis)
library(patchwork)

EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")
EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon
bering_sea <- map_data("world")

crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds'))[, c(26, 25, 6, 1, 2)]

crab_matrix_cold <- crab_summary[crab_summary$year == 2010,]
crab_matrix_warm <- crab_summary[crab_summary$year == 2019,]

# Get average RMSE per station
df_cold <- merge(crab_matrix_cold,
                 EBS_trans,
                 by.x = "station",
                 by.y = "STATIONID") # make spatial
sf_cold <- sf::st_as_sf(df_cold, # turn into sf object to plot
                        crs = st_crs(4269))
# Plot RMSE
cold_plot <- ggplot() +
  geom_polygon(aes(long, lat, group = group),
               data = bering_sea,
               fill = "lightyellow4",
               colour = "black") +
  geom_sf(data = sf_cold,
          aes(fill = temperature),
          inherit.aes = FALSE) +
  coord_sf(xlim = c(-179.5, -156),  # cannot extend out to -180, cuts off latitude labels
           ylim = c(54, 66),
           expand = FALSE) +
  scale_size_area() +
  scale_fill_viridis(option = "mako",
                     limits = c(-1.6, 15.3)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 24, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 16),
        axis.title = element_text(family = "serif", size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 20),
        legend.title = element_text(family = "serif", size = 18),
        legend.text = element_text(family = "serif", size = 16)) +
  labs(x = "Longitude \u00B0W",
       y = "Latitude \u00B0N",
       fill = "Temperature",
       title = "Cold Year") 

# Get average temperature
df_warm <- merge(crab_matrix_warm,
                 EBS_trans,
                 by.x = "station",
                 by.y = "STATIONID") # make spatial
sf_warm <- sf::st_as_sf(df_warm, # turn into sf object to plot
                        crs = st_crs(4269))
# Plot Warm
warm_plot <- ggplot() +
  geom_polygon(aes(long, lat, group = group),
               data = bering_sea,
               fill = "lightyellow4",
               colour = "black") +
  geom_sf(data = sf_warm,
          aes(fill = temperature),
          inherit.aes = FALSE) +
  coord_sf(xlim = c(-179.5, -156),  # cannot extend out to -180, cuts off latitude labels
           ylim = c(54, 66), 
           expand = FALSE) +
  scale_size_area() +
  scale_fill_viridis(option = "mako",
                     limits = c(-1.6, 15.3)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 24, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 16),
        axis.title = element_text(family = "serif", size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 20),
        legend.title = element_text(family = "serif", size = 18),
        legend.text = element_text(family = "serif", size = 16)) +
  labs(x = "Longitude \u00B0W",
       y = " ",
       fill = "Temperature",
       title = "Warm Year") 

windows(width = 18, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
cold_plot + warm_plot + plot_layout(guides = 'collect')
dev.copy(jpeg,
         here('results',
              'cold_warm_year.jpg'),
         height = 10,
         width = 18,
         res = 200,
         units = 'in')
dev.off()