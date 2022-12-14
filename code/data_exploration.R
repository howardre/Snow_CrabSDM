# Libraries
library(here)
library(dplyr)
library(ggplot2)

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds'))
bering_sea <- map_data("world")

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncpue_mat_female = log(mature_female + 1),
                     lncpue_imm_female = log(immature_female + 1),
                     lncpue_leg_male = log(legal_male + 1),
                     lncpue_sub_male = log(immature_male + 1),
                     lncpue_obs_female = log(obs_female + 1),
                     lncpue_obs_sub_male = log(obs_male_sub + 1),
                     lncpue_obs_leg_male = log(obs_male_legal + 1))

# Map of survey data
# Legal males
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 size = lncpue_leg_male),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Legal Males",
       x = "Longitude",
       y = "Latitude",
       size = "ln(catch+1)") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncpue_leg_male),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Legal Males",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Sublegal Males
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 size = lncpue_sub_male),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Sublegal Males",
       x = "Longitude",
       y = "Latitude",
       size = "ln(catch+1)") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncpue_sub_male),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Sublegal Males",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Mature Females
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 size = lncpue_mat_female),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Mature Females",
       x = "Longitude",
       y = "Latitude",
       size = "ln(catch+1)") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncpue_mat_female),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Mature Females",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Immature Females
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 size = lncpue_imm_female),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Immature Females",
       x = "Longitude",
       y = "Latitude",
       size = "ln(catch+1)") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncpue_imm_female),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Immature Females",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Observer Legal Males
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = subset(crab_trans, !is.na(lncpue_obs_leg_male)), 
             aes(longitude, latitude,
                 color = lncpue_obs_leg_male),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Legal Males (Observer)",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Observer Sublegal Males
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = subset(crab_trans, !is.na(lncpue_obs_sub_male)), 
             aes(longitude, latitude,
                 color = lncpue_obs_sub_male),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Sublegal Males (Observer)",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# Observer Females
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = subset(crab_trans, !is.na(lncpue_obs_female)), 
             aes(longitude, latitude,
                 color = lncpue_obs_female),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Females (Observer)",
       x = "Longitude",
       y = "Latitude",
       color = "ln(catch+1)") +
  facet_wrap(~ year)

# PCA
# Group values by bins
crab_grouped <- crab_trans %>%
  group_by(station) %>%
  summarise(legal_male_loading = mean(legal_male_loading),
            sublegal_male_loading = mean(sublegal_male_loading),
            female_loading = mean(female_loading),
            latitude = mean(latitude),
            longitude = mean(longitude))

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = subset(crab_grouped, !is.na(female_loading)), 
             aes(longitude, latitude,
                 color = female_loading),
             alpha = 0.2,
             size = 3) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 62)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(title = "Legal Male Loadings",
       x = "Longitude",
       y = "Latitude",
       color = "loading values")