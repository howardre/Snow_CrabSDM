# Libraries
library(here)
library(dplyr)
library(ggplot2)

# Load data ----
crab_summary <- readRDS(here('data', 'crab_pca.rds'))
bering_sea <- map_data("world")

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncount_mat_female = log(mature_female + 1),
                     lncount_imm_female = log(immature_female + 1),
                     lncount_leg_male = log(legal_male + 1),
                     lncount_sub_male = log(immature_male + 1),
                     pres_imm_female = ifelse(immature_female > 0, 1, 0),
                     pres_mat_female = ifelse(mature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(immature_male > 0, 1, 0),
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1))

# Map of survey data
# Legal males
ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans[crab_trans$legal_male > 0,], 
             aes(longitude, latitude,
                 size = legal_male),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 66)) +
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
       size = "Count") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncount_leg_male),
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
  geom_point(data = crab_trans[crab_trans$immature_male > 0, ], 
             aes(longitude, latitude,
                 size = immature_male),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 66)) +
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
       size = "Count") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncount_sub_male),
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
  geom_point(data = crab_trans[crab_trans$mature_female > 0, ], 
             aes(longitude, latitude,
                 size = mature_female),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 66)) +
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
       size = "Count") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncount_mat_female),
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
  geom_point(data = crab_trans[crab_trans$immature_female > 0, ], 
             aes(longitude, latitude,
                 size = immature_female),
             alpha = 0.2,
             color = "salmon") +
  scale_size_area() +
  coord_quickmap(xlim = c(-180, -156), ylim = c(54, 66)) +
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
       size = "Count") 

ggplot() +  
  geom_polygon(aes(long, lat, group = group), data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_point(data = crab_trans, 
             aes(longitude, latitude,
                 color = lncount_imm_female),
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
  geom_point(data = subset(crab_trans, !is.na(lncount_obs_leg_male)), 
             aes(longitude, latitude,
                 color = lncount_obs_leg_male),
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
  geom_point(data = subset(crab_trans, !is.na(lncount_obs_sub_male)), 
             aes(longitude, latitude,
                 color = lncount_obs_sub_male),
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
  geom_point(data = subset(crab_trans, !is.na(lncount_obs_female)), 
             aes(longitude, latitude,
                 color = lncount_obs_female),
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