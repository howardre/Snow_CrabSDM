# Libraries
library(here)
library(dplyr)
library(ggplot2)

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))
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
                 size = lncpue_leg_male),
             alpha = 0.2,
             color = "salmon") +
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
       size = "ln(catch+1)") +
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