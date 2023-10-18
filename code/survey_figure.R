# Libraries
library(ggplot2)
library(sf)
library(akgfmaps)

# Prepare data
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds')) %>%
  dplyr::select(-geometry)

crab_trans <- mutate(crab_summary,
                     lncount_mat_female = log(mature_female + 1),
                     lncount_imm_female = log(immature_female + 1),
                     lncount_leg_male = log(legal_male + 1),
                     lncount_sub_male = log(sublegal_male + 1),
                     pres_imm_female = ifelse(immature_female > 0, 1, 0),
                     pres_mat_female = ifelse(mature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(sublegal_male > 0, 1, 0),
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean),
         year_f != 2022,
         year > 2010) 

legal_males <- crab_trans %>%
  dplyr::select(longitude, latitude, 
                lncount_leg_male, legal_male, 
                pres_leg_male, year, station) %>%
  tidyr::drop_na(lncount_leg_male)

EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")
EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon
bering_sea <- map_data("world")

legal_male_df <- merge(legal_males,
                       EBS_trans,
                       by.x = "station",
                       by.y =  "STATIONID")
legal_male_sf <- st_as_sf(legal_male_df, 
                          crs = st_crs(4269))

# Make map
ggplot() +
  geom_polygon(aes(long, lat, group = group), 
               data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  geom_sf(data = legal_male_sf) +
  geom_jitter(data = legal_male_sf,
             aes(geometry = geometry,
                 size = lncount_leg_male),
             color = "royalblue3",
             alpha = 0.2,
             stat = "sf_coordinates",
             width = 0.05,
             height = 0.05) +
  coord_sf(xlim = c(-179.5, -156),
           ylim = c(54, 66), 
           expand = FALSE) +
  scale_size_area() +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 16),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 14)) +
  labs(x = "Longitude \u00B0W",
       y = "Latitude \u00B0N",
       size = "ln(count + 1)",
       title = "Legal Male Crab Survey Catches") 