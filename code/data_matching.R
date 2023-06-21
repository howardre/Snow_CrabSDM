# Title: Snow Crab Data Exploration
# Purpose: Investigate the three data sources
# Data created: 08/19/2022

# Load libraries ----
library(tidyverse)
library(here)
library(lubridate)
library(date)
library(RANN)
library(raster)
library(sp)
library(mice)
library(multilevel)
library(lattice)
library(sf)
library(akgfmaps)

# Functions ----
clean_data <- function(data){
  data$date <- as.Date(data$sampdate, "%m-%d-%Y")
  data_frame <- mutate(data,
                       year = lubridate::year(date),
                       month = lubridate::month(date),
                       day = lubridate::day(date))
  data_frame$doy <- as.numeric(mdy.date(data_frame$month, data_frame$day, 1960))
  data_frame[data_frame == -9] <- NA
  data_frame <- data_frame[!is.na(data_frame$latitude), ]
  data_frame <- data_frame[!is.na(data_frame$longitude), ]
  return(data_frame)
}

clean_offload <- function(data){
  data$date <- as.Date(data$sampdate, "%m-%d-%Y")
  data_frame <- mutate(data,
                       year = lubridate::year(date),
                       month = lubridate::month(date),
                       day = lubridate::day(date))
  data_frame$doy <- as.numeric(mdy.date(data_frame$month, data_frame$day, 1960))
  data_frame[data_frame == -9] <- NA
  return(data_frame)
}

# Convert to UTM
lon2UTM <- function(longitude){
  (floor((longitude + 180)/6) %% 60) + 1
}

# Get BCS proportions
bcs_calc <- function(data, stage){
  calc_data <- data %>%
    filter(mat_sex == stage) %>%
    group_by(station, year) %>%
    dplyr::mutate(count_sum = sum(count),
                  count_bcs = ifelse(bcs == 'N', 0, count)) %>%
    dplyr::mutate(bcs_prop = count_bcs / count_sum) %>%
    dplyr::mutate(bcs_prop = ifelse(is.nan(bcs_prop), 0, bcs_prop)) %>%
    as.data.frame
  filtered_data <- calc_data %>%
    dplyr::select(!c(count, cpue, count_bcs, bcs))
  
  merge_data <- filtered_data %>%
    group_by(station, year) %>%
    dplyr::mutate(bcs_prop = sum(bcs_prop))
  
  final_data <- merge_data[!duplicated(merge_data), ] 
  
  return(final_data)
}

combine_data <- function(data){
  mat_female <- bcs_calc(data, "Mature Female")
  imm_female <- bcs_calc(data, "Immature Female")
  leg_male <- bcs_calc(data, "Legal Male")
  sub_male <- bcs_calc(data, "Immature Male")
  mat_male <- bcs_calc(data, "Mature Male")
  
  df <- rbind(mat_female, imm_female, leg_male, sub_male, mat_male)
  
  return(df)
}

# Load data ----
# Catch of snow crab in directed fishery
crab_dump_most <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-1995-2020_crab_dump.csv'))
crab_potsum_most <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-1995-2020_potsum.csv'))
crab_retained_most <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-1995-2020_retained_size_freq.csv'))
crab_dump_21 <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-2021_crab_dump.csv'))
crab_potsum_21 <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-2021_potsum.csv'))
crab_retained_21 <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-2021_retained_size_freq.csv'))

crab_dump <- rbind(crab_dump_most, crab_dump_21)
crab_potsum <- rbind(crab_potsum_most, crab_potsum_21)
crab_retained <- rbind(crab_retained_most, crab_retained_21)

rm(crab_dump_most, crab_dump_21, crab_potsum_most, crab_potsum_21, crab_retained_most, crab_retained_21)

# Bycatch of snow crab in other crab fisheries
bycatch_dump <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2021', 'snowcrab_bycatch-1995-2021_crab_dump.csv'))
bycatch_potsum <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2021', 'snowcrab_bycatch-1995-2021_potsum.csv'))
bycatch_retained <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2021', 'snowcrab_bycatch-1995-2021_retained_size_freq.csv'))

# EBS survey data
crab_survey <- read_csv(here('data/Snow_CrabData', 'station_cpue_BCS_snowEBSNBS.csv'), col_select = -c(1))
crab_env <- read_csv(here('data/Snow_CrabData', 'SAPsurvey_tempdepth.csv'), col_select = -c(1))
crab_julian <- read_csv(here('data/Snow_CrabData', 'year_station_julian_day.csv'))
names(crab_survey) <- tolower(names(crab_survey))
names(crab_env) <- tolower(names(crab_env))
crab_survey <- as.data.frame(crab_survey)

# Add the julian days
crab_env$date <- as.Date(with(crab_env, paste(survey_year, survey_month, survey_day, sep = "-")),
                         "%Y-%m-%d")
crab_env$julian <- lubridate::yday(crab_env$date)
crab_all <- crab_survey %>% left_join(crab_env,
                                       by = c("akfin_survey_year" = "survey_year",
                                              "gis_station" = "gis_station"))
crab_all <- crab_all[-c(13:16)]

crab_fix <- crab_all %>% left_join(crab_julian,
                                   by = c("akfin_survey_year" = "year",
                                          "gis_station" = "station"))
crab_fill <- crab_fix %>% 
  mutate(julian = coalesce(julian.x, julian.y)) %>%
  dplyr::select(-c(julian.x, julian.y))

crab_reduced <- crab_fill %>%
  filter(akfin_survey_year >= 1995) %>% # observer data only available until this date
  rename(year = akfin_survey_year,
         station = gis_station,
         latitude = mid_latitude.x,
         longitude = mid_longitude.x,
         depth = bottom_depth,
         temperature = gear_temperature)

# Calculate proportion of crab with BCS and get total count
bcs_calc <- function(data, stage){
  # calculate the proportion of crab with BCS in a haul for sex/stage
  calc_data <- data %>%
    filter(mat_sex == stage) %>%
    group_by(station, year) %>%
    dplyr::mutate(count_sum = sum(count),
                  count_bcs = ifelse(bcs == 'N', 0, count)) %>%
    dplyr::mutate(bcs_prop = count_bcs / count_sum) %>%
    dplyr::mutate(bcs_prop = ifelse(is.nan(bcs_prop), 0, bcs_prop)) %>%
    as.data.frame
  
  # Remove unnecessary columns
  filtered_data <- calc_data %>%
    dplyr::select(!c(count, cpue, count_bcs, bcs))
  
  # Apply BCS proportion to each station/year for both duplicate rows
  merge_data <- filtered_data %>%
    group_by(station, year) %>%
    dplyr::mutate(bcs_prop = sum(bcs_prop))
  
  # Get rid of separate Y and N columns (duplicates at this point)
  final_data <- merge_data[!duplicated(merge_data), ] # removes the 
  
  return(final_data)
}

combine_data <- function(data){
  # Do the BCS calculation for each sex/stage
  mat_female <- bcs_calc(data, "Mature Female")
  imm_female <- bcs_calc(data, "Immature Female")
  leg_male <- bcs_calc(data, "Legal Male")
  sub_male <- bcs_calc(data, "Immature Male")
  mat_male <- bcs_calc(data, "Mature Male")
  
  # Put the data back together
  # This should be smaller than the input data because the Y/N rows are now gone
  df <- rbind(mat_female, imm_female, leg_male, sub_male, mat_male)
  
  return(df)
}

# Test the bcs_calc function
# mat_female <- bcs_calc(crab_reduced, "Mature Female")
# leg_male <- bcs_calc(crab_reduced, "Legal Male")

# Calculate the proportion of crab with BCS at a given station/year combo
crab_bcs <- combine_data(crab_reduced)

# Pivot to wide format
survey_wide <- crab_bcs %>%
  pivot_wider(names_from = mat_sex, 
              values_from = c(count_sum, bcs_prop))
survey_wide <- as.data.frame(survey_wide)
survey_jdate <- survey_wide %>% mutate(jdate = as.Date(strptime(paste(as.numeric(survey_wide$year),
                                                                      as.numeric(survey_wide$julian)),
                                                                format = "%Y%j"))) %>%
  mutate(date = coalesce(date, jdate))
survey_rename <- survey_jdate %>% 
  rename(sublegal_male = 'count_sum_Immature Male',
         mature_male = 'count_sum_Mature Male',
         legal_male = 'count_sum_Legal Male',
         immature_female = 'count_sum_Immature Female',
         mature_female = 'count_sum_Mature Female',
         bcs_sublegal_male = 'bcs_prop_Immature Male',
         bcs_mature_male = 'bcs_prop_Mature Male',
         bcs_legal_male = 'bcs_prop_Legal Male',
         bcs_immature_female = 'bcs_prop_Immature Female',
         bcs_mature_female = 'bcs_prop_Mature Female') %>%
  mutate(sublegal_male = ifelse(is.na(sublegal_male), 0, sublegal_male), # fill NA values with 0a
         mature_male = ifelse(is.na(mature_male), 0, mature_male),
         legal_male = ifelse(is.na(legal_male), 0, legal_male),
         immature_female = ifelse(is.na(immature_female), 0, immature_female),
         mature_female = ifelse(is.na(mature_female), 0, mature_female),
         bcs_sublegal_male = ifelse(is.na(bcs_sublegal_male), 0, bcs_sublegal_male), # fill NA values with 0a
         bcs_mature_male = ifelse(is.na(bcs_mature_male), 0, bcs_mature_male),
         bcs_legal_male = ifelse(is.na(bcs_legal_male), 0, bcs_legal_male),
         bcs_immature_female = ifelse(is.na(bcs_immature_female), 0, bcs_immature_female),
         bcs_mature_female = ifelse(is.na(bcs_mature_female), 0, bcs_mature_female))
survey_rename <- survey_rename[-23]

# Reformat data (summarize males, split up dates)
crab_detailed <- clean_data(crab_dump)
crab_summary <- clean_data(crab_potsum)
crab_summary$total <- crab_summary$female + crab_summary$tot_legal + crab_summary$sublegal
crab_summary$male <- crab_summary$tot_legal + crab_summary$sublegal
crab_offload <- clean_offload(crab_retained)

bycatch_detailed <- clean_data(bycatch_dump)
bycatch_summary <- clean_data(bycatch_potsum)
bycatch_summary$total <- bycatch_summary$female + bycatch_summary$tot_legal + bycatch_summary$sublegal
bycatch_summary$male <- bycatch_summary$tot_legal + bycatch_summary$sublegal
bycatch_offload <- clean_offload(bycatch_retained)

# Need to match up observer data with the survey data or make PCA
# Group observer data for each season (Nov - Mar)
# Match closest observer data from preceding season to a station for a year
# Get the survey grid
EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
class(EBS)

# Plot the survey grid
ggplot() +
  geom_sf(data = EBS$survey.grid) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks) + 
  theme_bw()

# Convert to lat/lon, match catches to station polygons
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")

EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon

crab_sf <- st_as_sf(crab_summary, 
                    coords = c("longitude", "latitude"), 
                    crs = 4269)
bycatch_sf <- st_as_sf(bycatch_summary, 
                       coords = c("longitude", "latitude"), 
                       crs = 4269)
survey_sf <- st_as_sf(survey_wide,
                      coords = c("longitude", "latitude"), 
                      crs = 4269)

sf_use_s2(FALSE)
observer_df <- st_join(crab_sf, EBS_trans, left = FALSE)
bycatch_df <- st_join(bycatch_sf, EBS_trans, left = FALSE)

ggplot() +
  geom_sf(data = EBS$survey.grid) +
  geom_sf(data = observer_df,
          aes(color = STATIONID)) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks)  # shows that they are grouped into the polygons

# Check for zeros
observer_df$legal_presence <- ifelse(observer_df$tot_legal > 0, 1, 0)
observer_df$sublegal_presence <- ifelse(observer_df$sublegal > 0, 1, 0)
observer_df$female_presence <- ifelse(observer_df$female > 0, 1, 0)

ggplot() +
  geom_sf(data = EBS$survey.grid) +
  geom_sf(data = observer_df,
          aes(color = factor(legal_presence)),
          size = 1.5) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Legal Male Presence",
       color = "Presence") +
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
        legend.text = element_text(family = "serif", size = 14))

ggplot() +
  geom_sf(data = EBS$survey.grid) +
  geom_sf(data = observer_df,
          aes(color = factor(sublegal_presence)),
          size = 1.5) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Sublegal Male Presence",
       color = "Presence") +
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
        legend.text = element_text(family = "serif", size = 14))

ggplot() +
  geom_sf(data = EBS$survey.grid) +
  geom_sf(data = observer_df,
          aes(color = factor(female_presence)),
          size = 1.5) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Female Presence",
       color = "Presence") +
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
        legend.text = element_text(family = "serif", size = 14))

# Combine the observer data
# observer_combined <- rbind(observer_df, bycatch_df)

# Average the female and male CPUE by season and nearest station
# Need to include previous year when grouping by season (months 11 and 12)
observer_dates <- mutate(observer_df,
                         date_lag = ymd(date) %m+% - months(-2),
                         year_lag = year(date_lag)) # lag to group by season
observer_summarized <- observer_dates %>%
  group_by(year_lag, STATIONID) %>%
  summarise(obs_male_legal = mean(tot_legal),
            obs_male_sub = mean(sublegal),
            obs_female = mean(female))

saveRDS(observer_dates, file = here('data/Snow_CrabData', 'observer_all.rds'))
saveRDS(observer_summarized, file = here('data/Snow_CrabData', 'observer_summarized.rds'))
#saveRDS(bycatch_summarized, file = here('data/Snow_CrabData', 'bycatch_df.rds'))

# Match to the survey data
survey_combined <- merge(survey_rename, observer_summarized,
                         by.x = c("year", "station"),
                         by.y = c("year_lag", "STATIONID"),
                         all.x = T)[-26] # only 1455 rows with observer data

# Add environmental data ----
# Load data
sst_data <- as.data.frame(read_csv(here('data', 'sst_latlon.csv'), col_select = -c(1)))
sst_data$month <- match(sst_data$month, month.abb)
ice_data <- as.data.frame(read_csv(here('data', 'ice_latlon.csv'), col_select = -c(1)))
ice_data$month <- match(ice_data$month, month.abb)

ice_data_filtered <- dplyr::filter(ice_data, month %in% c(1:4))

ice_sf <- st_as_sf(ice_data_filtered,
                   coords = c("lon", "lat"), 
                   crs = 4269)

ice_df <- ice_sf %>% st_join(EBS_trans, join = st_intersects, left = FALSE)

ice_means <- ice_df %>%
  group_by(year, STATIONID) %>%
  summarise(ice_early = mean(ice[month %in% c(1, 2)]),
            ice_late = mean(ice[month %in% c(3, 4)]))

ice_final <- ice_means %>%
  group_by(year, STATIONID) %>%
  summarise(ice_mean = mean(c(ice_early, ice_late)))

ggplot() +
  geom_sf(data = EBS$survey.grid) +
  geom_sf(data = ice_final,
          aes(color = factor(ice_mean)),
          size = 1.5) +
  coord_sf(xlim = EBS$plot.boundary$x,
           ylim = EBS$plot.boundary$y) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks)

phi_raster <- raster(here('data', 'EBS_phi_1km.gri'))
phi_pts <- rasterToPoints(phi_raster, spatial = T)
proj4string(phi_pts) # warning here related to the retirement of rgdal package
phi_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
phi_trans <- spTransform(phi_pts, CRS(phi_prj))
proj4string(phi_trans)
phi_data <- as.data.frame(phi_trans)

# Needs to be matched to the survey data 
# Currently matching the SST and ice to the time period of the survey
z <- lon2UTM(survey_combined[1, "longitude"])
survey_combined <- tibble::rowid_to_column(survey_combined, "index")
survey_data <- survey_combined

# Add index to summary to match back lat, lon
coordinates(survey_data) <- c("longitude", "latitude")
proj4string(survey_data) <- CRS("+proj=longlat +datum=WGS84")

data_xy <- spTransform(survey_data, CRS(paste0("+proj=utm +zone=", z, " ellps=WGS84")))
data_xy <- as.data.frame(data_xy)

coordinates(phi_data) <- c("x", "y")
proj4string(phi_data) <- CRS("+proj=longlat +datum=WGS84")

phi_data_xy <- spTransform(phi_data, CRS(paste0("+proj=utm +zone=", z, " ellps=WGS84")))
phi_data_xy <- as.data.frame(phi_data_xy)

# coordinates(sst_data) <- c("lon", "lat")
# proj4string(sst_data) <- CRS("+proj=longlat +datum=WGS84")

# sst_data_xy <- spTransform(sst_data, CRS(paste0("+proj=utm +zone=", z, " ellps=WGS84")))
# sst_data_xy <- as.data.frame(sst_data_xy)

# Everything must be a data frame, tables do not work
data_xy[, c(27, 28)] <- as.data.frame(RANN::nn2(phi_data_xy[, c('y', 'x')],
                                                data_xy[, c('latitude', 'longitude')],
                                                k = 1))
data_xy$phi <- phi_data_xy[c(data_xy$nn.idx), 1] # Match nearest phi value
data_xy <- data_xy[-c(27, 28)]

# Remove SST match for now, not using and must eliminate many rows
# data_xy <- data_xy[!is.na(data_xy$month), ]
# data_xy[, c(23, 24)] <- as.data.frame(RANN::nn2(sst_data_xy[, c('lat', 'lon', 'month', 'year')],
#                                                 data_xy[, c('latitude', 'longitude', 'month', 'year')],
#                                                 k = 1))
# data_xy$sst <- sst_data_xy[c(data_xy$nn.idx), 2] # Match nearest temperature value
# data_xy <- data_xy[-c(23, 24)]

# Match ice to survey stations
crab_ice <- merge(data_xy, ice_final,
                  by.x = c("year", "station"),
                  by.y = c("year", "STATIONID"),
                  all.x = T)

# Convert back to lat, lon
crab_final <- crab_ice[-c(6, 7, 25, 26, 29)]
crab_final$latitude <- survey_combined$latitude[match(survey_combined$index, crab_final$index)]
crab_final$longitude <- survey_combined$longitude[match(survey_combined$index, crab_final$index)]

pcod_catch <- read_csv(here('data/Snow_CrabData', 'GAP_survey_pcod.csv'), col_select = -c(1))
names(pcod_catch) <- tolower(names(pcod_catch))

survey_pcod <- merge(crab_final, pcod_catch,
                     by.x = c("year", "station"),
                     by.y = c("year", "stationid"),
                     all.x = T)[-c(27:35, 38)]

names(survey_pcod)[names(survey_pcod) == "cpue_kgha"] <- "pcod_cpue"
names(survey_pcod)[names(survey_pcod) == "latitude.x"] <- "latitude"
names(survey_pcod)[names(survey_pcod) == "longitude.x"] <- "longitude"

# Match ROMS bottom temperature
temp_data <- readRDS(here('data', 'hindcast_temps.rds'))

temp_sf <- st_as_sf(temp_data,
                    coords = c("lon", "lat"), 
                    crs = 4269)

temp_df <- temp_sf %>% st_join(EBS_trans, join = st_intersects, left = FALSE) # match points to stations

temp_df$day <- lubridate::day(temp_df$date)
temp_df$julian <- as.numeric(mdy.date(temp_df$month, temp_df$day, 1960))

temp_means <- temp_df %>%
  group_by(year, STATIONID) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE)) # calculate mean temperature based on station

crab_temp <- merge(survey_pcod, temp_means,
                   by.x = c("year", "station"),
                   by.y = c("year", "STATIONID"),
                   all.x = TRUE)

# survey_temp <- crab_temp %>% mutate(temperature = coalesce(temperature, mean_temp)) 

# Interpolate depths
# Use other data to fill in missing depth values
depth_loess <- loess(depth ~ longitude * latitude,
                     span = 0.01,
                     degree = 2,
                     data = crab_temp)
crab_temp$depth_pred <- predict(depth_loess, newdata = crab_temp)

survey_depth <- crab_temp %>% mutate(depth = coalesce(depth, depth_pred))

crab_summary <- survey_depth %>% drop_na(temperature, depth, phi, julian)

# Save final data set
saveRDS(crab_summary, file = here('data/Snow_CrabData', 'crab_summary.rds'))
