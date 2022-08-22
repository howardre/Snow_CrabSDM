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

# Match data for pot level data
match_data <- function(data, survey, phi_data, ice_data, sst_data){
  z <- lon2UTM(data[1, "longitude"])
  
  coordinates(survey) <- c("mid_longitude", "mid_latitude")
  proj4string(survey) <- CRS("+proj=longlat +datum=WGS84")
  
  survey_xy <- spTransform(survey, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
  survey_xy <- as.data.frame(survey_xy)

  # Add index to summary to match back lat, lon
  data <- tibble::rowid_to_column(data, "index")
  coordinates(data) <- c("longitude", "latitude")
  proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
  
  data_xy <- spTransform(data, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
  data_xy <- as.data.frame(data_xy)
  
  coordinates(phi_data) <- c("x", "y")
  proj4string(phi_data) <- CRS("+proj=longlat +datum=WGS84")
  
  phi_data_xy <- spTransform(phi_data, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
  phi_data_xy <- as.data.frame(phi_data_xy)
  
  coordinates(ice_data) <- c("lon", "lat")
  proj4string(ice_data) <- CRS("+proj=longlat +datum=WGS84")
  
  ice_data_xy <- spTransform(ice_data, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
  ice_data_xy <- as.data.frame(ice_data_xy)
  
  coordinates(sst_data) <- c("lon", "lat")
  proj4string(sst_data) <- CRS("+proj=longlat +datum=WGS84")
  
  sst_data_xy <- spTransform(sst_data, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
  sst_data_xy <- as.data.frame(sst_data_xy)
  
  # Everything must be a data frame, tables do not work
  data_xy[, c(30, 31)] <- as.data.frame(RANN::nn2(phi_data_xy[, c('y', 'x')],
                                                  data_xy[, c('latitude', 'longitude')],
                                                  k = 1))
  data_xy$phi <- phi_data_xy[c(data_xy$nn.idx), 1] # Match nearest phi value
  data_xy <- data_xy[-c(30, 31)]
  
  data_xy[, c(31, 32)] <- as.data.frame(RANN::nn2(ice_data_xy[, c('lat', 'lon', 'month', 'year')],
                                                  data_xy[, c('latitude', 'longitude', 'month', 'year')],
                                                  k = 1))
  data_xy$ice <- ice_data_xy[c(data_xy$nn.idx), 2] # Match nearest ice value
  data_xy <- data_xy[-c(31, 32)]
  
  data_xy[, c(32, 33)] <- as.data.frame(RANN::nn2(sst_data_xy[, c('lat', 'lon', 'month', 'year')],
                                                  data_xy[, c('latitude', 'longitude', 'month', 'year')],
                                                  k = 1))
  data_xy$sst <- sst_data_xy[c(data_xy$nn.idx), 2] # Match nearest temperature value
  data_xy <- data_xy[-c(32, 33)]
  
  data_xy[, c(33, 34)] <- as.data.frame(RANN::nn2(survey_wide[, c('mid_latitude', 
                                                                  'mid_longitude', 
                                                                  'akfin_survey_year')],
                                                  data_xy[, c('latitude', 
                                                              'longitude', 
                                                              'year')],
                                                  k = 1))
  data_xy$female_immature <- log(survey_wide[c(data_xy$nn.idx), 7] + 1) # Match the cpue values and transform
  data_xy$male_immature <- log(survey_wide[c(data_xy$nn.idx), 8] + 1)
  data_xy$male_mature <- log(survey_wide[c(data_xy$nn.idx), 10] + 1)
  data_xy$female_mature <- log(survey_wide[c(data_xy$nn.idx), 11] + 1)
  data_xy$bottom_temp <- survey_wide[c(data_xy$nn.idx), 4] # add bottom temperature (may end up ultimately using ROMS)
  data_xy <- data_xy[-c(33, 34)]
  
  # Convert back to lat, lon
  data_xy <- data_xy[-c(28, 29)]
  data_xy$latitude <- data$latitude[match(data$index, data_xy$index)]
  data_xy$longitude <- data$longitude[match(data$index, data_xy$index)]
  return(data_xy)
}

# Load data ----
# Catch of snow crab in directed fishery
crab_dump <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_crab_dump.csv'))
crab_potsum <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_potsum.csv'))
crab_retained <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_retained_size_freq.csv'))

# Bycatch of snow crab in other crab fisheries
# bycatch_dump <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_crab_dump.csv'))
# bycatch_potsum <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_potsum.csv'))
# bycatch_retained <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_retained_size_freq.csv'))

# EBS survey data
crab_survey <- read_csv(here('data/Snow_CrabData', 'station_cpue_snow.csv'), col_select = -c(1))
crab_dates <- read_csv(here('data/Snow_CrabData', 'year_station_julian_day.csv'))
names(crab_survey) <- tolower(names(crab_survey))
crab_survey <- as.data.frame(crab_survey)
# Add the julian days
crab_survey <- merge(crab_survey, 
                     crab_dates,
                     by.x = c("akfin_survey_year", "gis_station"),
                     by.y = c("year", "station"))
crab_survey <- crab_survey %>%
  filter(akfin_survey_year >= 1995) %>% # observer data only available until this date
  rename(year = akfin_survey_year,
         station = gis_station,
         latitude = mid_latitude,
         longitude = mid_longitude,
         depth = gear_depth,
         temperature = gear_temperature)

# Pivot to wide format
survey_wide <- crab_survey %>%
  pivot_wider(names_from = mat_sex, 
              values_from = cpue)
survey_wide <- as.data.frame(survey_wide)
survey_wide$date <- date.mmddyy(survey_wide$julian)
survey_wide$date <- as.Date(survey_wide$date, "%m/%d/%Y")
survey_wide <- survey_wide %>% 
  mutate(month = lubridate::month(date)) %>% # needed to match sst & ice
  rename(immature_male = 'Immature Male',
         mature_male = 'Mature Male',
         legal_male = 'Legal Male',
         immature_female = 'Immature Female',
         mature_female = 'Mature Female')
survey_wide <- survey_wide[, -13]

# Reformat data (summarize males, split up dates)
crab_detailed <- clean_data(crab_dump)
crab_summary <- clean_data(crab_potsum)
crab_summary$total <- crab_summary$female + crab_summary$tot_legal + crab_summary$sublegal
crab_summary$male <- crab_summary$tot_legal + crab_summary$sublegal
crab_offload <- clean_offload(crab_retained)

# bycatch_detailed <- clean_data(bycatch_dump)
# bycatch_summary <- clean_data(bycatch_potsum)
# bycatch_summary$total <- bycatch_summary$female + bycatch_summary$tot_legal + bycatch_summary$sublegal
# bycatch_summary$male <- bycatch_summary$tot_legal + bycatch_summary$sublegal
# bycatch_offload <- clean_offload(bycatch_retained)

# Need to match up observer data with the survey data
# Group observer data for each season (Nov - Mar)
# Match closest observer data from preceding season to a station for a year
# Get the survey grid
EBS <- get_base_layers(select.region = 'sebs', set.crs = 'auto')
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
                    coords = c("longitude", "latitude"),crs = 4269)
observer_df <- st_join(crab_sf, EBS_trans, left = FALSE)

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


# Separate male and female specimen data ----
# Add index to the specimen data
# Calculate average observer CPUE for each station for each season for each group
# Males 102 and above, < 102 for second group
  
  
  
# Females all together
crab_sorted <- crab_detailed %>%
  group_by(year, doy, adfg, trip) %>%
  mutate(index = cur_group_id()) %>%
  ungroup() %>%
  arrange(index)

crab_filtered <- crab_sorted[-c(1:6, 9, 12:15, 18:22)]

crab_male <- crab_filtered %>%
  filter(sex == 1) %>%
  group_by(index, latitude, longitude, depth, 
           soaktime, maturity, date, year, month, day, doy) %>%
  summarise(immature = sum(size <= 101), 
            mature = sum(size > 101),
            total = immature + mature) 

crab_female <- crab_filtered %>%
  filter(sex != 1) %>%
  group_by(index, latitude, longitude, depth, 
           soaktime, date, year, month, day, doy)  %>%
  summarise(immature = sum(size <= 50), 
            mature = sum(size > 50),
            total = immature + mature)




# Environmental data
# Needs to be matched to the survey data 
# Currently matching the SST and ice to the time period of the survey
sst_data <-as.data.frame(read_csv(here('data', 'sst_latlon.csv'), col_select = -c(1)))
sst_data$month <- match(sst_data$month, month.abb)
ice_data <- as.data.frame(read_csv(here('data', 'ice_latlon.csv'), col_select = -c(1)))
ice_data$month <- match(ice_data$month, month.abb)

phi_raster <- raster(here('data', 'EBS_phi_1km.gri'))
phi_pts <- rasterToPoints(phi_raster, spatial = T)
proj4string(phi_pts) # warning here related to the retirement of rgdal package
phi_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
phi_trans <- spTransform(phi_pts, CRS(phi_prj))
proj4string(phi_trans)
phi_data <- as.data.frame(phi_trans)



