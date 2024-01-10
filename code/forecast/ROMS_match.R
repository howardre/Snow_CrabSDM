# Match Bering10K Forecasts with Survey
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

lon2UTM <- function(longitude){ # Convert to UTM
  (floor((longitude + 180) / 6) %% 60) + 1
}

# Load data ----
# Catch of snow crab in directed fishery
crab_potsum_most <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-1995-2020_potsum.csv'))
crab_potsum_21 <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2021', 'snowcrab-2021_potsum.csv'))

crab_potsum <- rbind(crab_potsum_most, crab_potsum_21)

rm(crab_potsum_most, crab_potsum_21)

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

# Pivot to wide format
survey_wide <- crab_reduced %>%
  pivot_wider(names_from = mat_sex, 
              values_from = c(count))
survey_wide <- as.data.frame(survey_wide)
survey_jdate <- survey_wide %>% mutate(jdate = as.Date(strptime(paste(as.numeric(survey_wide$year),
                                                                      as.numeric(survey_wide$julian)),
                                                                format = "%Y%j"))) %>%
  mutate(date = coalesce(date, jdate))

survey_rename <- survey_jdate %>% 
  rename(sublegal_male = 'Immature Male',
         mature_male = 'Mature Male',
         legal_male = 'Legal Male',
         immature_female = 'Immature Female',
         mature_female = 'Mature Female') %>%
  mutate(sublegal_male = ifelse(is.na(sublegal_male), 0, sublegal_male), # fill NA values with 0a
         mature_male = ifelse(is.na(mature_male), 0, mature_male),
         legal_male = ifelse(is.na(legal_male), 0, legal_male),
         immature_female = ifelse(is.na(immature_female), 0, immature_female),
         mature_female = ifelse(is.na(mature_female), 0, mature_female))
survey_rename <- survey_rename[-20]

# Reformat data (summarize males, split up dates)
crab_summary <- clean_data(crab_potsum)
crab_summary$total <- crab_summary$female + crab_summary$tot_legal + crab_summary$sublegal
crab_summary$male <- crab_summary$tot_legal + crab_summary$sublegal

# Group observer data for each season (Nov - Mar)
# Match closest observer data from preceding season to a station for a year
# Get the survey grid
EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
class(EBS)

# Convert to lat/lon, match catches to station polygons
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")

EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon

crab_sf <- st_as_sf(crab_summary, 
                    coords = c("longitude", "latitude"), 
                    crs = 4269)
survey_sf <- st_as_sf(survey_wide,
                      coords = c("longitude", "latitude"), 
                      crs = 4269)

sf_use_s2(FALSE)
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

# Check for zeros
observer_df$legal_presence <- ifelse(observer_df$tot_legal > 0, 1, 0)
observer_df$sublegal_presence <- ifelse(observer_df$sublegal > 0, 1, 0)
observer_df$female_presence <- ifelse(observer_df$female > 0, 1, 0)

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

# Match to the survey data
survey_combined <- merge(survey_rename, observer_summarized,
                         by.x = c("year", "station"),
                         by.y = c("year_lag", "STATIONID"),
                         all.x = T)[-26] # only 1455 rows with observer data

# Add environmental data ----
# Load data
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
z <- lon2UTM(survey_combined[1, "longitude"])
survey_combined <- tibble::rowid_to_column(survey_combined, "index")

# Add index to summary to match back lat, lon
coordinates(survey_combined) <- c("longitude", "latitude")
proj4string(survey_combined) <- CRS("+proj=longlat +datum=WGS84")

data_xy <- spTransform(survey_combined, CRS(paste0("+proj=utm +zone=", z, " ellps=WGS84")))
data_xy <- as.data.frame(data_xy)

coordinates(phi_data) <- c("coords.x1", "coords.x2")
proj4string(phi_data) <- CRS("+proj=longlat +datum=WGS84")

phi_data_xy <- spTransform(phi_data, CRS(paste0("+proj=utm +zone=", z, " ellps=WGS84")))
phi_data_xy <- as.data.frame(phi_data_xy)

# Match phi, must be a data frame
data_xy[, c(24, 25)] <- as.data.frame(RANN::nn2(phi_data_xy[, c('coords.x2', 'coords.x1')],
                                                data_xy[, c('coords.x2', 'coords.x1')],
                                                k = 1))
data_xy$phi <- phi_data_xy[c(data_xy$nn.idx), 1] # Match nearest phi value
data_xy <- data_xy[-c(24, 25)]

# Match ice to survey stations
crab_ice <- merge(data_xy, ice_final,
                  by.x = c("year", "station"),
                  by.y = c("year", "STATIONID"),
                  all.x = T)

rm(crab_all, crab_env, crab_fill, crab_fix, crab_julian,
   crab_potsum, crab_reduced, crab_summary, crab_survey,
   EBS, EBS_grid, ice_data, ice_data_filtered, ice_df,
   ice_final, ice_means, ice_sf, observer_dates, observer_df,
   observer_summarized, phi_data, phi_data_xy, phi_pts,
   phi_raster, phi_trans, survey_combined, survey_jdate,
   survey_rename, survey_wide, EBS_poly)

# Convert back to lat, lon
crab_final <- crab_ice[-c(5, 7, 9, 22, 23)]
crab_final$latitude <- survey_combined$latitude[match(survey_combined$index, crab_final$index)]
crab_final$longitude <- survey_combined$longitude[match(survey_combined$index, crab_final$index)]

# Match ROMS bottom temperature from hindcast
temp_data <- readRDS(here('data', 'hindcast_temps.rds'))

temp_sf <- st_as_sf(temp_data,
                    coords = c("lon", "lat"), 
                    crs = 4269)

temp_df <- temp_sf %>% st_join(EBS_trans, join = st_intersects, left = FALSE) # match points to stations

temp_df$day <- lubridate::day(temp_df$date)
temp_df$julian <- as.numeric(mdy.date(temp_df$month, temp_df$day, 1960))

temp_means <- temp_df %>%
  group_by(year, STATIONID, julian) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE)) # calculate mean temperature based on station

crab_temp <- merge(crab_final, temp_means,
                   by.x = c("year", "station", "julian"),
                   by.y = c("year", "STATIONID", "julian"),
                   all.x = TRUE)

rm(temp_data, temp_sf, temp_df, temp_means)

# Import bias corrected forecast
forecast_match <- function(roms){

}

cesm_roms <- readRDS(here('data/ROMS', 'cesm_forecast_temp1.rds'))
cesm_126 <- cesm_roms %>% filter(projection == "ssp126", 
                                 year <= format(Sys.Date(), "%Y")) %>% # only need up to current year (the new survey year)
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon)) %>%
  rename(temp_126 = bc) %>%
  dplyr::select(-mo_baseline_value, -projection, -mo_avg_proj, -mean_proj_baseline, -delta)
cesm_585 <- cesm_roms %>% filter(projection == "ssp585", 
                                 year <= format(Sys.Date(), "%Y")) %>% 
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon)) %>%
  rename(temp_585 = bc) %>%
  dplyr::select(-mo_baseline_value, -projection, -mo_avg_proj, -mean_proj_baseline, -delta)

cesm_merge <- merge(as.data.frame(cesm_126), 
                    as.data.frame(cesm_585),
                    by = c("year", "month", "lat", "lon"))
cesm_mean <- cesm_merge %>%
  mutate(cesm_mean = mean(c(temp_126, temp_585))) %>%
  dplyr::select(-temp_126, -temp_585)

cesm_sf <- st_as_sf(cesm_mean,
                    coords = c("lon", "lat"), 
                    crs = 4269)

cesm_df <- cesm_sf %>% st_join(EBS_trans, join = st_intersects, left = FALSE) # match points to stations

cesm_station <- cesm_df %>%
  group_by(year, STATIONID, month) %>%
  summarise(avg_temp = mean(cesm_mean, na.rm = TRUE)) # calculate mean temperature based on station

crab_final$month <- lubridate::month(crab_final$date)

crab_cesm <- merge(crab_final, cesm_station,
                   by.x = c("year", "station", "month"),
                   by.y = c("year", "STATIONID", "month"),
                   all.x = TRUE)


gfdl_roms <- readRDS(here('data/ROMS', 'gfdl_forecast_temp1.rds'))
gfdl_sf <- st_as_sf(gfdl_roms,
                    coords = c("lon", "lat"), 
                    crs = 4269)

miroc_roms <- readRDS(here('data/ROMS', 'miroc_forecast_temp1.rds'))
miroc_sf <- st_as_sf(miroc_roms,
                     coords = c("lon", "lat"), 
                     crs = 4269)




# Interpolate depths
# Use other data to fill in missing depth values
depth_loess <- loess(depth ~ longitude * latitude,
                     span = 0.01,
                     degree = 2,
                     data = crab_temp)
crab_temp$depth_pred <- predict(depth_loess, newdata = crab_temp)

survey_depth <- crab_temp %>% 
  mutate(depth = coalesce(depth, depth_pred))

crab_roms <- survey_depth %>% 
  drop_na(mean_temp, depth, phi, julian) %>% 
  dplyr::select(-geometry.x, -geometry.y)

# Save final data set
saveRDS(crab_roms, file = here('data/Snow_CrabData', 'crab_roms.rds'))
