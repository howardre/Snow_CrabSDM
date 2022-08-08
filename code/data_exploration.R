# Title: Snow Crab Data Exploration
# Purpose: Investigate the three data sources
# Data created: 07/18/2022

# Load libraries ----
library(tidyverse)
library(here)
library(lubridate)
library(date)
library(RANN)
library(raster)
library(sp)

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

lon2UTM <- function(longitude){
  (floor((longitude + 180)/6) %% 60) + 1
}

# Load data ----
# Catch of snow crab in directed fishery
crab_dump <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_crab_dump.csv'))
crab_potsum <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_potsum.csv'))
crab_retained <- read_csv(here('data/Snow_CrabData/snowcrab-1995-2020', 'snowcrab-1995-2020_retained_size_freq.csv'))

# Bycatch of snow crab in other crab fisheries
bycatch_dump <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_crab_dump.csv'))
bycatch_potsum <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_potsum.csv'))
bycatch_retained <- read_csv(here('data/Snow_CrabData/snowcrab_bycatch-1995-2020', 'snowcrab_bycatch-1995-2020_retained_size_freq.csv'))

# EBS survey data
crab_survey <- read_csv(here('data/Snow_CrabData', 'station_cpue_snow.csv'), col_select = -c(1))

# Environmental data
sst_data <-as.data.frame(read_csv(here('data', 'sst_latlon.csv'), col_select = -c(1)))
sst_data$month <- match(sst_data$month, month.abb)
ice_data <- as.data.frame(read_csv(here('data', 'ice_latlon.csv'), col_select = -c(1)))
ice_data$month <- match(ice_data$month, month.abb)

phi_raster <- raster(here('data', 'EBS_phi_1km.gri'))
phi_pts <- rasterToPoints(phi_raster, spatial = T)
proj4string(phi_pts)
phi_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
phi_trans <- spTransform(phi_pts, CRS(phi_prj))
proj4string(phi_trans)
phi_data <- as.data.frame(phi_trans)

# Reformat data (summarize males, split up dates)
crab_detailed <- clean_data(crab_dump)
crab_summary <- clean_data(crab_potsum)
crab_summary$total <- crab_summary$female + crab_summary$tot_legal + crab_summary$sublegal
crab_summary$male <- crab_summary$tot_legal + crab_summary$sublegal
crab_offload <- clean_data(crab_retained)

bycatch_detailed <- clean_data(bycatch_dump)
bycatch_summary <- clean_data(bycatch_potsum)
bycatch_summary$total <- bycatch_summary$female + bycatch_summary$tot_legal + bycatch_summary$sublegal
bycatch_summary$male <- bycatch_summary$tot_legal + bycatch_summary$sublegal
bycatch_offload <- clean_data(bycatch_retained)

names(crab_survey) <- tolower(names(crab_survey))
crab_survey <- as.data.frame(crab_survey)

# Match environmental data ----
# Convert to UTM in order to properly match 
z <- lon2UTM(crab_summary[1, "longitude"])

crab_survey1 <- crab_survey
coordinates(crab_survey1) <- c("mid_longitude", "mid_latitude")
proj4string(crab_survey1) <- CRS("+proj=longlat +datum=WGS84")

crab_survey_xy <- spTransform(crab_survey1, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
crab_survey_xy <- as.data.frame(crab_survey_xy)

crab_summary1 <- crab_summary
coordinates(crab_summary1) <- c("longitude", "latitude")
proj4string(crab_summary1) <- CRS("+proj=longlat +datum=WGS84")

crab_summary_xy <- spTransform(crab_summary1, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
crab_summary_xy <- as.data.frame(crab_summary_xy)

phi_data1 <- phi_data
coordinates(phi_data1) <- c("x", "y")
proj4string(phi_data1) <- CRS("+proj=longlat +datum=WGS84")

phi_data_xy <- spTransform(phi_data1, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
phi_data_xy <- as.data.frame(phi_data_xy)

ice_data1 <- ice_data
coordinates(ice_data1) <- c("lon", "lat")
proj4string(ice_data1) <- CRS("+proj=longlat +datum=WGS84")

ice_data_xy <- spTransform(ice_data1, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
ice_data_xy <- as.data.frame(ice_data_xy)

sst_data1 <- sst_data
coordinates(sst_data1) <- c("lon", "lat")
proj4string(sst_data1) <- CRS("+proj=longlat +datum=WGS84")

sst_data_xy <- spTransform(sst_data1, CRS(paste0("+proj=utm +zone=", z, "ellps=WGS84")))
sst_data_xy <- as.data.frame(sst_data_xy)

# Everything must be a data frame, tables do not work
crab_summary_xy[, c(29, 30)] <- as.data.frame(RANN::nn2(phi_data_xy[, c('y', 'x')],
                                                     crab_summary_xy[, c('latitude', 'longitude')],
                                                     k = 1))
crab_summary_xy$phi <- phi_data_xy[c(crab_summary_xy$nn.idx), 1] # Match nearest phi value
crab_summary_xy <- crab_summary_xy[-c(29, 30)]

crab_summary_xy[, c(30, 31)] <- as.data.frame(RANN::nn2(ice_data_xy[, c('lat', 'lon', 'month', 'year')],
                                                     crab_summary_xy[, c('latitude', 'longitude', 'month', 'year')],
                                                     k = 1))
crab_summary_xy$ice <- ice_data_xy[c(crab_summary_xy$nn.idx), 2] # Match nearest ice value
crab_summary_xy <- crab_summary_xy[-c(30, 31)]

crab_summary_xy[, c(31, 32)] <- as.data.frame(RANN::nn2(sst_data_xy[, c('lat', 'lon', 'month', 'year')],
                                                     crab_summary_xy[, c('latitude', 'longitude', 'month', 'year')],
                                                     k = 1))
crab_summary_xy$sst <- sst_data_xy[c(crab_summary_xy$nn.idx), 2] # Match nearest temperature value
crab_summary_xy <- crab_summary_xy[-c(31, 32)]

# Match survey data ----
survey_wide <- crab_survey_xy %>%
  pivot_wider(names_from = mat_sex, 
              values_from = cpue)
survey_wide <- as.data.frame(survey_wide)

crab_summary_xy[, c(32, 33)] <- as.data.frame(RANN::nn2(survey_wide[, c('mid_latitude', 'mid_longitude', 'akfin_survey_year')],
                                                     crab_summary_xy[, c('latitude', 'longitude', 'year')],
                                                     k = 1))
crab_summary_xy$female_immature <- log(survey_wide[c(crab_summary_xy$nn.idx), 7] + 1) # Match the cpue values and transform
crab_summary_xy$male_immature <- log(survey_wide[c(crab_summary_xy$nn.idx), 8] + 1)
crab_summary_xy$male_mature <- log(survey_wide[c(crab_summary_xy$nn.idx), 10] + 1)
crab_summary_xy$female_mature <- log(survey_wide[c(crab_summary_xy$nn.idx), 11] + 1)
crab_summary_xy$bottom_temp <- survey_wide[c(crab_summary_xy$nn.idx, na.omit = T), 6] # add bottom temperature (may end up ultimately using ROMS)
crab_summary_xy <- crab_summary_xy[-c(32, 33)]


saveRDS(crab_summary, file = here('data/Snow_CrabData', 'crab_summary.rds'))

# Visualize ----
par(mfrow = c(3, 4))
plot(table(crab_summary$year[crab_summary$total > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'All')
plot(table(crab_summary$month[crab_summary$total > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'All')
plot(table(crab_summary$soaktime[crab_summary$total > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'All')
plot(table(crab_summary$gearcode[crab_summary$total > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'All')
plot(table(crab_summary$year[crab_summary$female > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Females')
plot(table(crab_summary$month[crab_summary$female > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Females')
plot(table(crab_summary$soaktime[crab_summary$female > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'Females')
plot(table(crab_summary$gearcode[crab_summary$female > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Females')
plot(table(crab_summary$year[crab_summary$male > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Males')
plot(table(crab_summary$month[crab_summary$male > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Males')
plot(table(crab_summary$soaktime[crab_summary$male > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'Males')
plot(table(crab_summary$gearcode[crab_summary$male > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Males')
