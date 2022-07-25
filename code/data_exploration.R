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
# Everything must be a data frame, tables do not work
crab_summary[, c(29, 30)] <- as.data.frame(RANN::nn2(phi_data[, c('y', 'x')],
                                                     crab_summary[, c('latitude', 'longitude')],
                                                     k = 1))
crab_summary$phi <- phi_data[c(crab_summary$nn.idx), 1] # Match nearest phi value
crab_summary <- crab_summary[-c(29, 30)]

crab_summary[, c(30, 31)] <- as.data.frame(RANN::nn2(ice_data[, c('lat', 'lon', 'month', 'year')],
                                                     crab_summary[, c('latitude', 'longitude', 'month', 'year')],
                                                     k = 1))
crab_summary$ice <- ice_data[c(crab_summary$nn.idx), 4] # Match nearest ice value
crab_summary <- crab_summary[-c(30, 31)]

crab_summary[, c(31, 32)] <- as.data.frame(RANN::nn2(sst_data[, c('lat', 'lon', 'month', 'year')],
                                                     crab_summary[, c('latitude', 'longitude', 'month', 'year')],
                                                     k = 1))
crab_summary$sst <- sst_data[c(crab_summary$nn.idx), 4] # Match nearest temperature value
crab_summary <- crab_summary[-c(31, 32)]

# Match survey data ----
survey_wide <- crab_survey %>%
  pivot_wider(names_from = mat_sex, 
              values_from = cpue)
survey_wide <- as.data.frame(survey_wide)

crab_summary[, c(32, 33)] <- as.data.frame(RANN::nn2(survey_wide[, c('mid_latitude', 'mid_longitude', 'akfin_survey_year')],
                                                     crab_summary[, c('latitude', 'longitude', 'year')],
                                                     k = 1))
crab_summary$female_immature <- survey_wide[c(crab_summary$nn.idx), 7] # Match the cpue values
crab_summary$male_immature <- survey_wide[c(crab_summary$nn.idx), 8]
crab_summary$male_mature <- survey_wide[c(crab_summary$nn.idx), 10]
crab_summary$female_mature <- survey_wide[c(crab_summary$nn.idx), 11]
crab_summary <- crab_summary[-c(32, 33)]

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
