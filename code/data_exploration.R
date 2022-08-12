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

# Convert to UTM
lon2UTM <- function(longitude){
  (floor((longitude + 180)/6) %% 60) + 1
}

# Match data for pot level data
match_data <- function(data, survey, phi_data, ice_data,
                       sst_data, survey_wide){
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

survey_wide <- crab_survey_xy %>%
  pivot_wider(names_from = mat_sex, 
              values_from = cpue)
survey_wide <- as.data.frame(survey_wide)

# Match environmental data ----
crab_final <- match_data(crab_summary, crab_survey, 
                         phi_data, ice_data,
                         sst_data, survey_wide)

bycatch_final <- match_data(bycatch_summary, crab_survey,
                            phi_data, ice_data,
                            sst_data, survey_wide)


saveRDS(crab_final, file = here('data/Snow_CrabData', 'crab_summary.rds'))
saveRDS(bycatch_final, file = here('data/Snow_CrabData', 'bycatch_summary.rds'))

# Separate male and female specimen data ----
# Add index to the specimen data
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
  summarise(immature = sum(size <= 90), 
            mature = sum(size > 90),
            total = immature + mature) 

crab_female <- crab_filtered %>%
  filter(sex != 1) %>%
  group_by(index, latitude, longitude, depth, 
           soaktime, date, year, month, day, doy) %>%
  summarise(immature = sum(size <= 50), 
            mature = sum(size > 50),
            total = immature + mature) # do by size or maturity? maturity often not recorded in these data

saveRDS(crab_male, file = here('data/Snow_CrabData', 'crab_male.rds'))
saveRDS(crab_female, file = here('data/Snow_CrabData', 'crab_female.rds'))


# Plot survey data bottom temps to identify gaps ----
# (777 missing values)
bering_sea <- map_data("world")

# Maps of temperatures
ggplot(data = crab_survey, 
       aes(mid_longitude, mid_latitude)) +  
  geom_point(aes(colour = cut(gear_temperature, c(-2.1, 12.5, NA))),
             size = 3) +
  labs(y = "Latitude",
       x = "Longitude") +
  geom_polygon(aes(long, lat, group = group), 
               data = bering_sea,
               fill = "lightyellow4", 
               colour = "black") +
  coord_quickmap(xlim = c(-178, -156), ylim = c(53, 64)) +
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
  facet_wrap(~ akfin_survey_year) +
  labs(colour = "Temperature Records")

# Temperatures over time
survey_wide %>%
  mutate(year_cut = cut(akfin_survey_year, breaks = seq(1975, 2019, by = 1))) %>% 
  ggplot() +
  geom_boxplot(aes(year_cut, gear_temperature),
               fill = "goldenrod1", 
               outlier.alpha = 0.2) +  
  labs(y = "Temperature",
       x = "Year",
       title = "Temperature by Year") +
  scale_x_discrete(labels = c("1975", "1976", "1977", "1978", "1979", "1980", "1981", "1982", 
                              "1983", "1984", "1985", "1986", "1987", "1988", "1989", "1990",
                              "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998",
                              "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006",
                              "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014",
                              "2015", "2016", "2017", "2018", "2019")) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, family = "serif", face = "bold"),
        axis.text = element_text(family = "serif", size = 14),
        axis.title = element_text(family = "serif", size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 18))


# Visualize ----
par(mfrow = c(3, 4))
plot(table(crab_final$year[crab_final$total > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'All')
plot(table(crab_final$month[crab_final$total > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'All')
plot(table(crab_final$soaktime[crab_final$total > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'All')
plot(table(crab_final$gearcode[crab_final$total > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'All')
plot(table(crab_final$year[crab_final$female > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Females')
plot(table(crab_final$month[crab_final$female > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Females')
plot(table(crab_final$soaktime[crab_final$female > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'Females')
plot(table(crab_final$gearcode[crab_final$female > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Females')
plot(table(crab_final$year[crab_final$male > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Males')
plot(table(crab_final$month[crab_final$male > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Males')
plot(table(crab_final$soaktime[crab_final$male > 0]),
     ylab = 'Frequency',
     xlab = 'Soak Time',
     main = 'Males')
plot(table(crab_final$gearcode[crab_final$male > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Males')


test_gam <- gam(mature + 1 ~ factor(year) +
                  s(latitude, longitude) +
                  s(doy) +
                  s(depth),
                data = crab_male,
                family = tw(link = "log"))
summary(test_gam)
