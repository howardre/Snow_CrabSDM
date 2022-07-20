# Title: Snow Crab Data Exploration
# Purpose: Investigate the three data sources
# Data created: 07/18/2022

# Load libraries ----
library(tidyverse)
library(here)
library(lubridate)
library(date)

# Functions ----
clean_data <- function(data){
  data$date <- as.Date(data$sampdate, "%m-%d-%Y")
  data_frame <- mutate(data,
                       year = lubridate::year(date),
                       month = lubridate::month(date),
                       day = lubridate::day(date))
  data_frame$doy <- as.numeric(mdy.date(data_frame$month, data_frame$day, 1960))
  data_frame[data_frame == -9] <- NA
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
crab_survey <- read_csv(here('data/Snow_CrabData', 'station_cpue_snow.csv'))

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
crab_survey <- crab_survey[-c(1)]

# Visualize
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
