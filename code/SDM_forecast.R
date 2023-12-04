# Title: Next-season forecast SDMs
# Purpose: Used to forecast abundance during survey for following summer
# Data created: 12/04/2023

# Load libraries and functions ----
library(gbm)
library(dismo)
library(here)
library(dplyr)
library(colorspace)
library(maps)
library(mapdata)
library(fields)
library(ggplot2)
library(gplots) # for heatmap
library(enmSdmX) # use for grid search, wrapper for dismo
source(here('code/functions', 'grid_search.R'))
source(here('code/functions', 'rel_inf.R'))
source(here('code/functions', 'rel_inf_fig.R'))
source(here('code/functions', 'part_depen.R'))
source(here('code/functions', 'grid_development.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'brt_grid_preds.R'))

# Load and prepare data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds')) %>%
  dplyr::select(-geometry)

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncount_mat_female = log(mature_female + 1),
                     lncount_imm_female = log(immature_female + 1),
                     lncount_leg_male = log(legal_male + 1),
                     lncount_sub_male = log(sublegal_male + 1),
                     pres_imm_female = ifelse(immature_female > 0, 1, 0),
                     pres_mat_female = ifelse(mature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(sublegal_male > 0, 1, 0),
                     year_f = as.factor(year)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean),
         year_f != 2022) # removed because no observer data

# Create train and test datasets
crab_train <- as.data.frame(crab_trans %>% 
                              filter(year < 2020))
crab_test <- as.data.frame(crab_trans %>% 
                             filter(year == 2021))

# Split into each sex/stage
# Training data
mat_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year)

imm_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year) %>%
  tidyr::drop_na(lncount_sub_male)

# Test data
mat_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year) %>%
  tidyr::drop_na(lncount_mat_female) 

imm_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, longitude, 
                latitude, julian, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year) %>%
  tidyr::drop_na(lncount_sub_male)

# Use LOESS for prediction grid
# Bottom depth based on survey
depth_loess <- loess(depth ~ longitude * latitude,
                     data = crab_trans,
                     span = 0.7,
                     degree = 2,
                     control = loess.control(surface = "interpolate"))
summary(lm(depth_loess$fitted ~ crab_trans$depth)) # check R2

# Grain size based on EBSSED
phi_loess <- loess(phi ~ longitude * latitude,
                   span = 0.05,
                   degree = 2,
                   data = crab_trans,
                   family = "symmetric",
                   control = loess.control(surface = "interpolate"))
summary(lm(phi_loess$fitted ~ crab_trans$phi)) # check R2

# Custom colors
contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
