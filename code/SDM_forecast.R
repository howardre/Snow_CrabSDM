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
source(here('code/functions', 'grid_search_f.R'))
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

# Mature females ----
# Get best models using training data
brt_mat_female_base <- grid_search_f(data = mat_female_train, # uses the trainBRT function in enmSdmX
                                     response = 10, 
                                     family = 'bernoulli')
brt_mat_female_base # check to make sure grid search produced best model, record hyperparameters if needed

brt_mat_female_abun <- grid_search_f(data = mat_female_train[mat_female_train$lncount_mat_female > 0, ],
                                   response = 8, 
                                   family = 'gaussian')
brt_mat_female_abun

# Predict on test data
mat_female_test$pred_base <- predict.gbm(brt_mat_female_base$model,
                                         mat_female_test,
                                         n.trees = brt_mat_female_base$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred_abun <- predict.gbm(brt_mat_female_abun$model,
                                         mat_female_test,
                                         n.trees = brt_mat_female_abun$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred_brt <- mat_female_test$pred_base * mat_female_test$pred_abun

rmse_mat_female_brt <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_brt)^2))
rmse_mat_female_brt

# Calculate deviance explained
brt_deviance(brt_mat_female_abun)
brt_deviance(brt_mat_female_base)

# Spearman correlation coefficient
cor.test(mat_female_test$lncount_mat_female, 
         mat_female_test$pred_brt, 
         method = 'spearman',
         exact = FALSE)

# Save models for future use
saveRDS(brt_mat_female_abun, file = here('data', 'brt_mat_female_abun_forecast.rds'))
saveRDS(brt_mat_female_base, file = here('data', 'brt_mat_female_base_forecast.rds'))

# Mature females ----
# Get best models using training data
brt_imm_female_base <- grid_search_f(data = imm_female_train, # uses the trainBRT function in enmSdmX
                                     response = 10, 
                                     family = 'bernoulli')
brt_imm_female_base # check to make sure grid search produced best model, record hyperparameters if needed

brt_imm_female_abun <- grid_search_f(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                                     response = 8, 
                                     family = 'gaussian')
brt_imm_female_abun

# Predict on test data
imm_female_test$pred_base <- predict.gbm(brt_imm_female_base$model,
                                         imm_female_test,
                                         n.trees = brt_imm_female_base$model$gbm.call$best.trees,
                                         type = "response")

imm_female_test$pred_abun <- predict.gbm(brt_imm_female_abun$model,
                                         imm_female_test,
                                         n.trees = brt_imm_female_abun$model$gbm.call$best.trees,
                                         type = "response")

imm_female_test$pred_brt <- imm_female_test$pred_base * imm_female_test$pred_abun

rmse_imm_female_brt <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_brt)^2))
rmse_imm_female_brt

# Calculate deviance explained
brt_deviance(brt_imm_female_abun)
brt_deviance(brt_imm_female_base)

# Spearman correlation coefficient
cor.test(imm_female_test$lncount_imm_female, 
         imm_female_test$pred_brt, 
         method = 'spearman',
         exact = FALSE)

# Save models for future use
saveRDS(brt_imm_female_abun, file = here('data', 'brt_imm_female_abun_forecast.rds'))
saveRDS(brt_imm_female_base, file = here('data', 'brt_imm_female_base_forecast.rds'))
