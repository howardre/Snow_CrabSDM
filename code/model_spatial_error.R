# Title: Snow Crab Forecast Models
# Purpose: Reformat to use variables for forecast
# Data created: 05/04/2023

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(here)
library(dplyr)
library(colorspace)
library(fields)
library(ggplot2)
library(akgfmaps)
library(enmSdmX) # use for grid search, wrapper for dismo
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'grid_search.R'))
source(here('code/functions', 'map_rmse.R'))
source(here('code/functions', 'brt_grid_preds.R'))

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

# Load data ----
# Get survey grid
EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")
EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon
bering_sea <- map_data("world")

# Make sure to run PCA first if updating the data matching script
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
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean),
         year_f != 2022)

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
# Need to filter out stations without observer data in order to get accurate comparison
crab_train <- as.data.frame(crab_trans %>% 
                              filter(year < 2015))
crab_test <- as.data.frame(crab_trans %>% 
                             filter(year > 2014))

crab_train_warm <- as.data.frame(crab_trans %>% 
                              filter(year < 2015 | year == 2018))
crab_test_warm <- as.data.frame(crab_trans %>% 
                             filter(year > 2014 & year != 2018))

# Split into each sex/stage
# Training data
mat_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_mat_female) 

mat_female_train_warm <- crab_train_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_mat_female)

imm_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_imm_female)

imm_female_train_warm <- crab_train_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_leg_male)

leg_male_train_warm <- crab_train_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

sub_male_train_warm <- crab_train_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

# Test data
mat_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_mat_female)

mat_female_test_warm <- crab_test_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_mat_female) 

imm_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_imm_female)

imm_female_test_warm <- crab_test_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station,
                station) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_leg_male)

leg_male_test_warm <- crab_test_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

sub_male_test_warm <- crab_test_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station,
                station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

# Boosted regression trees ----
# Adjust the bag fraction to a value between 0.5-0.75 as suggested by Elith et al. (2008)
# The learning rate could range from 0.1-0.0001, higher value usually means less trees
# Depending on the number of samples, want tree complexity to be high enough (likely using 5)
# Want at least 1000 trees, but don't need to go way beyond it
vars <- c(1:8, 10, 16)

## Mature females ----
# Get best models
# Read in BRTs
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base.rds'))

brt_mat_female_base_warm <- grid_search(mat_female_train_warm, 13, 'bernoulli')
brt_mat_female_base_warm

brt_mat_female_abun_warm <- grid_search(mat_female_train_warm[mat_female_train_warm$lncount_mat_female > 0,],
                                        11, 'gaussian')
brt_mat_female_abun_warm

# Predict on test data
# Original model
mat_female_test$pred_base <- predict.gbm(brt_mat_female_base$model,
                                         mat_female_test,
                                         n.trees = brt_mat_female_base$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred_abun <- predict.gbm(brt_mat_female_abun$model,
                                         mat_female_test,
                                         n.trees = brt_mat_female_abun$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred_brt <- mat_female_test$pred_base * mat_female_test$pred_abun

# Model trained with warm year
mat_female_test_warm$pred_base <- predict.gbm(brt_mat_female_base_warm$model,
                                              mat_female_test_warm,
                                              n.trees = brt_mat_female_base_warm$model$gbm.call$best.trees,
                                              type = "response")

mat_female_test_warm$pred_abun <- predict.gbm(brt_mat_female_abun_warm$model,
                                              mat_female_test_warm,
                                              n.trees = brt_mat_female_abun_warm$model$gbm.call$best.trees,
                                              type = "response")

mat_female_test_warm$pred_brt <- mat_female_test_warm$pred_base * mat_female_test_warm$pred_abun

# Calculate RMSE
# 1.62 for base model
rmse_mat_female_brt <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_brt)^2))
rmse_mat_female_brt # 1.68
rmse_mat_female_brt_warm <- sqrt(mean((mat_female_test_warm$lncount_mat_female - mat_female_test_warm$pred_brt)^2))
rmse_mat_female_brt_warm # 1.55

# Calculate deviance explained
dev_mat_female_abun_warm <- brt_deviance(brt_mat_female_abun_warm)
dev_mat_female_pres_warm <- brt_deviance(brt_mat_female_base_warm)

dev_mat_female_abun_warm # 39.8% deviance explained
dev_mat_female_pres_warm # 56.0% deviance explained

# Save models for future use
saveRDS(brt_mat_female_abun_warm, file = here('data', 'brt_mat_female_abun_warm.rds'))
saveRDS(brt_mat_female_base_warm, file = here('data', 'brt_mat_female_pres_warm.rds'))

# Spatial RMSE
# Map RMSE
# Need to get average RMSE per station
# Then plot by station to get spatial error
brt_mat_female_abun_warm <- readRDS(file = here('data', 'brt_mat_female_abun_warm.rds'))
brt_mat_female_base_warm <- readRDS(file = here('data', 'brt_mat_female_pres_warm.rds'))

mat_female_rmse_warm <- mat_female_test_warm %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_mat_female, pred_brt)) # calculate by station

mat_female_rmse <- mat_female_test %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_mat_female, pred_brt))

mat_female_rmse$change <- ((mat_female_rmse$rmse - mat_female_rmse_warm$rmse) / 
  mat_female_rmse$rmse) * 100

map_rmse(mat_female_rmse_warm, "Spatial Error for Mature Female Crab - Warm Year")
dev.copy(jpeg,
         here('results/RMSE',
              'mat_female_rmse_warm.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_rmse(mat_female_rmse, "Spatial Error for Mature Female Crab")
dev.copy(jpeg,
         here('results/RMSE',
              'mat_female_rmse_base.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_change(mat_female_rmse)
dev.copy(jpeg,
         here('results/RMSE',
              'mat_female_rmse_change.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

## Immature females ----
# Get best models
# Read in BRTs
brt_imm_female_abun <- readRDS(file = here('data', 'brt_imm_female_abun.rds'))
brt_imm_female_base <- readRDS(file = here('data', 'brt_imm_female_base.rds'))

brt_imm_female_base_warm <- grid_search(imm_female_train_warm, 13, 'bernoulli')
brt_imm_female_base_warm

brt_imm_female_abun_warm <- grid_search(imm_female_train_warm[imm_female_train_warm$lncount_imm_female > 0,], 
                                        11, 'gaussian')
brt_imm_female_abun_warm

# using the converged model hyperparameters to re-run
# No issues with minimum trees or convergence so not clear why this is happening
# brt_imm_female_model <- gbm.step(data = imm_female_train_warm[imm_female_train_warm$lncount_imm_female > 0,],
#                                  gbm.x = vars,
#                                  gbm.y = 11,
#                                  family = "gaussian",
#                                  learning.rate = 0.01,
#                                  tree.complexity = 11,
#                                  bag.fraction = 0.75,
#                                  max.trees = 3000,
#                                  step.size = 40)
# brt_imm_female_abun_warm <- list(brt_imm_female_model, brt_imm_female_abun_warm$tuning)
# names(brt_imm_female_abun_warm)[[1]] <- "model"
# names(brt_imm_female_abun_warm)[[2]] <- "tuning"

# Predict on test data
# Original model
imm_female_test$pred_base <- predict.gbm(brt_imm_female_base$model,
                                         imm_female_test,
                                         n.trees = brt_imm_female_base$model$gbm.call$best.trees,
                                         type = "response")

imm_female_test$pred_abun <- predict.gbm(brt_imm_female_abun$model,
                                         imm_female_test,
                                         n.trees = brt_imm_female_abun$model$gbm.call$best.trees,
                                         type = "response")

imm_female_test$pred_brt <- imm_female_test$pred_base * imm_female_test$pred_abun

# Model trained with warm year
imm_female_test_warm$pred_base <- predict.gbm(brt_imm_female_base_warm$model,
                                              imm_female_test_warm,
                                              n.trees = brt_imm_female_base_warm$model$gbm.call$best.trees,
                                              type = "response")

imm_female_test_warm$pred_abun <- predict.gbm(brt_imm_female_abun_warm$model,
                                              imm_female_test_warm,
                                              n.trees = brt_imm_female_abun_warm$model$gbm.call$best.trees,
                                              type = "response")

imm_female_test_warm$pred_brt <- imm_female_test_warm$pred_base * imm_female_test_warm$pred_abun


# Calculate RMSE
rmse_imm_female_brt <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_brt)^2))
rmse_imm_female_brt # 1.41
rmse_imm_female_brt_warm <- sqrt(mean((imm_female_test_warm$lncount_imm_female - imm_female_test_warm$pred_brt)^2))
rmse_imm_female_brt_warm # 1.43

# Calculate deviance explained
dev_imm_female_abun_warm <- brt_deviance(brt_imm_female_abun_warm)
dev_imm_female_pres_warm <- brt_deviance(brt_imm_female_base_warm)

dev_imm_female_abun_warm # 46.3% deviance explained
dev_imm_female_pres_warm # 46.6% deviance explained

# Save models for future use
saveRDS(brt_imm_female_abun_warm, file = here('data', 'brt_imm_female_abun_warm.rds'))
saveRDS(brt_imm_female_base_warm, file = here('data', 'brt_imm_female_pres_warm.rds'))

# Spatial RMSE
# Map RMSE
# Need to get average RMSE per station
# Then plot by station to get spatial error
brt_imm_female_abun_warm <- readRDS(file = here('data', 'brt_imm_female_abun_warm.rds'))
brt_imm_female_base_warm <- readRDS(file = here('data', 'brt_imm_female_pres_warm.rds'))

imm_female_rmse_warm <- imm_female_test_warm %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_imm_female, pred_brt)) # calculate by station

imm_female_rmse <- imm_female_test %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_imm_female, pred_brt))

imm_female_rmse$change <- ((imm_female_rmse$rmse - imm_female_rmse_warm$rmse) / 
                             imm_female_rmse$rmse) * 100

map_rmse(imm_female_rmse_warm, "Spatial Error for Immature Female Crab - Warm Year")
dev.copy(jpeg,
         here('results/RMSE',
              'imm_female_rmse_warm.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_rmse(imm_female_rmse, "Spatial Error for Immature Female Crab")
dev.copy(jpeg,
         here('results/RMSE',
              'imm_female_rmse_base.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_change(imm_female_rmse)
dev.copy(jpeg,
         here('results/RMSE',
              'imm_female_rmse_change.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

## Legal Males ----
# Get best models
# Read in BRTs
brt_leg_male_abun <- readRDS(file = here('data', 'brt_leg_male_abun.rds'))
brt_leg_male_base <- readRDS(file = here('data', 'brt_leg_male_base.rds'))

brt_leg_male_base_warm <- grid_search(leg_male_train_warm, 13, 'bernoulli')
brt_leg_male_base_warm

brt_leg_male_abun_warm <- grid_search(leg_male_train_warm[leg_male_train_warm$lncount_leg_male > 0,],
                                        11, 'gaussian')
brt_leg_male_abun_warm

# Predict on test data
# Original model
leg_male_test$pred_base <- predict.gbm(brt_leg_male_base$model,
                                         leg_male_test,
                                         n.trees = brt_leg_male_base$model$gbm.call$best.trees,
                                         type = "response")

leg_male_test$pred_abun <- predict.gbm(brt_leg_male_abun$model,
                                         leg_male_test,
                                         n.trees = brt_leg_male_abun$model$gbm.call$best.trees,
                                         type = "response")

leg_male_test$pred_brt <- leg_male_test$pred_base * leg_male_test$pred_abun

# Model trained with warm year
leg_male_test_warm$pred_base <- predict.gbm(brt_leg_male_base_warm$model,
                                              leg_male_test_warm,
                                              n.trees = brt_leg_male_base_warm$model$gbm.call$best.trees,
                                              type = "response")

leg_male_test_warm$pred_abun <- predict.gbm(brt_leg_male_abun_warm$model,
                                              leg_male_test_warm,
                                              n.trees = brt_leg_male_abun_warm$model$gbm.call$best.trees,
                                              type = "response")

leg_male_test_warm$pred_brt <- leg_male_test_warm$pred_base * leg_male_test_warm$pred_abun


# Calculate RMSE
rmse_leg_male_brt <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_brt)^2))
rmse_leg_male_brt # 1.25
rmse_leg_male_brt_warm <- sqrt(mean((leg_male_test_warm$lncount_leg_male - leg_male_test_warm$pred_brt)^2))
rmse_leg_male_brt_warm # 1.15

# Calculate deviance explained
dev_leg_male_abun_warm <- brt_deviance(brt_leg_male_abun_warm)
dev_leg_male_pres_warm <- brt_deviance(brt_leg_male_base_warm)

dev_leg_male_abun_warm # 48.5% deviance explained
dev_leg_male_pres_warm # 53.8% deviance explained

# Save models for future use
saveRDS(brt_leg_male_abun_warm, file = here('data', 'brt_leg_male_abun_warm.rds'))
saveRDS(brt_leg_male_base_warm, file = here('data', 'brt_leg_male_pres_warm.rds'))

# Spatial RMSE
# Map RMSE
# Need to get average RMSE per station
# Then plot by station to get spatial error
brt_leg_male_abun_warm <- readRDS(file = here('data', 'brt_leg_male_abun_warm.rds'))
brt_leg_male_base_warm <- readRDS(file = here('data', 'brt_leg_male_pres_warm.rds'))

leg_male_rmse_warm <- leg_male_test_warm %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_leg_male, pred_brt)) # calculate by station

leg_male_rmse <- leg_male_test %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_leg_male, pred_brt))

leg_male_rmse$change <- ((leg_male_rmse$rmse - leg_male_rmse_warm$rmse) / 
                             leg_male_rmse$rmse) * 100

map_rmse(leg_male_rmse_warm, "Spatial Error for Legal Male Crab - Warm Year")
dev.copy(jpeg,
         here('results/RMSE',
              'leg_male_rmse_warm.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_rmse(leg_male_rmse, "Spatial Error for Legal Male Crab")
dev.copy(jpeg,
         here('results/RMSE',
              'leg_male_rmse_base.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_change(leg_male_rmse)
dev.copy(jpeg,
         here('results/RMSE',
              'leg_male_rmse_change.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

## Sublegal Males ----
# Get best models
# Read in BRTs
brt_sub_male_abun <- readRDS(file = here('data', 'brt_sub_male_abun.rds'))
brt_sub_male_base <- readRDS(file = here('data', 'brt_sub_male_base.rds'))

brt_sub_male_base_warm <- grid_search(sub_male_train_warm, 13, 'bernoulli')
brt_sub_male_base_warm

brt_sub_male_abun_warm <- grid_search(sub_male_train_warm[sub_male_train_warm$lncount_sub_male > 0,],
                                        11, 'gaussian')
brt_sub_male_abun_warm

# Predict on test data
# Original model
sub_male_test$pred_base <- predict.gbm(brt_sub_male_base$model,
                                         sub_male_test,
                                         n.trees = brt_sub_male_base$model$gbm.call$best.trees,
                                         type = "response")

sub_male_test$pred_abun <- predict.gbm(brt_sub_male_abun$model,
                                         sub_male_test,
                                         n.trees = brt_sub_male_abun$model$gbm.call$best.trees,
                                         type = "response")

sub_male_test$pred_brt <- sub_male_test$pred_base * sub_male_test$pred_abun

# Model trained with warm year
sub_male_test_warm$pred_base <- predict.gbm(brt_sub_male_base_warm$model,
                                              sub_male_test_warm,
                                              n.trees = brt_sub_male_base_warm$model$gbm.call$best.trees,
                                              type = "response")

sub_male_test_warm$pred_abun <- predict.gbm(brt_sub_male_abun_warm$model,
                                              sub_male_test_warm,
                                              n.trees = brt_sub_male_abun_warm$model$gbm.call$best.trees,
                                              type = "response")

sub_male_test_warm$pred_brt <- sub_male_test_warm$pred_base * sub_male_test_warm$pred_abun


# Calculate RMSE
rmse_sub_male_brt <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_brt)^2))
rmse_sub_male_brt # 1.60
rmse_sub_male_brt_warm <- sqrt(mean((sub_male_test_warm$lncount_sub_male - sub_male_test_warm$pred_brt)^2))
rmse_sub_male_brt_warm # 1.50

# Calculate deviance explained
dev_sub_male_abun_warm <- brt_deviance(brt_sub_male_abun_warm)
dev_sub_male_pres_warm <- brt_deviance(brt_sub_male_base_warm)

dev_sub_male_abun_warm # 64.9% deviance explained
dev_sub_male_pres_warm # 59.6% deviance explained

# Save models for future use
saveRDS(brt_sub_male_abun_warm, file = here('data', 'brt_sub_male_abun_warm.rds'))
saveRDS(brt_sub_male_base_warm, file = here('data', 'brt_sub_male_pres_warm.rds'))

# Spatial RMSE
# Map RMSE
# Need to get average RMSE per station
# Then plot by station to get spatial error
brt_sub_male_abun_warm <- readRDS(file = here('data', 'brt_sub_male_abun_warm.rds'))
brt_sub_male_base_warm <- readRDS(file = here('data', 'brt_sub_male_pres_warm.rds'))

sub_male_rmse_warm <- sub_male_test_warm %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_sub_male, pred_brt)) # calculate by station

sub_male_rmse <- sub_male_test %>%
  group_by(station) %>%
  summarize(rmse = Metrics::rmse(lncount_sub_male, pred_brt))

sub_male_rmse$change <- ((sub_male_rmse$rmse - sub_male_rmse_warm$rmse) / 
                           sub_male_rmse$rmse) * 100

map_rmse(sub_male_rmse_warm, "Spatial Error for Sublegal Male Crab - Warm Year")
dev.copy(jpeg,
         here('results/RMSE',
              'sub_male_rmse_warm.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_rmse(sub_male_rmse, "Spatial Error for Sublegal Male Crab")
dev.copy(jpeg,
         here('results/RMSE',
              'sub_male_rmse_base.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

map_change(sub_male_rmse)
dev.copy(jpeg,
         here('results/RMSE',
              'sub_male_rmse_change.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()