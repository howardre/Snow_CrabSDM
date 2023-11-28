# Title: Snow Crab Sex- and Maturity-Specific SDMs
# Purpose: Final models to be used for Bering Sea snow crab, to be rerun if updated models needed
# Data created: 11/28/2023

# Load libraries and functions ----
library(gbm)
library(dismo)
# library(scales)
library(here)
library(dplyr)
library(colorspace)
library(maps)
library(mapdata)
library(fields)
library(ggplot2)
# library(data.table)
# library(RColorBrewer)
# library(gplots) # for heatmap
# library(enmSdmX) # use for grid search, wrapper for dismo
# library(fastshap) # calculate SHAP values quickly
# library(mshap) # combine SHAP values for two-part models
# library(shapviz) # visualize SHAP values
# library(doParallel)
source(here('code/functions', 'grid_search.R'))
source(here('code/functions', 'rel_inf.R'))
source(here('code/functions', 'part_depen.R'))
source(here('code/functions', 'grid_development.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'brt_grid_preds.R'))
source(here('code/functions', 'map_pred_brt.R'))

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
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean),
         year_f != 2022) # removed because no observer data

# Create train and test datasets
crab_train <- as.data.frame(crab_trans %>% 
                              filter(year < 2015))
crab_test <- as.data.frame(crab_trans %>% 
                             filter(year > 2014))

# Split into each sex/stage
# Training data
mat_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station)

imm_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

# Test data
mat_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_mat_female) 

imm_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station) %>%
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

# Boosted regression trees ----
# See model_evaluation.R script for more details on hyperparameters, tuning
vars <- c(1:8, 10, 16) # model covariates for each model, in the grid_search function

## Mature females ----
# Get best models using training data
brt_mat_female_base <- grid_search(mat_female_train, 13, 'bernoulli')
brt_mat_female_base

brt_mat_female_abun <- grid_search(mat_female_train[mat_female_train$lncount_mat_female > 0, ],
                                   11, 'gaussian')
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


# Calculate RMSE
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
saveRDS(brt_mat_female_abun, file = here('data', 'brt_mat_female_abun.rds'))
saveRDS(brt_mat_female_base, file = here('data', 'brt_mat_female_base.rds'))

# Read in BRTs
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base.rds'))

# Variable importance
rel_inf(brt_mat_female_base, 'Variable Influence on Mature Female Snow Crab')

# Plot partial dependence
part_depen(brt_mat_female_base)

# Plot map of the predicted distribution
spatial_grid_mat_female <- grid_development(mat_female_train) # create a prediction grid
spatial_grid_mat_female$female_loading_station <- median(mat_female_train$female_loading_station, na.rm = T)
spatial_grid_mat_female$bcs_mature_female <- median(mat_female_train$bcs_mature_female, na.rm = T)
mat_female_preds <- brt_grid_preds(spatial_grid_mat_female, 
                                   brt_mat_female_abun,
                                   brt_mat_female_base)

map_pred_brt(mat_female_preds, mat_female_train, "Mature Female Snow Crab") # works better when saved as .jpeg


## Legal males ----
# Get best models using training data
brt_leg_male_base <- grid_search(leg_male_train, 13, 'bernoulli')
brt_leg_male_base

brt_leg_male_abun <- grid_search(leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                                 11, 'gaussian')
brt_leg_male_abun

# Predict on test data
leg_male_test$pred_base <- predict.gbm(brt_leg_male_base$model,
                                       leg_male_test,
                                       n.trees = brt_leg_male_base$model$gbm.call$best.trees,
                                       type = "response")

leg_male_test$pred_abun <- predict.gbm(brt_leg_male_abun$model,
                                       leg_male_test,
                                       n.trees = brt_leg_male_abun$model$gbm.call$best.trees,
                                       type = "response")

leg_male_test$pred_brt <- leg_male_test$pred_base * leg_male_test$pred_abun

# Calculate RMSE
rmse_leg_male_brt <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_brt)^2))
rmse_leg_male_brt

# Calculate deviance explained
brt_deviance(brt_leg_male_abun)
brt_deviance(brt_leg_male_base)

# Spearman correlation coefficient
cor.test(leg_male_test$lncount_leg_male, 
         leg_male_test$pred_brt, 
         method = 'spearman',
         exact = FALSE)

# Save models for future use
saveRDS(brt_leg_male_abun, file = here('data', 'brt_leg_male_abun.rds'))
saveRDS(brt_leg_male_base, file = here('data', 'brt_leg_male_base.rds'))

# Read in BRTs
brt_leg_male_abun <- readRDS(file = here('data', 'brt_leg_male_abun.rds'))
brt_leg_male_base <- readRDS(file = here('data', 'brt_leg_male_base.rds'))

# Variable importance
rel_inf(brt_leg_male_base, 'Variable Influence on Legal Male Snow Crab')

# Plot partial dependence
part_depen(brt_leg_male_base)

# Plot map of the predicted distribution
spatial_grid_leg_male <- grid_development(leg_male_train) # create a prediction grid
spatial_grid_leg_male$legal_male_loading_station <- median(leg_male_train$legal_male_loading_station, na.rm = T)
spatial_grid_leg_male$bcs_legal_male <- median(leg_male_train$bcs_legal_male, na.rm = T)
leg_male_preds <- brt_grid_preds(spatial_grid_leg_male,
                                 brt_leg_male_abun,
                                 brt_leg_male_base)

map_pred_brt(leg_male_preds, leg_male_train, "Legal Male Snow Crab") # works better when saved as .jpeg


## Immature females ----
# Get best models using training data
brt_imm_female_base <- grid_search(imm_female_train, 13, 'bernoulli')
brt_imm_female_base

brt_imm_female_abun <- grid_search(imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                                   11, 'gaussian')
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

# Calculate RMSE
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
saveRDS(brt_imm_female_abun, file = here('data', 'brt_imm_female_abun.rds'))
saveRDS(brt_imm_female_base, file = here('data', 'brt_imm_female_base.rds'))

# Read in BRTs
brt_imm_female_abun <- readRDS(file = here('data', 'brt_imm_female_abun.rds'))
brt_imm_female_base <- readRDS(file = here('data', 'brt_imm_female_base.rds'))

# Variable importance
rel_inf(brt_imm_female_base, 'Variable Influence on Immature Female Snow Crab')

# Plot partial dependence
part_depen(brt_imm_female_base)

# Plot map of the predicted distribution
spatial_grid_imm_female <- grid_development(imm_female_train) # create a prediction grid
spatial_grid_imm_female$female_loading_station <- median(imm_female_train$female_loading_station, na.rm = T)
spatial_grid_imm_female$bcs_immature_female <- median(imm_female_train$bcs_immature_female, na.rm = T)
imm_female_preds <- brt_grid_preds(spatial_grid_imm_female, 
                                   brt_imm_female_abun,
                                   brt_imm_female_base)

map_pred_brt(imm_female_preds, imm_female_train, "Immature Female Snow Crab") # works better when saved as .jpeg


## Sublegal males ----
# Get best models using training data
brt_sub_male_base <- grid_search(sub_male_train, 13, 'bernoulli')
brt_sub_male_base

brt_sub_male_abun <- grid_search(sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                                 11, 'gaussian')
brt_sub_male_abun

# Predict on test data
sub_male_test$pred_base <- predict.gbm(brt_sub_male_base$model,
                                       sub_male_test,
                                       n.trees = brt_sub_male_base$model$gbm.call$best.trees,
                                       type = "response")

sub_male_test$pred_abun <- predict.gbm(brt_sub_male_abun$model,
                                       sub_male_test,
                                       n.trees = brt_sub_male_abun$model$gbm.call$best.trees,
                                       type = "response")

sub_male_test$pred_brt <- sub_male_test$pred_base * sub_male_test$pred_abun

# Calculate RMSE
rmse_sub_male_brt <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_brt)^2))
rmse_sub_male_brt

# Calculate deviance explained
brt_deviance(brt_sub_male_abun)
brt_deviance(brt_sub_male_base)

# Spearman correlation coefficient
cor.test(sub_male_test$lncount_sub_male, 
         sub_male_test$pred_brt, 
         method = 'spearman',
         exact = FALSE)

# Save models for future use
saveRDS(brt_sub_male_abun, file = here('data', 'brt_sub_male_abun.rds'))
saveRDS(brt_sub_male_base, file = here('data', 'brt_sub_male_base.rds'))

# Read in BRTs
brt_sub_male_abun <- readRDS(file = here('data', 'brt_sub_male_abun.rds'))
brt_sub_male_base <- readRDS(file = here('data', 'brt_sub_male_base.rds'))

# Variable importance
rel_inf(brt_sub_male_base, 'Variable Influence on Sublegal Male Snow Crab')

# Plot partial dependence
part_depen(brt_sub_male_base)

# Plot map of the predicted distribution
spatial_grid_sub_male <- grid_development(sub_male_train) # create a prediction grid
spatial_grid_sub_male$sublegal_male_loading_station <- median(sub_male_train$sublegal_male_loading_station, na.rm = T)
spatial_grid_sub_male$bcs_sublegal_male <- median(sub_male_train$bcs_sublegal_male, na.rm = T)
sub_male_preds <- brt_grid_preds(spatial_grid_sub_male,
                                 brt_sub_male_abun,
                                 brt_sub_male_base)

map_pred_brt(sub_male_preds, sub_male_train, "Sublegal Male Snow Crab") # works better when saved as .jpeg


## Heatmap of relative influence ----
mat_female_abun_inf <- summary(brt_mat_female_abun$model)[-1]
colnames(mat_female_abun_inf)[1] <- "mature female"
rownames(mat_female_abun_inf)[rownames(mat_female_abun_inf) == "female_loading_station"] <- "fishery loading"
rownames(mat_female_abun_inf)[rownames(mat_female_abun_inf) == "bcs_mature_female"] <- "BCS"
mat_female_pres_inf <- summary(brt_mat_female_base$model)[-1]
colnames(mat_female_pres_inf)[1] <- "mature female"
rownames(mat_female_pres_inf)[rownames(mat_female_pres_inf) == "female_loading_station"] <- "fishery loading"
rownames(mat_female_pres_inf)[rownames(mat_female_pres_inf) == "bcs_mature_female"] <- "BCS"
imm_female_abun_inf <- summary(brt_imm_female_abun$model)[-1]
colnames(imm_female_abun_inf)[1] <- "immature female"
rownames(imm_female_abun_inf)[rownames(imm_female_abun_inf) == "female_loading_station"] <- "fishery loading"
rownames(imm_female_abun_inf)[rownames(imm_female_abun_inf) == "bcs_immature_female"] <- "BCS"
imm_female_pres_inf <- summary(brt_imm_female_base$model)[-1]
colnames(imm_female_pres_inf)[1] <- "immature female"
rownames(imm_female_pres_inf)[rownames(imm_female_pres_inf) == "female_loading_station"] <- "fishery loading"
rownames(imm_female_pres_inf)[rownames(imm_female_pres_inf) == "bcs_immature_female"] <- "BCS"
leg_male_abun_inf <- summary(brt_leg_male_abun$model)[-1]
colnames(leg_male_abun_inf)[1] <- "legal male"
rownames(leg_male_abun_inf)[rownames(leg_male_abun_inf) == "legal_male_loading_station"] <- "fishery loading"
rownames(leg_male_abun_inf)[rownames(leg_male_abun_inf) == "bcs_legal_male"] <- "BCS"
leg_male_pres_inf <- summary(brt_leg_male_base$model)[-1]
colnames(leg_male_pres_inf)[1] <- "legal male"
rownames(leg_male_pres_inf)[rownames(leg_male_pres_inf) == "legal_male_loading_station"] <- "fishery loading"
rownames(leg_male_pres_inf)[rownames(leg_male_pres_inf) == "bcs_legal_male"] <- "BCS"
sub_male_abun_inf <- summary(brt_sub_male_abun$model)[-1]
colnames(sub_male_abun_inf)[1] <- "sublegal male"
rownames(sub_male_abun_inf)[rownames(sub_male_abun_inf) == "sublegal_male_loading_station"] <- "fishery loading"
rownames(sub_male_abun_inf)[rownames(sub_male_abun_inf) == "bcs_sublegal_male"] <- "BCS"
sub_male_pres_inf <- summary(brt_sub_male_base$model)[-1]
colnames(sub_male_pres_inf)[1] <- "sublegal male"
rownames(sub_male_pres_inf)[rownames(sub_male_pres_inf) == "sublegal_male_loading_station"] <- "fishery loading"
rownames(sub_male_pres_inf)[rownames(sub_male_pres_inf) == "bcs_sublegal_male"] <- "BCS"

all_abun_list <- list(mat_female_abun_inf,
                      imm_female_abun_inf,
                      leg_male_abun_inf,
                      sub_male_abun_inf)

all_pres_list <- list(mat_female_pres_inf,
                      imm_female_pres_inf,
                      leg_male_pres_inf,
                      sub_male_pres_inf)

# Create figures
rel_inf_fig(all_abun_list)
title("Abundance", line = -1.7)


rel_inf_fig(all_pres_list)
title("Presence/Absence", line = -1.7)