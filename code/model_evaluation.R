# Title: Snow Crab Model Exploration
# Purpose: Investigate potential models
# Data created: 07/20/2022

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(here)
library(mgcv)
library(dplyr)
library(colorspace)
library(maps)
library(mapdata)
library(fields)
library(ggplot2)
library(enmSdmX) # use for grid search, wrapper for dismo
library(fastshap) # calculate SHAP values quickly
library(mshap) # combine SHAP values for two-part models
library(shapviz) # visualize SHAP values
library(doParallel)
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'brt_grid_preds.R'))
source(here('code/functions', 'calculate_shap.R'))
source(here('code/functions', 'gbm.plot2.R'))
source(here('code/functions', 'grid_development.R'))
source(here('code/functions', 'grid_search.R'))
source(here('code/functions', 'map_pred_brt.R'))
source(here('code/functions', 'part_depen.R'))
source(here('code/functions', 'plot_variable.R'))
source(here('code/functions', 'rel_inf.R'))
source(here('code/functions', 'variable_figure.R'))
source(here('code/functions', 'sv_dependence2D2.R'))
source(here('code/functions', 'sv_dependence2.R'))

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
options(shapviz.brewer_args = list(low = "darkslateblue",
                                   mid = "gainsboro",
                                   high = "darkred",
                                   midpoint = 0)) # use to change color in sv_dependence2D2 function

# Load data ----
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
         !is.na(depth)) 

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
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

# Functions and LOESS ----
# Use LOESS for prediction grid, works better than fixed values for depth and phi
depth_loess <- loess(depth ~ longitude * latitude,
                     data = crab_trans,
                     span = 0.7,
                     degree = 2,
                     control = loess.control(surface = "interpolate"))
summary(lm(depth_loess$fitted ~ crab_trans$depth)) # check R2


phi_loess <- loess(phi ~ longitude * latitude,
                   span = 0.05,
                   degree = 2,
                   data = crab_trans,
                   family = "symmetric",
                   control = loess.control(surface = "interpolate"))
summary(lm(phi_loess$fitted ~ crab_trans$phi)) # check R2

# GAMs ----
# Need to remove year and rerun
## Mature Female ----
# Gaussian
# Base model with presence/absence
mat_female_gam_base <- gam(pres_mat_female ~ s(longitude, latitude) +
                             s(julian),
                           data = mat_female_train,
                           family = "binomial")
summary(mat_female_gam_base) # 51.5% explained

# Abundance model
mat_female_gam_abun <- gam(lncount_mat_female ~ s(longitude, latitude) +
                             s(julian) +
                             s(depth, k = 5) +
                             s(phi, k = 5) +
                             s(temperature, k = 5) +
                             s(ice_mean, k = 5) +
                             s(female_loading_station, k = 5) +
                             s(log_pcod_cpue, k = 5) +
                             s(bcs_mature_female, k = 5),
                           data = mat_female_train[mat_female_train$lncount_mat_female > 0, ])
summary(mat_female_gam_abun) # 32%

par(mfrow = c(2, 2))
gam.check(mat_female_gam_abun)

par(mfrow = c(3, 3))
plot(mat_female_gam_abun)

# Tweedie
mat_female_tweedie <- gam(mature_female + 1 ~ s(longitude, latitude) +
                             s(julian) +
                             s(depth, k = 5) +
                             s(phi, k = 5) +
                             s(temperature, k = 5) +
                             s(ice_mean, k = 5) +
                             s(female_loading_station, k = 5) +
                             s(log_pcod_cpue, k = 5) +
                             s(bcs_mature_female, k = 5),
                           data = mat_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie) # 66%

par(mfrow = c(2, 2))
gam.check(mat_female_tweedie)

par(mfrow = c(3, 3))
plot(mat_female_tweedie)

# Predict on test data
mat_female_test$pred_gam <- predict(mat_female_tweedie,
                                    mat_female_test,
                                    type = "link") # remove to allow predictions on new years
mat_female_test$pred_gam_base <- predict(mat_female_gam_base,
                                         mat_female_test,
                                         type = "response")
mat_female_test$pred_gam_abun <- predict(mat_female_gam_abun,
                                         mat_female_test,
                                         type = "response")

mat_female_test$pred_gam_delta <- mat_female_test$pred_gam_base * mat_female_test$pred_gam_abun

rmse_mat_female_tweedie <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_gam)^2, na.rm = T))
rmse_mat_female_tweedie # 2.65

rmse_mat_female_delta <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_gam_delta)^2, na.rm = T))
rmse_mat_female_delta # 2.58

# Variable plots
tiff(here('results/GAM',
          'gaussian_mat_female_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(mat_female_gam_abun, c(-7.5, 4))
dev.off()

tiff(here('results/GAM',
          'tweedie_mat_female_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(mat_female_tweedie, c(-7.5, 4))
dev.off()

## Immature Female ----
# Gaussian
# Base model with presence/absence
imm_female_gam_base <- gam(pres_imm_female ~ s(longitude, latitude) +
                             s(julian),
                           data = imm_female_train,
                           family = "binomial")
summary(imm_female_gam_base) # 43.3% explained

# Abundance model
imm_female_gam_abun <- gam(lncount_imm_female ~ s(longitude, latitude) +
                             s(julian) +
                             s(depth, k = 5) +
                             s(phi, k = 5) +
                             s(temperature, k = 5) +
                             s(ice_mean, k = 5) +
                             s(female_loading_station, k = 5) +
                             s(log_pcod_cpue, k = 5) +
                             s(bcs_immature_female, k = 5),
                           data = imm_female_train[imm_female_train$lncount_imm_female > 0, ])
summary(imm_female_gam_abun) # 34.2%

par(mfrow = c(2, 2))
gam.check(imm_female_gam_abun)

par(mfrow = c(3, 4))
plot(imm_female_gam_abun)

# Tweedie
imm_female_tweedie <- gam(immature_female + 1 ~ s(longitude, latitude) +
                            s(julian) +
                            s(depth, k = 5) +
                            s(phi, k = 5) +
                            s(temperature, k = 5) +
                            s(ice_mean, k = 5) +
                            s(female_loading_station, k = 5) +
                            s(log_pcod_cpue, k = 5) +
                            s(bcs_immature_female, k = 5),
                          data = imm_female_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(imm_female_tweedie) # 69.2%

par(mfrow = c(2, 2))
gam.check(imm_female_tweedie)

par(mfrow = c(3, 3))
plot(imm_female_tweedie)

# Predict on test data
imm_female_test$pred_gam <- predict(imm_female_tweedie,
                                    imm_female_test,
                                    type = "link") 
imm_female_test$pred_gam_base <- predict(imm_female_gam_base,
                                         imm_female_test,
                                         type = "response")
imm_female_test$pred_gam_abun <- predict(imm_female_gam_abun,
                                         imm_female_test,
                                         type = "response")

imm_female_test$pred_gam_delta <- imm_female_test$pred_gam_base * imm_female_test$pred_gam_abun

rmse_imm_female_tweedie <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_gam)^2, na.rm = T))
rmse_imm_female_tweedie # 2.01

rmse_imm_female_delta <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_gam_delta)^2, na.rm = T))
rmse_imm_female_delta # 1.05

# Variable plots
tiff(here('results/GAM',
          'gaussian_imm_female_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(imm_female_gam_abun, c(-2.5, 2.5))
dev.off()

tiff(here('results/GAM',
          'tweedie_imm_female_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(imm_female_tweedie, c(-2.5, 2.5))
dev.off()

## Legal Male ----
# Gaussian
# Base model with presence/absence
leg_male_gam_base <- gam(pres_leg_male ~ s(longitude, latitude) +
                           s(julian),
                         data = leg_male_train,
                         family = "binomial")
summary(leg_male_gam_base) # 56.1% explained

# Abundance model
leg_male_gam_abun <- gam(lncount_leg_male ~ s(longitude, latitude) +
                           s(julian) +
                           s(depth, k = 5) +
                           s(phi, k = 5) +
                           s(temperature, k = 5) +
                           s(ice_mean, k = 5) +
                           s(legal_male_loading_station, k = 5) +
                           s(log_pcod_cpue, k = 5) +
                           s(bcs_legal_male, k = 5),
                         data = leg_male_train[leg_male_train$lncount_leg_male > 0, ])
summary(leg_male_gam_abun) # 43.0%

par(mfrow = c(2, 2))
gam.check(leg_male_gam_abun)

par(mfrow = c(3, 3))
plot(leg_male_gam_abun)

# Tweedie
leg_male_tweedie <- gam(legal_male + 1 ~ s(longitude, latitude) +
                          s(julian) +
                          s(depth, k = 5) +
                          s(phi, k = 5) +
                          s(temperature, k = 5) +
                          s(ice_mean, k = 5) +
                          s(legal_male_loading_station, k = 5) +
                          s(log_pcod_cpue, k = 5) +
                          s(bcs_legal_male, k = 5),
                        data = leg_male_train,
                        family = tw(link = "log"),
                        method = "REML")
summary(leg_male_tweedie) # 64.1%

par(mfrow = c(2, 2))
gam.check(leg_male_tweedie)

par(mfrow = c(3, 4))
plot(leg_male_tweedie)

# Predict on test data
leg_male_test$pred_gam <- predict(leg_male_tweedie,
                                  leg_male_test,
                                  type = "link") 
leg_male_test$pred_gam_base <- predict(leg_male_gam_base,
                                       leg_male_test,
                                       type = "response")
leg_male_test$pred_gam_abun <- predict(leg_male_gam_abun,
                                       leg_male_test,
                                       type = "response")

leg_male_test$pred_gam_delta <- leg_male_test$pred_gam_base * leg_male_test$pred_gam_abun

rmse_leg_male_tweedie <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_gam)^2, na.rm = T))
rmse_leg_male_tweedie # 1.74

rmse_leg_male_delta <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_gam_delta)^2, na.rm = T))
rmse_leg_male_delta # 1.25

# Variable plots
tiff(here('results/GAM',
          'gaussian_leg_male_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(leg_male_gam_abun, c(-3.5, 4))
dev.off()

tiff(here('results/GAM',
          'tweedie_leg_male_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(leg_male_tweedie, c(-3, 2.5))
dev.off()

## Sublegal Male ----
# Gaussian
# Base model with presence/absence
sub_male_gam_base <- gam(pres_sub_male ~ s(longitude, latitude) +
                           s(julian),
                         data = sub_male_train,
                         family = "binomial")
summary(sub_male_gam_base) # 60.4% explained

# Abundance model
sub_male_gam_abun <- gam(lncount_sub_male ~ s(longitude, latitude) +
                           s(julian) +
                           s(depth, k = 5) +
                           s(phi, k = 5) +
                           s(temperature, k = 5) +
                           s(ice_mean, k = 5) +
                           s(sublegal_male_loading_station, k = 5) +
                           s(log_pcod_cpue, k = 5) +
                           s(bcs_sublegal_male, k = 5),
                         data = sub_male_train[sub_male_train$lncount_sub_male > 0, ])
summary(sub_male_gam_abun) # 53.3%

par(mfrow = c(2, 2))
gam.check(sub_male_gam_abun)

par(mfrow = c(3, 4))
plot(sub_male_gam_abun)

# Tweedie
sub_male_tweedie <- gam(sublegal_male + 1 ~ s(longitude, latitude) +
                          s(julian) +
                          s(depth, k = 5) +
                          s(phi, k = 5) +
                          s(temperature, k = 5) +
                          s(ice_mean, k = 5) +
                          s(sublegal_male_loading_station, k = 5) +
                          s(log_pcod_cpue, k = 5) +
                          s(bcs_sublegal_male, k = 5),
                        data = sub_male_train,
                        family = tw(link = "log"),
                        method = "REML")
summary(sub_male_tweedie) # 69.4%

par(mfrow = c(2, 2))
gam.check(sub_male_tweedie)

par(mfrow = c(3, 4))
plot(sub_male_tweedie)

# Predict on test data
sub_male_test$pred_gam <- predict(sub_male_tweedie,
                                  sub_male_test,
                                  type = "link") 
sub_male_test$pred_gam_base <- predict(sub_male_gam_base,
                                       sub_male_test,
                                       type = "response")
sub_male_test$pred_gam_abun <- predict(sub_male_gam_abun,
                                       sub_male_test,
                                       type = "response")

sub_male_test$pred_gam_delta <- sub_male_test$pred_gam_base * sub_male_test$pred_gam_abun

rmse_sub_male_tweedie <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_gam)^2, na.rm = T))
rmse_sub_male_tweedie # 1.84

rmse_sub_male_delta <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_gam_delta)^2, na.rm = T))
rmse_sub_male_delta # 1.29

# Variable plots
tiff(here('results/GAM',
          'gaussian_sub_male_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(sub_male_gam_abun, c(-7, 4.5))
dev.off()

tiff(here('results/GAM',
          'tweedie_sub_male_plots.jpeg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(sub_male_tweedie, c(-6.5, 2))
dev.off()


# Boosted regression trees ----
# Adjust the bag fraction to a value between 0.5-0.75 as suggested by Elith et al. (2008)
# The learning rate could range from 0.1-0.0001, higher value usually means less trees
# Depending on the number of samples, want tree complexity to be high enough (likely using 5)
# Want at least 1000 trees, but don't need to go way beyond it
vars <- c(1:8, 10, 16)

## Mature females ----
# Get best models
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
rmse_mat_female_brt # 1.62

# Calculate deviance explained
dev_mat_female_abun <- brt_deviance(brt_mat_female_abun)
dev_mat_female_pres <- brt_deviance(brt_mat_female_base)

dev_mat_female_abun # 33.5% deviance explained
dev_mat_female_pres # 56.6% deviance explained

# Save models for future use
saveRDS(brt_mat_female_abun, file = here('data', 'brt_mat_female_abun.rds'))
saveRDS(brt_mat_female_base, file = here('data', 'brt_mat_female_base.rds'))

# Read in BRTs
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base.rds'))

# Variable names
column_names <- c("depth", "temperature", "phi", "ice_mean", "bcs_mature_female", "longitude",
                  "latitude", "julian", "female_loading_station", "log_pcod_cpue")
final_names <- c("depth", "temperature", "phi", "ice concentration",
                 "proportion BCS", "longitude", "latitude", "julian",
                 "female loading", "log(cod cpue + 1)")

column_labels <- data.frame(column_names, final_names)
match_labels <- column_labels[match(colnames(mat_female_train)[vars], column_labels$column_names), ]

brt_mat_female_summary <- summary(brt_mat_female_abun$model)

brt_labels <- match_labels[order(match(names(mat_female_train)[vars], brt_mat_female_summary$var)), ]
labels <- brt_labels$final_names

# Variable importance
rel_inf(brt_mat_female_base, 'Variable Influence on Mature Female Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'female_mat_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_mat_female_base)
dev.copy(jpeg,
         here('results/BRT',
              'female_mat_plots.jpg'),
         height = 12,
         width = 15,
         res = 200,
         units = 'in')
dev.off()


# Plot the fits
females_mat_int <- gbm.interactions(brt_mat_female_abun$model)
females_mat_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(brt_mat_female_abun$model,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_mat_female_abun$model,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_mat_female_abun$model,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Plot map of the predicted distribution
# Prediction grid map
spatial_grid_mat_female <- grid_development(mat_female_train)
spatial_grid_mat_female$female_loading_station <- median(mat_female_train$female_loading_station, na.rm = T)
spatial_grid_mat_female$bcs_mature_female <- median(mat_female_train$bcs_mature_female, na.rm = T)
mat_female_preds <- brt_grid_preds(spatial_grid_mat_female, 
                                   brt_mat_female_abun,
                                   brt_mat_female_base)

map_pred_brt(mat_female_preds, mat_female_train, "Distribution of Mature Female Snow Crab (BRT)")
dev.copy(jpeg,
         here('results/BRT',
              'female_mat_map.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

## Immature females ----
# Get best models
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
rmse_imm_female_brt # 1.40

# Calculate deviance
dev_imm_female_abun <- brt_deviance(brt_imm_female_abun)
dev_imm_female_pres <- brt_deviance(brt_imm_female_base)

dev_imm_female_abun # 50.6% deviance explained
dev_imm_female_pres # 47.1% deviance explained

# Save models for future use
saveRDS(brt_imm_female_abun, file = here('data', 'brt_imm_female_abun.rds'))
saveRDS(brt_imm_female_base, file = here('data', 'brt_imm_female_base.rds'))

# Read in BRTs
brt_imm_female_abun <- readRDS(file = here('data', 'brt_imm_female_abun.rds'))
brt_imm_female_base <- readRDS(file = here('data', 'brt_imm_female_base.rds'))

# Variable importance
rel_inf(brt_imm_female_base, 'Variable Influence on Immature Female Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'female_imm_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_imm_female_base)
dev.copy(jpeg,
         here('results/BRT',
              'female_imm_plots.jpg'),
         height = 12,
         width = 15,
         res = 200,
         units = 'in')
dev.off()


# Plot the fits
females_imm_int <- gbm.interactions(brt_imm_female_abun$model)
females_imm_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(brt_imm_female_abun$model,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_imm_female_abun$model,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_imm_female_abun$model,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Plot map of the predicted distribution
# Prediction grid map
# Variable names
match_labels <- column_labels[match(colnames(imm_female_train)[vars], column_labels$column_names), ]

brt_imm_female_summary <- summary(brt_imm_female_abun$model)

brt_labels <- match_labels[order(match(names(imm_female_train)[vars], brt_imm_female_summary$var)), ]
labels <- brt_labels$final_names

column_names <- c("depth", "temperature", "phi", "ice_mean", "bcs_immature_female", "longitude",
                  "latitude", "julian", "female_loading_station", "log_pcod_cpue", "year_f")
final_names <- c("depth", "temperature", "phi", "ice concentration",
                 "proportion BCS", "longitude", "latitude", "julian",
                 "female loading", "log(cod cpue + 1)", "year")
column_labels <- data.frame(column_names, final_names)

spatial_grid_imm_female <- grid_development(imm_female_train)
spatial_grid_imm_female$female_loading_station <- median(imm_female_train$female_loading_station, na.rm = TRUE)
spatial_grid_imm_female$bcs_immature_female <- median(imm_female_train$bcs_immature_female, na.rm = TRUE)

imm_female_preds <- brt_grid_preds(spatial_grid_imm_female, 
                                   brt_imm_female_abun,
                                   brt_imm_female_base)

map_pred_brt(imm_female_preds, imm_female_train, "Distribution of Immature Female Snow Crab (BRT)")
dev.copy(jpeg,
         here('results/BRT',
              'female_imm_map.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()


##Legal Males ----
# Get best models
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
rmse_leg_male_brt # 1.22

# Calculate deviance
dev_leg_male_abun <- brt_deviance(brt_leg_male_abun)
dev_leg_male_pres <- brt_deviance(brt_leg_male_base)

dev_leg_male_abun # 64% deviance explained
dev_leg_male_pres # 62% deviance explained

# Save models for future use
saveRDS(brt_leg_male_abun, file = here('data', 'brt_leg_male_abun.rds'))
saveRDS(brt_leg_male_base, file = here('data', 'brt_leg_male_base.rds'))

# Read in BRTs
brt_leg_male_abun <- readRDS(file = here('data', 'brt_leg_male_abun.rds'))
brt_leg_male_base <- readRDS(file = here('data', 'brt_leg_male_base.rds'))

# Variable importance
rel_inf(brt_leg_male_base, 'Variable Influence on Legal Male Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'male_leg_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_leg_male_base)
dev.copy(jpeg,
         here('results/BRT',
              'male_leg_plots.jpg'),
         height = 12,
         width = 15,
         res = 200,
         units = 'in')
dev.off()

# Plot the fits
males_leg_int <- gbm.interactions(brt_leg_male_abun$model)
males_leg_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(brt_leg_male_abun$model,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_leg_male_abun$model,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_leg_male_abun$model,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Plot map of the predicted distribution
# Prediction grid map
# Variable names
match_labels <- column_labels[match(colnames(leg_male_train)[vars], column_labels$column_names), ]

brt_leg_male_summary <- summary(brt_leg_male_abun$model)

brt_labels <- match_labels[order(match(names(leg_male_train)[vars], brt_leg_male_summary$var)), ]
labels <- brt_labels$final_names

column_names <- c("depth", "temperature", "phi", "ice_mean", "bcs_legal_male", "longitude",
                  "latitude", "julian", "legal_male_loading", "log_pcod_cpue", "year_f")
final_names <- c("depth", "temperature", "phi", "ice concentration",
                 "proportion BCS", "longitude", "latitude", "julian",
                 "legal male loading", "log(cod cpue + 1)", "year")
column_labels <- data.frame(column_names, final_names)

spatial_grid_leg_male <- grid_development(leg_male_train)
spatial_grid_leg_male$legal_male_loading_station <- median(leg_male_train$legal_male_loading_station, na.rm = TRUE)
spatial_grid_leg_male$bcs_legal_male <- median(leg_male_train$bcs_legal_male, na.rm = TRUE)

leg_male_preds <- brt_grid_preds(spatial_grid_leg_male,
                                 brt_leg_male_abun,
                                 brt_leg_male_base)

map_pred_brt(leg_male_preds, leg_male_train, "Distribution of Legal Male Snow Crab (BRT)")
dev.copy(jpeg,
         here('results/BRT',
              'male_leg_map.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

## Sublegal Males ----
# Get best models
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
rmse_sub_male_brt # 1.57

# Calculate deviance explained
dev_sub_male_abun <- brt_deviance(brt_sub_male_abun)
dev_sub_male_pres <- brt_deviance(brt_sub_male_base)

dev_sub_male_abun # 65.3%% deviance explained
dev_sub_male_pres # 62.1% deviance explained

# Save models for future use
saveRDS(brt_sub_male_abun, file = here('data', 'brt_sub_male_abun.rds'))
saveRDS(brt_sub_male_base, file = here('data', 'brt_sub_male_base.rds'))

# Read in BRTs
brt_sub_male_abun <- readRDS(file = here('data', 'brt_sub_male_abun.rds'))
brt_sub_male_base <- readRDS(file = here('data', 'brt_sub_male_base.rds'))

# Variable names
match_labels <- column_labels[match(colnames(sub_male_train)[vars], column_labels$column_names), ]

brt_sub_male_summary <- summary(brt_sub_male_abun$model)

brt_labels <- match_labels[order(match(names(sub_male_train)[vars], brt_sub_male_summary$var)), ]
labels <- brt_labels$final_names

column_names <- c("depth", "temperature", "phi", "ice_mean", "bcs_sublegal_male", "longitude",
                  "latitude", "julian", "sublegal_male_loading", "log_pcod_cpue", "year_f")
final_names <- c("depth", "temperature", "phi", "ice concentration",
                 "proportion BCS", "longitude", "latitude", "julian",
                 "sublegal male loading", "log(cod cpue + 1)", "year")
column_labels <- data.frame(column_names, final_names)

# Variable importance
rel_inf(brt_sub_male_abun, 'Variable Influence on Sublegal Male Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'male_sub_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_sub_male_abun)
dev.copy(jpeg,
         here('results/BRT',
              'male_sub_plots.jpg'),
         height = 12,
         width = 15,
         res = 200,
         units = 'in')
dev.off()

# Plot the fits
males_sub_int <- gbm.interactions(brt_sub_male_abun$model)
males_sub_int$interactions

windows()
par(mfrow = c(1, 3))
gbm.perspec(brt_sub_male_abun$model,
            2, 7,
            z.range = c(0, 8.3),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_sub_male_abun$model,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(brt_sub_male_abun$model,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Plot map of the predicted distribution
# Prediction grid map
spatial_grid_sub_male <- grid_development(sub_male_train)
spatial_grid_sub_male$sublegal_male_loading <- median(sub_male_train$sublegal_male_loading, na.rm = TRUE)
spatial_grid_sub_male$bcs_sublegal_male <- median(sub_male_train$bcs_sublegal_male, na.rm = TRUE)

sub_male_preds <- brt_grid_preds(spatial_grid_sub_male, 
                                 brt_sub_male_abun,
                                 brt_sub_male_base)

map_pred_brt(sub_male_preds, sub_male_train, "Distribution of Sublegal Male Snow Crab (BRT)")
dev.copy(jpeg,
         here('results/BRT',
              'male_sub_map.jpg'),
         height = 10,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# SHAP values ----
# SHapley Additive exPlanations
leg_male_names <- c("depth", "temperature", "phi", "ice_mean", 
                    "julian", "legal_male_loading_station",
                    "bcs_legal_male", "log_pcod_cpue")
sub_male_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                    "latitude", "julian", "sublegal_male_loading_station",
                    "bcs_sublegal_male", "log_pcod_cpue")
mat_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", "female_loading_station",
                      "bcs_mature_female", "log_pcod_cpue")
imm_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", "female_loading_station",
                      "bcs_immature_female", "log_pcod_cpue")
pred_fun <- function(X_model, newdata){
  gbm::predict.gbm(X_model, 
                   newdata,
                   n.trees = X_model$model$gbm.call$best.trees,
                   type = "response")
}

num_cores <- detectCores() - 2

## Legal Males ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
leg_male_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male, ice_mean) # full data set with just important variables

leg_male_shaps <- calculate_shap(brt_leg_male_abun, brt_leg_male_base, leg_male_data)
stopCluster(cl)

saveRDS(leg_male_shaps, file = here('data', 'leg_male_shaps.rds'))

# Visualize
# Gives the global effect of variables (absolute value, not directional)
leg_male_shaps <- readRDS(file = here('data', 'leg_male_shaps.rds'))

leg_male_mshap <- as.matrix(leg_male_shaps[[3]]$shap_vals)
leg_male_mshap_sv <- shapviz(leg_male_mshap, X = leg_male_data)
leg_male_pshap <- as.matrix(leg_male_shaps[[2]]$shapley_values)
leg_male_pshap_sv <- shapviz(leg_male_pshap, X = leg_male_data)
leg_male_ashap <- as.matrix(leg_male_shaps[[1]]$shapley_values)
leg_male_ashap_sv <- shapviz(leg_male_ashap, X = leg_male_data)

# Variable importance
sv_importance(leg_male_mshap_sv)
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(leg_male_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(leg_male_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(leg_male_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(leg_male_mshap_sv, 
              v = leg_male_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
sv_dependence2D2(leg_male_mshap_sv, 
                x = "longitude", 
                y = "latitude",
                size = 2.5,
                jitter_width = 0.5,
                jitter_height = 0.5) +  
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Spatial SHAP Values for Legal Male Crab")
dev.copy(jpeg,
         here('results/BRT',
              'leg_male_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Force plots 
# Yellow means variable pushes prediction higher, purple means variable pushes prediction lower
# Scores close to the f(x) value have more of an impact (indicated by SHAP magnitude too)
# Would like to figure out how to separate by year without including in the model
sv_force(leg_male_sv, 2) # one observation

## Sublegal Males ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
sub_male_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station) %>%
  tidyr::drop_na(lncount_sub_male, ice_mean) 

sub_male_shaps <- calculate_shap(brt_sub_male_abun, brt_sub_male_base, sub_male_data)
stopCluster(cl)

saveRDS(sub_male_shaps, file = here('data', 'sub_male_shaps.rds'))

# Visualize
sub_male_shaps <- readRDS(file = here('data', 'sub_male_shaps.rds'))

sub_male_mshap <- as.matrix(sub_male_shaps[[3]]$shap_vals)
sub_male_mshap_sv <- shapviz(sub_male_mshap, X = sub_male_data)
sub_male_pshap <- as.matrix(sub_male_shaps[[2]]$shapley_values)
sub_male_pshap_sv <- shapviz(sub_male_pshap, X = sub_male_data)
sub_male_ashap <- as.matrix(sub_male_shaps[[1]]$shapley_values)
sub_male_ashap_sv <- shapviz(sub_male_ashap, X = sub_male_data)

# Variable importance
sv_importance(sub_male_mshap_sv)
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(sub_male_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(sub_male_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(sub_male_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(sub_male_mshap_sv, 
              v = sub_male_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
sv_dependence2D2(sub_male_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Spatial SHAP Values for Sublegal Male Crab")
dev.copy(jpeg,
         here('results/BRT',
              'sub_male_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## Mature Females ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
mat_female_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_mat_female, ice_mean) 

mat_female_shaps <- calculate_shap(brt_mat_female_abun, brt_mat_female_base, mat_female_data)
stopCluster(cl)

saveRDS(mat_female_shaps, file = here('data', 'mat_female_shaps.rds'))

# Visualize
mat_female_shaps <- readRDS(file = here('data', 'mat_female_shaps.rds'))

mat_female_mshap <- as.matrix(mat_female_shaps[[3]]$shap_vals)
mat_female_mshap_sv <- shapviz(mat_female_mshap, X = mat_female_data)
mat_female_pshap <- as.matrix(mat_female_shaps[[2]]$shapley_values)
mat_female_pshap_sv <- shapviz(mat_female_pshap, X = mat_female_data)
mat_female_ashap <- as.matrix(mat_female_shaps[[1]]$shapley_values)
mat_female_ashap_sv <- shapviz(mat_female_ashap, X = mat_female_data)

# Variable importance
sv_importance(mat_female_mshap_sv)
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(mat_female_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(mat_female_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(mat_female_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(mat_female_mshap_sv, 
              v = mat_female_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
sv_dependence2D2(mat_female_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Spatial SHAP Values for Mature Female Crab")
dev.copy(jpeg,
         here('results/BRT',
              'mat_female_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## Immature Females ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
imm_female_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_imm_female, ice_mean) 

imm_female_shaps <- calculate_shap(brt_imm_female_abun, brt_imm_female_base, imm_female_data)
stopCluster(cl)

saveRDS(imm_female_shaps, file = here('data', 'imm_female_shaps.rds'))

# Visualize
imm_female_shaps <- readRDS(file = here('data', 'imm_female_shaps.rds'))

imm_female_mshap <- as.matrix(imm_female_shaps[[3]]$shap_vals)
imm_female_mshap_sv <- shapviz(imm_female_mshap, X = imm_female_data)
imm_female_pshap <- as.matrix(imm_female_shaps[[2]]$shapley_values)
imm_female_pshap_sv <- shapviz(imm_female_pshap, X = imm_female_data)
imm_female_ashap <- as.matrix(imm_female_shaps[[1]]$shapley_values)
imm_female_ashap_sv <- shapviz(imm_female_ashap, X = imm_female_data)

# Variable importance
sv_importance(imm_female_mshap_sv)
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(imm_female_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(imm_female_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(imm_female_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black")
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(imm_female_mshap_sv, 
              v = imm_female_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
sv_dependence2D2(imm_female_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude", 
       y = "Latitude",
       title = "Spatial SHAP Values for Immature Female Crab")
dev.copy(jpeg,
         here('results/BRT',
              'imm_female_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()