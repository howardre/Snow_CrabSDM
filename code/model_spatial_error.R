# Title: Snow Crab Forecast Models
# Purpose: Reformat to use variables for forecast
# Data created: 05/04/2023

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
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'brt_grid_preds.R'))
source(here('code/functions', 'gbm.plot2.R'))
source(here('code/functions', 'grid_development.R'))
source(here('code/functions', 'grid_search.R'))
source(here('code/functions', 'map_pred_brt.R'))
source(here('code/functions', 'part_depen.R'))
source(here('code/functions', 'rel_inf.R'))

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

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
                pres_mat_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_mat_female) 

mat_female_train_warm <- crab_train_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_mat_female)

imm_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_imm_female)

imm_female_train_warm <- crab_train_warm %>%
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

leg_male_train_warm <- crab_train_warm %>%
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

sub_male_train_warm <- crab_train_warm %>%
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

mat_female_test_warm <- crab_test_warm %>%
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

imm_female_test_warm <- crab_test_warm %>%
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

leg_male_test_warm <- crab_test_warm %>%
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
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

sub_male_test_warm <- crab_test_warm %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year, sublegal_male_loading_station) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

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

# Boosted regression trees ----
# Adjust the bag fraction to a value between 0.5-0.75 as suggested by Elith et al. (2008)
# The learning rate could range from 0.1-0.0001, higher value usually means less trees
# Depending on the number of samples, want tree complexity to be high enough (likely using 5)
# Want at least 1000 trees, but don't need to go way beyond it
vars <- c(1:8, 10, 16)

## Mature females ----
# Get best models
# Read in BRTs
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun_base.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base_pres.rds'))

brt_mat_female_base_warm <- grid_search(mat_female_train_warm, 13, 'bernoulli')
brt_mat_female_base_warm

brt_mat_female_abun_warm <- grid_search(mat_female_train_warm[mat_female_train_warm$lncount_mat_female > 0, ],
                                   11, 'gaussian')
brt_mat_female_abun_warm

# Predict on test data
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
# 1.65 for base model
rmse_mat_female_brt_warm <- sqrt(mean((mat_female_test_warm$lncount_mat_female - mat_female_test_warm$pred_brt)^2))
rmse_mat_female_brt_warm # 1.57

# Calculate deviance explained
dev_mat_female_abun_warm <- brt_deviance(brt_mat_female_abun_warm)
dev_mat_female_pres_warm <- brt_deviance(brt_mat_female_base_warm)

dev_mat_female_abun_warm # 41.1% deviance explained
dev_mat_female_pres_warm # 56.2% deviance explained

# Save models for future use
saveRDS(brt_mat_female_abun_warm, file = here('data', 'brt_mat_female_abun_warm.rds'))
saveRDS(brt_mat_female_base_warm, file = here('data', 'brt_mat_female_pres_warm.rds'))

# Spatial RMSE
# Map RMSE
# Need to get average RMSE per station
# Then plot by station to get spatial error
EBS <- get_base_layers(select.region = 'ebs', set.crs = 'auto')
EBS_grid <- EBS$survey.grid
EBS_poly <- st_cast(EBS_grid, "MULTIPOLYGON")
EBS_trans <- st_transform(EBS_poly, "+proj=longlat +datum=NAD83") # change to lat/lon
mat_female_rmse1 <- mat_female_test1 %>%
  group_by(station) %>%
summarize(rmse = Metrics::rmse(lncount_mat_female, pred_brt)) # calculate by station
mat_female_df1 <- merge(mat_female_rmse1,  
                        EBS_trans, 
                        by.x = "station",
                        by.y = "STATIONID") # make spatial
mat_female_sf1 <- st_as_sf(mat_female_df1, # turn into sf object to plot
                           crs = 4269)

ggplot() +
  geom_sf(data = mat_female_sf1,
          aes(fill = rmse)) +
  scale_x_continuous(name = "Longitude", 
                     breaks = EBS$lon.breaks) + 
  scale_y_continuous(name = "Latitude", 
                     breaks = EBS$lat.breaks)

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

# Relative influence
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
part_depen(brt_sub_male_pres)
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