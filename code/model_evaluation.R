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
library(kernelshap)
library(shapviz)
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))

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
                pres_mat_female, year_f, year, female_loading_station) %>%
  tidyr::drop_na(lncount_mat_female) 

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

# Plot GAM output
plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam,
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = "#eb8055b2",
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6.5,
       cex.lab = 6.5,
       family = "serif",
       lwd = 2.5)
}

variable_figure <- function(gam, range){
  par(mfrow = c(3, 3),
      mar = c(17, 21, .5, 0.6) + 0.1,
      oma = c(3, 1, 1, 1),
      mgp =  c(9, 3, 0))
  plot_variable(gam,
                covariate = 3,
                bounds = range,
                "Day of Year",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 4,
                bounds = range,
                "Depth",
                " ",
                "n")
  plot_variable(gam,
                covariate = 5,
                bounds = range,
                "Phi",
                " ",
                "n")
  plot_variable(gam,
                covariate = 6,
                bounds = range,
                "Temperature",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 7,
                bounds = range,
                "Ice Concentration",
                " ",
                "n")
  plot_variable(gam,
                covariate = 8,
                bounds = range,
                "PCA Loading",
                " ",
                "n")
  plot_variable(gam,
                covariate = 9,
                bounds = range,
                "Cod Abundance",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 10,
                bounds = range,
                "BCS Prevalence",
                " ",
                "n")
  }

# Rewrite gbm.plot to use better variable names
gbm.plot2 <- function (gbm.object, variable.no = 0, smooth = FALSE, rug = TRUE, 
                       n.plots = length(pred.names), common.scale = TRUE, write.title = TRUE, 
                       y.label = "fitted function", x.label = NULL, show.contrib = TRUE, 
                       plot.layout = c(3, 4), ...) 
{
  if (!requireNamespace("gbm")) {
    stop("you need to install the gbm package to run this function")
  }
  requireNamespace("splines")
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  pred.names <- gbm.call$predictor.names
  response.name <- gbm.call$response.name
  data <- gbm.call$dataframe
  max.plots <- plot.layout[1] * plot.layout[2]
  plot.count <- 0
  n.pages <- 1
  if (length(variable.no) > 1) {
    stop("only one response variable can be plotted at a time")
  }
  if (variable.no > 0) {
    n.plots <- 1
  }
  max.vars <- length(gbm.object$contributions$var)
  if (n.plots > max.vars) {
    n.plots <- max.vars
    warning("reducing no of plotted predictors to maximum available (", 
            max.vars, ")")
  }
  predictors <- list(rep(NA, n.plots))
  responses <- list(rep(NA, n.plots))
  for (j in c(1:n.plots)) {
    if (n.plots == 1) {
      k <- variable.no
    }
    else {
      k <- match(gbm.object$contributions$var[j], pred.names)
    }
    if (is.null(x.label)) {
      var.name <- labels[j]
    }
    else {
      var.name <- labels[j]
    }
    pred.data <- data[, gbm.call$gbm.x[k]]
    response.matrix <- gbm::plot.gbm(gbm.object, k, return.grid = TRUE)
    predictors[[j]] <- response.matrix[, 1]
    if (is.factor(data[, gbm.call$gbm.x[k]])) {
      predictors[[j]] <- factor(predictors[[j]], levels = levels(data[, 
                                                                      gbm.call$gbm.x[k]]))
    }
    responses[[j]] <- response.matrix[, 2] - mean(response.matrix[, 
                                                                  2])
    if (j == 1) {
      ymin = min(responses[[j]])
      ymax = max(responses[[j]])
    }
    else {
      ymin = min(ymin, min(responses[[j]]))
      ymax = max(ymax, max(responses[[j]]))
    }
  }
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = plot.layout)
  for (j in c(1:n.plots)) {
    if (plot.count == max.plots) {
      plot.count = 0
      n.pages <- n.pages + 1
    }
    plot.count <- plot.count + 1
    if (n.plots == 1) {
      k <- match(pred.names[variable.no], gbm.object$contributions$var)
      if (show.contrib) {
        x.label <- paste(var.name, "  (", round(gbm.object$contributions[k, 
                                                                         2], 1), "%)", sep = "")
      }
    }
    else {
      k <- match(gbm.object$contributions$var[j], pred.names)
      var.name <- labels[j]
      if (show.contrib) {
        x.label <- paste(var.name, "  (", round(gbm.object$contributions[j, 
                                                                         2], 1), "%)", sep = "")
      }
      else x.label <- var.name
    }
    if (common.scale) {
      plot(predictors[[j]], responses[[j]], ylim = c(ymin, 
                                                     ymax), type = "l", xlab = x.label, ylab = y.label, 
           ...)
    }
    else {
      plot(predictors[[j]], responses[[j]], type = "l", 
           xlab = x.label, ylab = y.label, ...)
    }
    if (smooth & is.vector(predictors[[j]])) {
      temp.lo <- loess(responses[[j]] ~ predictors[[j]], 
                       span = 0.3)
      lines(predictors[[j]], fitted(temp.lo), lty = 2, 
            col = 2)
    }
    if (plot.count == 1 & n.plots == 1) {
      if (write.title) {
        title(paste(response.name, " - page ", n.pages, 
                    sep = ""))
      }
      if (rug & is.vector(data[, gbm.call$gbm.x[variable.no]])) {
        rug(quantile(data[, gbm.call$gbm.x[variable.no]], 
                     probs = seq(0, 1, 0.1), na.rm = TRUE))
      }
    }
    else {
      if (write.title & j == 1) {
        title(response.name)
      }
      if (rug & is.vector(data[, gbm.call$gbm.x[k]])) {
        rug(quantile(data[, gbm.call$gbm.x[k]], probs = seq(0, 
                                                            1, 0.1), na.rm = TRUE))
      }
    }
  }
  graphics::par(op)
}


grid_search <- function(data, response, family){
  trainBRT(data = data,
           preds = vars,
           resp = response,
           family = family,
           treeComplexity = c(1, 5, 10),
           learningRate = c(0.01, 0.05, 0.1),
           bagFraction = c(0.25, 0.5, 0.75),
           minTrees = 1000, # recommended minimum by Elith
           maxTrees = 2000,
           cores = 6, # increase speed
           out = c('model', 'tuning')) # should return model and table with hyperparameters
}

rel_inf <- function(abun_brt, title){
  windows()
  effects <- tibble::as_tibble(summary.gbm(abun_brt$model, plotit = FALSE))
  effects %>% arrange(desc(rel.inf)) %>%
    ggplot(aes(x = forcats::fct_reorder(.f = var,
                                        .x = rel.inf),
               y = rel.inf,
               fill = rel.inf)) +
    geom_col() +
    coord_flip() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = 'Variable',
         y = 'Relative Influence',
         title = title) +
    theme_minimal() +
    theme(legend.position = "none",       
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 16),
          axis.title = element_text(family = "serif", size = 20),
          strip.text = element_text(family = "serif", size = 20)) 
  }

part_depen <- function(abun_brt){
  windows()
  gbm.plot(abun_brt$model,
            plot.layout = c(3, 4),
            write.title = F,
            smooth = T,
            common.scale = T,
            cex.axis = 1.7,
            cex.lab = 1.7,
            lwd = 1.5,
            show.contrib = FALSE,
            family = "serif")
}

grid_development <- function(train_data){
  nlat = 40
  nlon = 60
  latd = seq(min(train_data$latitude), max(train_data$latitude), length.out = nlat)
  lond = seq(min(train_data$longitude), max(train_data$longitude), length.out = nlon)
  spatial_grid <- expand.grid(lond, latd) # create grid
  names(spatial_grid) <- c('longitude', 'latitude')
  spatial_grid$dist <- NA # calculate distance from nearest station
  for (k in 1:nrow(spatial_grid)) {
    dist <-  distance_function(spatial_grid$latitude[k],
                               spatial_grid$longitude[k],
                               train_data$latitude,
                               train_data$longitude)
    spatial_grid$dist[k] <- min(dist)
  }
  
  # Add in LOESS interpolated values for phi and depth
  spatial_grid$depth <- as.vector(predict(depth_loess,
                                          newdata = spatial_grid))
  spatial_grid$phi <- as.vector(predict(phi_loess,
                                        newdata = spatial_grid))
  
  # Add median values for all other variables
  spatial_grid$year_f <- as.factor('2010')
  spatial_grid$julian <- median(train_data$julian, na.rm = T)
  spatial_grid$temperature <- median(train_data$temperature, na.rm = T)
  spatial_grid$ice_mean <- median(train_data$ice_mean, na.rm = T)
  spatial_grid$log_pcod_cpue <- median(train_data$log_pcod_cpue, na.rm = T)
  return(spatial_grid)
}

map_pred_brt <- function(spatial_grid, train_data, title){
  # Create palette dependent on scale
  my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  color_levels = 100
  max_absolute_value = max(abs(c(min(spatial_grid$pred_brt, na.rm = T),
                                 max(spatial_grid$pred_brt, na.rm = T))))
  color_sequence = seq(max(spatial_grid$pred_brt, na.rm = T), 
                       min(spatial_grid$pred_brt, na.rm = T),
                       length.out = color_levels + 1)
  n_in_class = hist(spatial_grid$pred_brt, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  
  # Use for matrix
  nlat = 40
  nlon = 60
  latd = seq(min(train_data$latitude), max(train_data$latitude), length.out = nlat)
  lond = seq(min(train_data$longitude), max(train_data$longitude), length.out = nlon)
  
  # Make map
  windows(width = 12, height = 10)
  par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0),
      family = "serif")
  image(lond,
        latd,
        t(matrix(spatial_grid$pred_brt,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-181, -156),
        ylim = range(train_data$latitude, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(spatial_grid$pred_brt,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include],
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-181, -156),
        ylim = range(train_data$latitude, na.rm = TRUE) + c(-.4, .5),
        main = title,
        cex.main = 2,
        cex.lab = 2,
        cex.axis = 1.8)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.18, .21, .17, .38),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.6,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(spatial_grid$pred_brt, na.rm = T), 
                      max(spatial_grid$pred_brt, na.rm = T)),
             legend.args = list("log(count+1)",
                                side = 2,
                                cex = 1.8,
                                line = 1.3,
                                family =  "serif"))
  }

brt_grid_preds <- function(spatial_grid, abun_brt, base_brt){
  spatial_grid$pred_abun <- predict.gbm(abun_brt$model,
                                        spatial_grid,
                                        n.trees = abun_brt$model$gbm.call$best.trees,
                                        type =  "response")
  spatial_grid$pred_base <- predict.gbm(base_brt$model,
                                        spatial_grid,
                                        n.trees = base_brt$model$gbm.call$best.trees,
                                        type = "response")
  spatial_grid$pred_brt <- spatial_grid$pred_base * spatial_grid$pred_abun
  
  spatial_grid$pred_brt[spatial_grid$dist > 28000] <- NA # Remove predictions too far from samples
  return(spatial_grid)
}

brt_deviance <- function(brt){
  deviance <- (brt$model$self.statistics$mean.null - brt$model$cv.statistics$deviance.mean) /
    brt$model$self.statistics$mean.null
  return(deviance)
}

# GAMs ----
# Must use random effect for year in order to use train/test data set
# Unsure why not needed for BRTs

## Mature Female ----
# Gaussian
# Base model with presence/absence
mat_female_gam_base <- gam(pres_mat_female ~ s(year, bs = "re") +
                             s(longitude, latitude) +
                             s(julian),
                           data = mat_female_train,
                           family = "binomial")
summary(mat_female_gam_base) # 51.2% explained

# Abundance model
mat_female_gam_abun <- gam(lncount_mat_female ~ s(year, bs = "re") +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth, k = 5) +
                             s(phi, k = 5) +
                             s(temperature, k = 5) +
                             s(ice_mean, k = 5) +
                             s(female_loading_station, k = 5) +
                             s(log_pcod_cpue, k = 5) +
                             s(bcs_mature_female, k = 5),
                           data = mat_female_train[mat_female_train$lncount_mat_female > 0, ])
summary(mat_female_gam_abun) # 33.4%

par(mfrow = c(2, 2))
gam.check(mat_female_gam_abun)

par(mfrow = c(3, 4))
plot(mat_female_gam_abun)

# Tweedie
mat_female_tweedie <- gam(mature_female + 1 ~ s(year, bs = 're') +
                             s(longitude, latitude) +
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
summary(mat_female_tweedie) # 61.5%

par(mfrow = c(2, 2))
gam.check(mat_female_tweedie)

par(mfrow = c(3, 3))
plot(mat_female_tweedie)

# Predict on test data
mat_female_test$pred_gam <- predict(mat_female_tweedie,
                                    mat_female_test,
                                    type = "link",
                                    exclude = "s(year)") # remove to allow predictions on new years
mat_female_test$pred_gam_base <- predict(mat_female_gam_base,
                                         mat_female_test,
                                         type = "response",
                                         exclude = "s(year)")
mat_female_test$pred_gam_abun <- predict(mat_female_gam_abun,
                                         mat_female_test,
                                         type = "response",
                                         exclude = "s(year)")

mat_female_test$pred_gam_delta <- mat_female_test$pred_gam_base * mat_female_test$pred_gam_abun

rmse_mat_female_tweedie <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_gam)^2, na.rm = T))
rmse_mat_female_tweedie # 2.80

rmse_mat_female_delta <- sqrt(mean((mat_female_test$lncount_mat_female - mat_female_test$pred_gam_delta)^2, na.rm = T))
rmse_mat_female_delta # 3.59

# Variable plots
tiff(here('results/GAM',
          'gaussian_mat_female_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(mat_female_gam_abun, c(-5, 4))
dev.off()

tiff(here('results/GAM',
          'tweedie_mat_female_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(mat_female_tweedie, c(-6, 4))
dev.off()

## Immature Female ----
# Gaussian
# Base model with presence/absence
imm_female_gam_base <- gam(pres_imm_female ~ s(year, bs = "re") +
                             s(longitude, latitude) +
                             s(julian),
                           data = imm_female_train,
                           family = "binomial")
summary(imm_female_gam_base) # 44.1% explained

# Abundance model
imm_female_gam_abun <- gam(lncount_imm_female ~ s(year, bs = "re") +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth, k = 5) +
                             s(phi, k = 5) +
                             s(temperature, k = 5) +
                             s(ice_mean, k = 5) +
                             s(female_loading_station, k = 5) +
                             s(log_pcod_cpue, k = 5) +
                             s(bcs_immature_female, k = 5),
                           data = imm_female_train[imm_female_train$lncount_imm_female > 0, ])
summary(imm_female_gam_abun) # 48.1%

par(mfrow = c(2, 2))
gam.check(imm_female_gam_abun)

par(mfrow = c(3, 4))
plot(imm_female_gam_abun)

# Tweedie
imm_female_tweedie <- gam(immature_female + 1 ~ s(year, bs = 're') +
                            s(longitude, latitude) +
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
summary(imm_female_tweedie) # 69.7%

par(mfrow = c(2, 2))
gam.check(imm_female_tweedie)

par(mfrow = c(3, 3))
plot(imm_female_tweedie)

# Predict on test data
imm_female_test$pred_gam <- predict(imm_female_tweedie,
                                    imm_female_test,
                                    type = "link",
                                    exclude = "s(year)") 
imm_female_test$pred_gam_base <- predict(imm_female_gam_base,
                                         imm_female_test,
                                         type = "response",
                                         exclude = "s(year)")
imm_female_test$pred_gam_abun <- predict(imm_female_gam_abun,
                                         imm_female_test,
                                         type = "response",
                                         exclude = "s(year)")

imm_female_test$pred_gam_delta <- imm_female_test$pred_gam_base * imm_female_test$pred_gam_abun

rmse_imm_female_tweedie <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_gam)^2, na.rm = T))
rmse_imm_female_tweedie # 2.43

rmse_imm_female_delta <- sqrt(mean((imm_female_test$lncount_imm_female - imm_female_test$pred_gam_delta)^2, na.rm = T))
rmse_imm_female_delta # 2.59

# Variable plots
tiff(here('results/GAM',
          'gaussian_imm_female_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(imm_female_gam_abun, c(-4.5, 3.5))
dev.off()

tiff(here('results/GAM',
          'tweedie_imm_female_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(imm_female_tweedie, c(-5, 2.5))
dev.off()

## Legal Male ----
# Gaussian
# Base model with presence/absence
leg_male_gam_base <- gam(pres_leg_male ~ s(year, bs = "re") +
                           s(longitude, latitude) +
                           s(julian),
                         data = leg_male_train,
                         family = "binomial")
summary(leg_male_gam_base) # 56% explained

# Abundance model
leg_male_gam_abun <- gam(lncount_leg_male ~ s(year, bs = "re") +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth, k = 5) +
                           s(phi, k = 5) +
                           s(temperature, k = 5) +
                           s(ice_mean, k = 5) +
                           s(legal_male_loading_station, k = 5) +
                           s(log_pcod_cpue, k = 5) +
                           s(bcs_legal_male, k = 5),
                         data = leg_male_train[leg_male_train$lncount_leg_male > 0, ])
summary(leg_male_gam_abun) # 41.8%

par(mfrow = c(2, 2))
gam.check(leg_male_gam_abun)

par(mfrow = c(3, 3))
plot(leg_male_gam_abun)

# Tweedie
leg_male_tweedie <- gam(legal_male + 1 ~ s(year, bs = 're') +
                          s(longitude, latitude) +
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
summary(leg_male_tweedie) # 61.5%

par(mfrow = c(2, 2))
gam.check(leg_male_tweedie)

par(mfrow = c(3, 4))
plot(leg_male_tweedie)

# Predict on test data
leg_male_test$pred_gam <- predict(leg_male_tweedie,
                                  leg_male_test,
                                  type = "link",
                                  exclude = "s(year)") 
leg_male_test$pred_gam_base <- predict(leg_male_gam_base,
                                       leg_male_test,
                                       type = "response",
                                       exclude = "s(year)")
leg_male_test$pred_gam_abun <- predict(leg_male_gam_abun,
                                       leg_male_test,
                                       type = "response",
                                       exclude = "s(year)")

leg_male_test$pred_gam_delta <- leg_male_test$pred_gam_base * leg_male_test$pred_gam_abun

rmse_leg_male_tweedie <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_gam)^2, na.rm = T))
rmse_leg_male_tweedie # 1.74

rmse_leg_male_delta <- sqrt(mean((leg_male_test$lncount_leg_male - leg_male_test$pred_gam_delta)^2, na.rm = T))
rmse_leg_male_delta # 9.28

# Variable plots
tiff(here('results/GAM',
          'gaussian_leg_male_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(leg_male_gam_abun, c(-3.5, 4))
dev.off()

tiff(here('results/GAM',
          'tweedie_leg_male_plots.jpg'),
     units = "in",
     width = 45,
     height = 40,
     res = 200)
variable_figure(leg_male_tweedie, c(-3, 2.5))
dev.off()

## Sublegal Male ----
# Gaussian
# Base model with presence/absence
sub_male_gam_base <- gam(pres_sub_male ~ s(year, bs = "re") +
                           s(longitude, latitude) +
                           s(julian),
                         data = sub_male_train,
                         family = "binomial")
summary(sub_male_gam_base) # 59.9% explained

# Abundance model
sub_male_gam_abun <- gam(lncount_sub_male ~ s(year, bs = "re") +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth, k = 5) +
                           s(phi, k = 5) +
                           s(temperature, k = 5) +
                           s(ice_mean, k = 5) +
                           s(sublegal_male_loading_station, k = 5) +
                           s(log_pcod_cpue, k = 5) +
                           s(bcs_sublegal_male, k = 5),
                         data = sub_male_train[sub_male_train$lncount_sub_male > 0, ])
summary(sub_male_gam_abun) # 61%

par(mfrow = c(2, 2))
gam.check(sub_male_gam_abun)

par(mfrow = c(3, 4))
plot(sub_male_gam_abun)

# Tweedie
sub_male_tweedie <- gam(sublegal_male + 1 ~ s(year, bs = 're') +
                          s(longitude, latitude) +
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
summary(sub_male_tweedie) # 70.6%

par(mfrow = c(2, 2))
gam.check(sub_male_tweedie)

par(mfrow = c(3, 4))
plot(sub_male_tweedie)

# Predict on test data
sub_male_test$pred_gam <- predict(sub_male_tweedie,
                                  sub_male_test,
                                  type = "link",
                                  exclude = "s(year)") 
sub_male_test$pred_gam_base <- predict(sub_male_gam_base,
                                       sub_male_test,
                                       type = "response",
                                       exclude = "s(year)")
sub_male_test$pred_gam_abun <- predict(sub_male_gam_abun,
                                       sub_male_test,
                                       type = "response",
                                       exclude = "s(year)")

sub_male_test$pred_gam_delta <- sub_male_test$pred_gam_base * sub_male_test$pred_gam_abun

rmse_sub_male_tweedie <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_gam)^2, na.rm = T))
rmse_sub_male_tweedie # 2.23

rmse_sub_male_delta <- sqrt(mean((sub_male_test$lncount_sub_male - sub_male_test$pred_gam_delta)^2, na.rm = T))
rmse_sub_male_delta # 10.12


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
rmse_mat_female_brt # 1.5

# Calculate deviance explained
dev_mat_female_abun <- brt_deviance(brt_mat_female_abun)
dev_mat_female_pres <- brt_deviance(brt_mat_female_base)

dev_mat_female_abun # 53.4% deviance explained
dev_mat_female_pres # 58.6% deviance explained

# Save models for future use
saveRDS(brt_mat_female_abun, file = here('data', 'brt_mat_female_abun.rds'))
saveRDS(brt_mat_female_base, file = here('data', 'brt_mat_female_base.rds'))

# Read in BRTs
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base.rds'))

# Variable names
column_names <- c("depth", "temperature", "phi", "ice_mean", "bcs_mature_female", "longitude",
                  "latitude", "julian", "female_loading_station", "log_pcod_cpue", "year_f")
final_names <- c("depth", "temperature", "phi", "ice concentration",
                 "proportion BCS", "longitude", "latitude", "julian",
                 "female loading", "log(cod cpue + 1)", "year")

column_labels <- data.frame(column_names, final_names)
match_labels <- column_labels[match(colnames(mat_female_train)[vars], column_labels$column_names), ]

brt_mat_female_summary <- summary(brt_mat_female_abun$model)

brt_labels <- match_labels[order(match(names(mat_female_train)[vars], brt_mat_female_summary$var)), ]
labels <- brt_labels$final_names

# Variable importance
rel_inf(brt_mat_female_abun, 'Variable Influence on Mature Female Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'female_mat_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_mat_female_abun)
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
rmse_imm_female_brt # 1.42

# Calculate deviance
dev_imm_female_abun <- brt_deviance(brt_imm_female_abun)
dev_imm_female_pres <- brt_deviance(brt_imm_female_base)

dev_imm_female_abun # 63% deviance explained
dev_imm_female_pres # 54% deviance explained

# Save models for future use
saveRDS(brt_imm_female_abun, file = here('data', 'brt_imm_female_abun.rds'))
saveRDS(brt_imm_female_base, file = here('data', 'brt_imm_female_base.rds'))

# Read in BRTs
brt_imm_female_abun <- readRDS(file = here('data', 'brt_imm_female_abun.rds'))
brt_imm_female_base <- readRDS(file = here('data', 'brt_imm_female_base.rds'))

# Variable importance
rel_inf(brt_imm_female_abun, 'Variable Influence on Immature Female Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'female_imm_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_imm_female_abun)
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

# Coefficient of determination
brt_leg_male_abun$model

# Save models for future use
saveRDS(brt_leg_male_abun, file = here('data', 'brt_leg_male_abun.rds'))
saveRDS(brt_leg_male_base, file = here('data', 'brt_leg_male_base.rds'))

# Read in BRTs
brt_leg_male_abun <- readRDS(file = here('data', 'brt_leg_male_abun.rds'))
brt_leg_male_base <- readRDS(file = here('data', 'brt_leg_male_base.rds'))

# Variable importance
rel_inf(brt_leg_male_abun, 'Variable Influence on Legal Male Snow Crab')
dev.copy(jpeg,
         here('results/BRT',
              'male_leg_rel_inf.jpg'),
         height = 12,
         width = 9,
         res = 200,
         units = 'in')
dev.off()

# Plot the variables
part_depen(brt_leg_male_abun)
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
spatial_grid_leg_male$legal_male_loading <- median(leg_male_train$legal_male_loading, na.rm = TRUE)
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
rmse_sub_male_brt # 1.60

# Calculate deviance explained
dev_sub_male_abun <- brt_deviance(brt_sub_male_abun)
dev_sub_male_pres <- brt_deviance(brt_sub_male_base)

dev_sub_male_abun # 74% deviance explained
dev_sub_male_pres # 59% deviance explained

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
leg_male_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                    "latitude", "julian", "legal_male_loading_station",
                    "bcs_legal_male", "log_pcod_cpue")
sub_male_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                    "latitude", "julian", "sublegal_male_loading_station",
                    "bcs_sublegal_male", "log_pcod_cpue")
mat_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", "female_loading_station",
                      "bcs_mature_female", "log_pcod_cpue", "year")
imm_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", "female_loading_station",
                      "bcs_immature_female", "log_pcod_cpue")
## Legal Males ----
### TreeSHAP ----
library(treeshap)
library(waterfalls)
leg_male_data <- crab_trans %>% # can use full data set
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male, ice_mean)
leg_male_explain <- leg_male_data[vars] %>% tidyr::drop_na() # cannot use year factor unless hot coded
leg_male_gbm <- gbm.unify(brt_leg_male_abun$model, leg_male_explain)
leg_male_gbm_abun <- treeshap(leg_male_gbm, leg_male_explain)
plot_contribution(leg_male_gbm_temp)
leg_male_sv <- shapviz(leg_male_gbm_temp)

sv_importance(leg_male_sv)
sv_importance(leg_male_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
sv_waterfall(leg_male_sv, 1) # one observation
sv_force(leg_male_sv, 2) # one observation
sv_dependence(leg_male_sv, 
              v = "temperature", 
              color_var = "ice_mean") # specific variable relationships
sv_dependence(leg_male_sv, v = leg_male_names)

sv_dependence(leg_male_sv,
              v = "temperature",
              color_var = NULL)

# Try waterfall split out by year
leg_male_shap <- leg_male_gbm_temp$shaps
leg_male_df <- cbind(leg_male_shap, year = leg_male_data[, 14])
leg_male_final <- leg_male_df[, c(11, 1:10)]
leg_male_waterfall <- leg_male_final %>% 
  group_by(year) %>%
  summarise(across(everything(), list(sum)))

waterfall(leg_male_waterfall)
ggplot(leg_male_df, aes(x = temperature, y = temp_shap)) +
  geom_point(aes(color = ice_mean))

### mSHAP ----
library(mshap)
leg_male_data <- crab_trans %>% # can use full data set
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male, ice_mean)
leg_male_explain <- leg_male_data[vars] %>% tidyr::drop_na()
leg_male_pres <- leg_male_data[, 13]
leg_male_abun <- leg_male_data[, 11]

leg_male_gbm_abun <- gbm.unify(brt_leg_male_abun$model, leg_male_explain)
leg_male_gbm_pres <- gbm.unify(brt_leg_male_base$model, leg_male_explain)

leg_male_mshap_abun <- treeshap(leg_male_gbm_abun, leg_male_explain)
leg_male_mshap_pres <- treeshap(leg_male_gbm_pres, leg_male_explain)

leg_male_mshap_comb <- mshap(shap_1 = leg_male_mshap_pres$shaps,
                             shap_2 = leg_male_mshap_abun$shaps,
                             ex_1 = leg_male_mshap_pres, # need to figure out how to get expected values
                             ex_2 = leg_male_train[, 11])

mshap::summary_plot(variable_values = leg_male_explain,
                    shap_values = leg_male_mshap_comb)

### KernelSHAP ----
leg_male_explain <- leg_male_train[vars] # only use columns in model
leg_male_x <- leg_male_explain[sample(nrow(leg_male_explain), 500), ]
leg_male_shap_abun <- kernelshap(brt_leg_male_abun$model, 
                                 leg_male_explain, 
                                 bg_X = leg_male_x)
saveRDS(leg_male_shap_abun, file = here('data', 'leg_male_shap_full.rds'))

# Visualize
leg_male_sv <- shapviz(leg_male_shap_abun)
# Gives the global effect of variables (absolute value, not directional)
sv_importance(leg_male_sv)
sv_importance(leg_male_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
sv_waterfall(leg_male_sv, 1) # one observation
sv_waterfall(leg_male_sv, leg_male_sv$X$year == "1996") # observations in one year
# Force plots 
# Yellow means variable pushes prediction higher, purple means variable pushes prediction lower
# Scores close to the f(x) value have more of an impact (indicated by SHAP magnitude too)
sv_force(leg_male_sv, 2) # one observation
sv_force(leg_male_sv, leg_male_sv$X$year == "1995") # observations in one year
sv_dependence(leg_male_sv, 
              v = "temperature", 
              color_var = "ice_mean") # specific variable relationships
sv_dependence(leg_male_sv, v = leg_male_names)

# Could I make a heatmap for years by variable? Show change in SHAP over time
# Could just do heatmap by stage/sex

## Sublegal Males ----
sub_male_explain <- sub_male_train[vars] # only use columns in model
sub_male_x <- sub_male_explain[sample(nrow(sub_male_explain), 500), ]
sub_male_shap <- kernelshap(brt_sub_male_abun$model, 
                            sub_male_explain, 
                            bg_X = sub_male_x)
saveRDS(sub_male_shap, file = here('data', 'sub_male_shap.rds'))

## Mature Females ----
mat_female_explain <- mat_female_train[vars] # only use columns in model
mat_female_x <- mat_female_explain[sample(nrow(mat_female_explain), 500), ]
system.time(mat_female_shap <- kernelshap(brt_mat_female_abun$model,
                                          mat_female_explain,
                                          bg_X = mat_female_x))
saveRDS(mat_female_shap, file = here('data', 'mat_female_shap_full.rds'))

mat_female_shap <- readRDS(file = here('data', 'mat_female_shap_full.rds'))

# Visualize
mat_female_sv <- shapviz(mat_female_shap)
mat_female_subgroups <- c("2017" = mat_female_sv[mat_female_test$year_f == "2017"],
                          "2018" = mat_female_sv[mat_female_test$year_f == "2018"],
                          "2019" = mat_female_sv[mat_female_test$year_f == "2019"])
mat_female_levels <- setNames(mat_female_shap, 
                              levels(mat_female_explain$year_f))
# Gives the global effect of variables (absolute value, not directional)
sv_importance(mat_female_sv)
sv_importance(mat_female_sv, kind = "bee") 
sv_waterfall(mat_female_sv, 1) # one observation
sv_waterfall(mat_female_sv, mat_female_sv$X$year == "1999") # observations in one year
sv_force(mat_female_sv, mat_female_sv$X$year == "1999")

sv_force(mat_female_subgroups, v = "temperature") +
  plot_layout(nrow = 3,
              ncol = 1,
              widths = c(1))

sv_dependence(mat_female_sv, 
              v = "temperature", 
              color_var = "auto")
sv_dependence(mat_female_sv, v = mat_female_names)

## Immature Females ----
imm_female_explain <- imm_female_train[vars] # only use columns in model
imm_female_x <- imm_female_explain[sample(nrow(imm_female_explain), 500), ]
system.time(imm_female_shap <- kernelshap(brt_imm_female_abun$model,
                                          imm_female_explain,
                                          bg_X = imm_female_x))
saveRDS(imm_female_shap, file = here('data', 'imm_female_shap.rds'))