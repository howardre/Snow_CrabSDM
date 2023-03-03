# Title: Snow Crab Model Exploration
# Purpose: Investigate potential models
# Data created: 07/20/2022

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(randomForest)
library(here)
library(mgcv)
library(dplyr)
library(colorspace)
library(maps)
library(mapdata)
library(fields)
library(ggplot2)
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

# Load data ----
# Make sure to run PCA first if updating the data matching script
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds'))

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
         !is.na(depth)) # The missing depth and temperature values end up removing the ability to use BCS

# Look at correlation
pairs(crab_trans[, c(7, 8, 22:25, 39)], cex = 0.3) # potential issue with phi and depth

summary(lm(phi ~ depth, data = crab_trans)) # R2: 0.33
summary(lm(log_pcod_cpue ~ temperature, data = crab_trans)) # R2: 0.0006

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
# Need to filter out stations without observer data in order to get accurate comparison
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
         pres_mat_female, year_f) %>%
  tidyr::drop_na(lncount_mat_female) 

imm_female_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_train <- crab_train %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)

# Test data
mat_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f) %>%
  tidyr::drop_na(lncount_mat_female) 

imm_female_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian, female_loading,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f) %>%
  tidyr::drop_na(lncount_imm_female)

leg_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f) %>%
  tidyr::drop_na(lncount_leg_male)

sub_male_test <- crab_test %>%
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, sublegal_male_loading,
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f) %>%
  tidyr::drop_na(lncount_sub_male) %>%
  dplyr::rename(sublegal_male = sublegal_male)


# GAMs ----
## Mature Female ----
# Gaussian
hist(mat_female_train$lncount_mat_female,
     main = "Mature Female Catch",
     xlab = "log(count+1)") # zero-inflated
hist(mat_female_train$lncount_mat_female[mat_female_train$lncount_mat_female > 0],
     main = "Mature Female Catch",
     xlab = "log(count+1)")

# Base model with presence/absence
mat_female_gam_base <- gam(pres_mat_female ~ year_f +
                             s(longitude, latitude) +
                             s(julian),
                           data = mat_female_train,
                           family = "binomial")
summary(mat_female_gam_base) # 54.2% explained

# Add environmental data
mat_female_gam1 <- gam(lncount_mat_female ~ year_f + 
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean),
                       data = mat_female_train[mat_female_train$lncount_mat_female > 0, ])
summary(mat_female_gam1) # 37.1% explained

# Add obs data
mat_female_gam2 <- gam(lncount_mat_female ~ year_f +
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean) +
                         s(longitude, latitude, by = female_loading),
                       data = mat_female_train[mat_female_train$lncount_mat_female > 0, ])
summary(mat_female_gam2) # 38.8%

# Add pcod and bcs data
mat_female_gam3 <- gam(lncount_mat_female ~ year_f +
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean) +
                         s(longitude, latitude, by = female_loading) +
                         s(log_pcod_cpue) +
                         s(bcs_mature_female),
                       data = mat_female_train[mat_female_train$lncount_mat_female > 0, ])
summary(mat_female_gam3) # 40.8%

par(mfrow = c(2, 2))
gam.check(mat_female_gam3)

par(mfrow = c(3, 3))
plot(mat_female_gam3)

# Tweedie
# Base model
mat_female_tweedie <- gam(mature_female + 1 ~ year_f +
                            s(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = mat_female_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(mat_female_tweedie) # 63% explained

# Add environmental data
mat_female_tweedie1 <- gam(mature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean),
                           data = mat_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie1) # 65.5% explained

# Add obs data
mat_female_tweedie2 <- gam(mature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean) +
                             s(longitude, latitude, by = female_loading),
                           data = mat_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie2) # 67%

# Add pcod data
mat_female_tweedie3 <- gam(mature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean) +
                             s(longitude, latitude, by = female_loading) +
                             s(log_pcod_cpue) +
                             s(bcs_mature_female),
                           data = mat_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie3) # 67.6%

par(mfrow = c(2, 2))
gam.check(mat_female_tweedie3)

par(mfrow = c(3, 3))
plot(mat_female_tweedie3)

# Make map
windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(mat_female_tweedie3,
          view = c('longitude', 'latitude'),
          too.far = 0.07,
          plot.type = 'contour',
          contour.col = contour_col,
          color = "jet" ,
          type = 'link',
          xlim = c(-181, -156),
          ylim = range(mat_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
symbols(mat_female_train$longitude[mat_female_train$lncount_mat_female > 0],
        mat_female_train$latitude[mat_female_train$lncount_mat_female > 0],
        circles = log(mat_female_train$lncount_mat_female + 1)[mat_female_train$lncount_mat_female > 0],
        inches = 0.1,
        bg = alpha('grey', 0.3),
        fg = alpha('black', 0.1),
        add = T)
points(mat_female_train$longitude[mat_female_train$lncount_mat_female == 0],
       mat_female_train$latitude[mat_female_train$lncount_mat_female == 0],
       pch =  '')
maps::map('worldHires',
          add = T,
          col = 'antiquewhite4',
          fill = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(mat_female_tweedie3$linear.predictors),
                    max(mat_female_tweedie3$linear.predictors)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

# Prediction grid map
nlat = 40
nlon = 60
latd = seq(min(mat_female_train$latitude), max(mat_female_train$latitude), length.out = nlat)
lond = seq(min(mat_female_train$longitude), max(mat_female_train$longitude), length.out = nlon)
spatial_grid_mat_female <- expand.grid(lond, latd)
names(spatial_grid_mat_female) <- c('longitude', 'latitude')
spatial_grid_mat_female$dist <- NA
for (k in 1:nrow(spatial_grid_mat_female)) {
  dist <-  distance_function(spatial_grid_mat_female$latitude[k],
                             spatial_grid_mat_female$longitude[k],
                             mat_female_train$latitude,
                             mat_female_train$longitude)
  spatial_grid_mat_female$dist[k] <- min(dist)
}
spatial_grid_mat_female$year_f <- as.factor('2010')
spatial_grid_mat_female$bcs_prop_f <- as.factor("Y")
spatial_grid_mat_female$depth <- median(mat_female_train$depth, na.rm = T)
spatial_grid_mat_female$phi <- median(mat_female_train$phi, na.rm = T)
spatial_grid_mat_female$julian <- median(mat_female_train$julian, na.rm = T)
spatial_grid_mat_female$temperature <- median(mat_female_train$temperature, na.rm = T)
spatial_grid_mat_female$ice_mean <- median(mat_female_train$ice_mean, na.rm = T)
spatial_grid_mat_female$female_loading <- median(mat_female_train$female_loading, na.rm = T) 
spatial_grid_mat_female$log_pcod_cpue <- median(mat_female_train$log_pcod_cpue, na.rm = T)

# Tidy grid
spatial_grid_tidy <- mat_female_train %>%
  mutate(longitude = cut(longitude, breaks = seq(floor(min(mat_female_train$longitude)), 
                                                 ceiling(max(mat_female_train$longitude)), by = 0.02),
                         labels = seq(floor(min(mat_female_train$longitude)) + 0.01, 
                                      ceiling(max(mat_female_train$longitude)), by = 0.02)),
         latitude = cut(latitude, breaks = seq(floor(min(mat_female_train$latitude)), 
                                               ceiling(max(mat_female_train$latitude)), by = 0.02),
                        labels = seq(floor(min(mat_female_train$latitude)) + 0.01,
                                     ceiling(max(mat_female_train$latitude)), by = 0.02))) %>%
  group_by(longitude, latitude) %>%
  summarise(depth = median(depth, na.rm = T),
            phi = median(phi, na.rm = T),
            temperature = median(temperature, na.rm = T),
            ice_mean = median(ice_mean, na.rm = T),
            .groups = "drop") %>%
  mutate(across(where(is.factor), ~as.numeric(as.character(.x))))

spatial_grid_tidy$year_f <- as.factor('2010')
spatial_grid_tidy$bcs_prop_f <- as.factor('Y')
spatial_grid_tidy$julian <- median(mat_female_train$julian, na.rm = T)
spatial_grid_tidy$female_loading <- median(mat_female_train$female_loading, na.rm = T)
spatial_grid_tidy$log_pcod_cpue <- median(mat_female_train$log_pcod_cpue, na.rm = T)
spatial_grid_tidy$dist <- NA
for (k in 1:nrow(spatial_grid_tidy)) {
  dist <-  distance_function(spatial_grid_tidy$latitude[k],
                             spatial_grid_tidy$longitude[k],
                             mat_female_train$latitude,
                             mat_female_train$longitude)
  spatial_grid_tidy$dist[k] <- min(dist)
}

spatial_grid_tidy$pred <- predict(mat_female_tweedie3,
                                        spatial_grid_tidy,
                                        type = "response")

spatial_grid_tidy$pred[spatial_grid_tidy$dist > 30000] <- NA

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(spatial_grid_tidy$pred, na.rm = T),
                               max(spatial_grid_tidy$pred, na.rm = T))))
color_sequence = seq(max(spatial_grid_tidy$pred, na.rm = T), 
                     min(spatial_grid_tidy$pred, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(spatial_grid_tidy$pred, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(spatial_grid_tidy$longitude,
      spatial_grid_tidy$latitude,
      t(matrix(spatial_grid_tidy$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(mat_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(spatial_grid_tidy$longitude,
      spatial_grid_tidy$latitude,
      t(matrix(spatial_grid_tidy$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-181, -156),
      ylim = range(mat_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Mature Females (GAM)",
      cex.main = 1.7,
      cex.lab = 1.7,
      cex.axis = 1.5)
# symbols(mat_female_train$longitude[mat_female_train$lncount_mat_female > 0],
#         mat_female_train$latitude[mat_female_train$lncount_mat_female > 0],
#         circles = log(mat_female_train$lncount_mat_female + 1)[mat_female_train$lncount_mat_female > 0],
#         inches = 0.1,
#         bg = alpha('grey', 0.3),
#         fg = alpha('black', 0.1),
#         add = T)
# points(mat_female_train$longitude[mat_female_train$lncount_mat_female == 0],
#        mat_female_train$latitude[mat_female_train$lncount_mat_female == 0],
#        pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_tidy$pred, na.rm = T), 
                    max(spatial_grid_tidy$pred, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

## Immature Female ----
# Gaussian
hist(imm_female_train$lncount_imm_female,
     main = "Immature Female Catch",
     xlab = "log(count+1)") # zero-inflated
hist(imm_female_train$lncount_imm_female[imm_female_train$lncount_imm_female > 0])

imm_female_gam_base <- gam(pres_imm_female ~ year_f +
                             s(longitude, latitude) +
                             s(julian),
                           data = imm_female_train,
                           family = "binomial")
summary(imm_female_gam_base) # 46.4% explained

# Add environmental data
imm_female_gam1 <- gam(lncount_imm_female ~ year_f + 
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean),
                       data = imm_female_train[imm_female_train$lncount_imm_female > 0, ])
summary(imm_female_gam1) # 46.7% explained

# Add obs data
imm_female_gam2 <- gam(lncount_imm_female ~ year_f +
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean) +
                         s(longitude, latitude, by = female_loading),
                       data = imm_female_train[imm_female_train$lncount_imm_female > 0, ])
summary(imm_female_gam2) # 48.5%

# Add pcod and bcs data
imm_female_gam3 <- gam(lncount_imm_female ~ year_f +
                         s(longitude, latitude) +
                         s(julian, k = 4) +
                         s(depth, k = 4) +
                         s(phi, k = 4) +
                         s(temperature, k = 4) +
                         s(ice_mean) +
                         s(longitude, latitude, by = female_loading) +
                         s(log_pcod_cpue) +
                         s(bcs_immature_female, k = 4),
                       data = imm_female_train[imm_female_train$lncount_imm_female > 0, ])
summary(imm_female_gam3) # 50.2%

par(mfrow = c(2, 2))
gam.check(imm_female_gam3)

par(mfrow = c(3, 3))
plot(imm_female_gam3)

# Tweedie
# Base model
imm_female_tweedie <- gam(immature_female + 1 ~ year_f +
                            s(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = imm_female_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(imm_female_tweedie) # 70.1% explained

# Add environmental data
imm_female_tweedie1 <- gam(immature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean),
                           data = imm_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie1) # 70% explained

# Add obs data
imm_female_tweedie2 <- gam(immature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean) +
                             s(longitude, latitude, by = female_loading),
                           data = imm_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie2) # 71.3%

# Add pcod and bcs data
imm_female_tweedie3 <- gam(immature_female + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean) +
                             s(longitude, latitude, by = female_loading) +
                             s(log_pcod_cpue) +
                             s(bcs_immature_female),
                           data = imm_female_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie3) # 71.4%

par(mfrow = c(2, 2))
gam.check(imm_female_tweedie3)

par(mfrow = c(3, 3))
plot(imm_female_tweedie3)

# Make maps
windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(imm_female_tweedie3,
          view = c('longitude', 'latitude'),
          too.far = 0.07,
          plot.type = 'contour',
          contour.col = contour_col,
          color = "jet" ,
          type = 'link',
          xlim = c(-181, -156),
          ylim = range(imm_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
symbols(imm_female_train$longitude[imm_female_train$lncount_imm_female > 0],
        imm_female_train$latitude[imm_female_train$lncount_imm_female > 0],
        circles = log(imm_female_train$lncount_imm_female + 1)[imm_female_train$lncount_imm_female > 0],
        inches = 0.1,
        bg = alpha('grey', 0.3),
        fg = alpha('black', 0.1),
        add = T)
points(imm_female_train$longitude[imm_female_train$lncount_imm_female == 0],
       imm_female_train$latitude[imm_female_train$lncount_imm_female == 0],
       pch =  '')
maps::map('worldHires',
          add = T,
          col = 'antiquewhite4',
          fill = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(imm_female_tweedie3$linear.predictors),
                    max(imm_female_tweedie3$linear.predictors)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))


## Legal Male ----
# Gaussian
hist(leg_male_train$lncount_leg_male,
     main = "Legal Male Catch",
     xlab = "log(count+1)") 
hist(leg_male_train$lncount_leg_male[leg_male_train$lncount_leg_male > 0])

leg_male_gam_base <- gam(pres_leg_male ~ year_f +
                           s(longitude, latitude) +
                           s(julian),
                         data = leg_male_train,
                         family = "binomial")
summary(leg_male_gam_base) # 58.4% explained

# Add environmental data
leg_male_gam1 <- gam(lncount_leg_male ~ year_f + 
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean),
                       data = leg_male_train[leg_male_train$lncount_leg_male > 0, ])
summary(leg_male_gam1) # 49.2% explained

# Add obs data
leg_male_gam2 <- gam(lncount_leg_male ~ year_f +
                         s(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_mean) +
                         s(longitude, latitude, by = legal_male_loading),
                       data = leg_male_train[leg_male_train$lncount_leg_male > 0, ])
summary(leg_male_gam2) # 51.8% explained

# Add pcod and bcs data
leg_male_gam3 <- gam(lncount_leg_male ~ year_f +
                       s(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_mean) +
                       s(longitude, latitude, by = legal_male_loading) +
                       s(log_pcod_cpue) +
                       s(bcs_legal_male, k = 4),
                     data = leg_male_train[leg_male_train$lncount_leg_male > 0, ])
summary(leg_male_gam3) # 51.8% explained

par(mfrow = c(2, 2))
gam.check(leg_male_gam3)

par(mfrow = c(3, 3))
plot(leg_male_gam3) # less impact of pcod on legal males

# Tweedie
# Base model
leg_male_tweedie <- gam(legal_male + 1 ~ year_f +
                            s(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = leg_male_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(leg_male_tweedie) # 62.5% explained

# Add environmental data
leg_male_tweedie1 <- gam(legal_male + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean),
                           data = leg_male_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(leg_male_tweedie1) # 66% explained

# Add obs data
leg_male_tweedie2 <- gam(legal_male + 1 ~ year_f +
                             s(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_mean) +
                             s(longitude, latitude, by = legal_male_loading),
                           data = leg_male_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(leg_male_tweedie2) # 67.5%

# Add pcod data
leg_male_tweedie3 <- gam(legal_male + 1 ~ year_f +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_mean) +
                           s(longitude, latitude, by = legal_male_loading) +
                           s(log_pcod_cpue) +
                           s(bcs_legal_male),
                         data = leg_male_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(leg_male_tweedie3) # 67.3%

par(mfrow = c(2, 2))
gam.check(leg_male_tweedie3)

par(mfrow = c(3, 3))
plot(leg_male_tweedie3)

# Make maps
windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(leg_male_tweedie3,
          view = c('longitude', 'latitude'),
          too.far = 0.07,
          plot.type = 'contour',
          contour.col = contour_col,
          color = "jet" ,
          type = 'link',
          xlim = c(-181, -156),
          ylim = range(leg_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
symbols(leg_male_train$longitude[leg_male_train$lncount_leg_male > 0],
        leg_male_train$latitude[leg_male_train$lncount_leg_male > 0],
        circles = log(leg_male_train$lncount_leg_male + 1)[leg_male_train$lncount_leg_male > 0],
        inches = 0.1,
        bg = alpha('grey', 0.3),
        fg = alpha('black', 0.1),
        add = T)
points(leg_male_train$longitude[leg_male_train$lncount_leg_male == 0],
       leg_male_train$latitude[leg_male_train$lncount_leg_male == 0],
       pch =  '')
maps::map('worldHires',
          add = T,
          col = 'antiquewhite4',
          fill = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(leg_male_tweedie3$linear.predictors),
                    max(leg_male_tweedie3$linear.predictors)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

## Sublegal Male ----
# Gaussian
hist(sub_male_train$lncount_sub_male,
     main = "Sublegal Male Catch",
     xlab = "log(count+1)") 
hist(sub_male_train$lncount_sub_male[sub_male_train$lncount_sub_male > 0])

sub_male_gam_base <- gam(pres_sub_male ~ year_f +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth),
                         data = sub_male_train,
                         family = "binomial")
summary(sub_male_gam_base) # 64.7% explained

# Add environmental data
sub_male_gam1 <- gam(lncount_sub_male ~ year_f + 
                       s(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_mean),
                     data = sub_male_train[sub_male_train$lncount_sub_male > 0, ])
summary(sub_male_gam1) # 64.3% explained

# Add obs data
sub_male_gam2 <- gam(lncount_sub_male ~ year_f +
                       s(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_mean) +
                       s(longitude, latitude, by = sublegal_male_loading),
                     data = sub_male_train[sub_male_train$lncount_sub_male > 0, ])
summary(sub_male_gam2) # 65.6% explained

# Add pcod data
sub_male_gam3 <- gam(lncount_sub_male ~ year_f +
                       s(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_mean) +
                       s(longitude, latitude, by = sublegal_male_loading) +
                       s(log_pcod_cpue) +
                       s(bcs_sublegal_male),
                     data = sub_male_train[sub_male_train$lncount_sub_male > 0, ])
summary(sub_male_gam3) # 65.9% explained

par(mfrow = c(2, 2))
gam.check(sub_male_gam3)

par(mfrow = c(3, 3))
plot(sub_male_gam3)

# Tweedie
# Base model
sub_male_tweedie <- gam(sublegal_male + 1 ~ year_f +
                          s(longitude, latitude) +
                          s(julian) +
                          s(depth),
                        data = sub_male_train,
                        family = tw(link = "log"),
                        method = "REML")
summary(sub_male_tweedie) # 72% explained

# Add environmental data
sub_male_tweedie1 <- gam(sublegal_male + 1 ~ year_f +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_mean),
                         data = sub_male_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie1) # 73.4% explained

# Add obs data
sub_male_tweedie2 <- gam(sublegal_male + 1 ~ year_f +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_mean) +
                           s(longitude, latitude, by = sublegal_male_loading),
                         data = sub_male_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie2) # 74.5%

# Add pcod data
sub_male_tweedie3 <- gam(sublegal_male + 1 ~ year_f +
                           s(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_mean) +
                           s(longitude, latitude, by = sublegal_male_loading) +
                           s(log_pcod_cpue) +
                           s(bcs_sublegal_male),
                         data = sub_male_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie3) # 74.1%

par(mfrow = c(2, 2))
gam.check(sub_male_tweedie3)

par(mfrow = c(3, 3))
plot(sub_male_tweedie3)

# Make maps
windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(sub_male_tweedie3,
          view = c('longitude', 'latitude'),
          too.far = 0.07,
          plot.type = 'contour',
          contour.col = contour_col,
          color = "jet" ,
          type = 'link',
          xlim = c(-181, -156),
          ylim = range(sub_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
symbols(sub_male_train$longitude[sub_male_train$lncount_sub_male > 0],
        sub_male_train$latitude[sub_male_train$lncount_sub_male > 0],
        circles = log(sub_male_train$lncount_sub_male + 1)[sub_male_train$lncount_sub_male > 0],
        inches = 0.1,
        bg = alpha('grey', 0.3),
        fg = alpha('black', 0.1),
        add = T)
points(sub_male_train$longitude[sub_male_train$lncount_sub_male == 0],
       sub_male_train$latitude[sub_male_train$lncount_sub_male == 0],
       pch =  '')
maps::map('worldHires',
          add = T,
          col = 'antiquewhite4',
          fill = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(sub_male_tweedie3$linear.predictors),
                    max(sub_male_tweedie3$linear.predictors)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))



# Boosted regression trees ----
# Adjust the bag fraction to a value between 0.5-0.75 as suggested by Elith et al. (2008)
# The learning rate could range from 0.1-0.0001, higher value usually means less trees
# Depending on the number of samples, want tree complexity to be high enough (likely using 5)
# Want at least 1000 trees, but don't need to go way beyond it

library(enmSdmX) # use for grid search, wrapper for dismo

grid_search <- function(data, response, family){
  trainBRT(data = data,
           preds = c(1:10, 14),
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

## Mature females ----
# Get best models
brt_mat_females_base <- grid_search(mat_female_train, 13, 'bernoulli')
brt_mat_females_base

brt_mat_females_abun <- grid_search(mat_female_train[mat_female_train$lncount_mat_female > 0, ],
                                    11, 'gaussian')
brt_mat_females_abun

# Predict on test data
mat_female_test$pred_base <- predict.gbm(brt_mat_females_base$model,
                                         mat_female_test,
                                         n.trees = brt_mat_females_base$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred_abun <- predict.gbm(brt_mat_females_abun$model,
                                         mat_female_test,
                                         n.trees = brt_mat_females_abun$model$gbm.call$best.trees,
                                         type = "response")

mat_female_test$pred <- mat_female_test$pred_base * mat_female_test$pred_abun


# Calculate RMSE
rmse_mat_female <- sqrt(mean((mat_female_test$pres_mat_female - mat_female_test$pred_base)^2))
rmse_mat_female


# Attempt dropping variable
females_mat_simp <- gbm.simplify(brt_mat_females4, n.drops = 8) # this takes forever
summary(females_mat_simp)

# Choose final model
females_mat_final <- brt_mat_females4 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_mat_final,
         n.plots = 8,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
female_mat_effects <- tibble::as_tibble(summary.gbm(females_mat_final, plotit = F))
female_mat_effects %>% arrange(desc(rel.inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var,
                                      .x = rel.inf),
             y = rel.inf,
             fill = rel.inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.title = element_text()) +
  xlab('Variable') +
  ylab('Relative Influence') +
  ggtitle('Variable Influence on Mature Female Snow Crab')

# Plot the fits
females_mat_int <- gbm.interactions(females_mat_final)
females_mat_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(females_mat_final,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_mat_final,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_mat_final,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Training data map
nlat = 40
nlon = 60
latd = seq(min(mat_female_train$latitude), max(mat_female_train$latitude), length.out = nlat)
lond = seq(min(mat_female_train$longitude), max(mat_female_train$longitude), length.out = nlon)
spatial_grid_mat_female <- expand.grid(lond, latd)
names(spatial_grid_mat_female) <- c('longitude', 'latitude')
spatial_grid_mat_female$dist <- NA
for (k in 1:nrow(spatial_grid_mat_female)) {
  dist <-  distance_function(spatial_grid_mat_female$latitude[k],
                             spatial_grid_mat_female$longitude[k],
                             mat_female_train$latitude,
                             mat_female_train$longitude)
  spatial_grid_mat_female$dist[k] <- min(dist)
}
spatial_grid_mat_female$year_f <- as.factor('2010')
spatial_grid_mat_female$depth <- median(mat_female_train$depth, na.rm = T)
spatial_grid_mat_female$phi <- median(mat_female_train$phi, na.rm = T)
spatial_grid_mat_female$julian <- median(mat_female_train$julian, na.rm = T)
spatial_grid_mat_female$temperature <- median(mat_female_train$temperature, na.rm = T)
spatial_grid_mat_female$ice_mean <- median(mat_female_train$ice_mean, na.rm = T)
spatial_grid_mat_female$female_loading <- median(mat_female_train$female_loading, na.rm = T) 
spatial_grid_mat_female$log_pcod_cpue <- median(mat_female_train$log_pcod_cpue, na.rm = T)

# Need to add year as a constant because it is included as a factor (const = add)
# Must also be in the prediction grid as a factor 
# If neither of the above is done, RStudio will crash
year_f <- factor('2012', levels = levels(mat_female_train$year_f))
add <- data.frame(year_f)

preds_mat_female <- predict(females_mat_final,
                            spatial_grid_mat_female,
                            n.trees = females_mat_final$gbm.call$best.trees, 
                            type = "response",
                            const = add)
base_mat_female <- predict.gbm(brt_mat_females_base,
                               spatial_grid_mat_female,
                               n.trees = brt_mat_females_base$gbm.call$best.trees, 
                               type = "response",
                               const = add)
spatial_grid_mat_female$pred <- preds_mat_female * base_mat_female


spatial_grid_mat_female$pred[spatial_grid_mat_female$dist > 30000] <- NA

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(spatial_grid_mat_female$pred, na.rm = T),
                               max(spatial_grid_mat_female$pred, na.rm = T))))
color_sequence = seq(max(spatial_grid_mat_female$pred, na.rm = T), 
                     min(spatial_grid_mat_female$pred, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(spatial_grid_mat_female$pred, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(spatial_grid_mat_female$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(mat_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(lond,
      latd,
      t(matrix(spatial_grid_mat_female$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-181, -156),
      ylim = range(mat_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Mature Females (BRT)",
      cex.main = 1.7,
      cex.lab = 1.7,
      cex.axis = 1.5)
# symbols(mat_female_train$longitude[mat_female_train$lncount_mat_female > 0],
#         mat_female_train$latitude[mat_female_train$lncount_mat_female > 0],
#         circles = log(mat_female_train$lncount_mat_female + 1)[mat_female_train$lncount_mat_female > 0],
#         inches = 0.1,
#         bg = alpha('grey', 0.3),
#         fg = alpha('black', 0.1),
#         add = T)
# points(mat_female_train$longitude[mat_female_train$lncount_mat_female == 0],
#        mat_female_train$latitude[mat_female_train$lncount_mat_female == 0],
#        pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_mat_female$pred, na.rm = T), 
                    max(spatial_grid_mat_female$pred, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))


## Immature females ----
brt_imm_females_base <- gbm.step(data = imm_female_train,
                                gbm.x = c(1:10, 13),
                                gbm.y = 12,
                                family = 'bernoulli',
                                tree.complexity = 5,
                                learning.rate = 0.05,
                                bag.fraction = 0.5)

brt_imm_females1 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_imm_females1) 


brt_imm_females2 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females2)


brt_imm_females3 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 7,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_imm_females3)


brt_imm_females4 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10, 13),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 10,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_imm_females4)


brt_imm_females5 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.75)
summary(brt_imm_females5)


brt_imm_females6 <- gbm.step(data = imm_female_train[imm_female_train$lncount_imm_female > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.25)
summary(brt_imm_females6)

# Attempt dropping variable
females_imm_simp <- gbm.simplify(brt_imm_females4, n.drops = 8) # this takes forever
summary(females_imm_simp)

# Choose final model
females_imm_final <- brt_imm_females4 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_imm_final,
         n.plots = 8,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
female_imm_effects <- tibble::as_tibble(summary.gbm(females_imm_final, plotit = F))
female_imm_effects %>% arrange(desc(rel.inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var,
                                      .x = rel.inf),
             y = rel.inf,
             fill = rel.inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.title = element_text()) +
  xlab('Variable') +
  ylab('Relative Influence') +
  ggtitle('Variable Influence on Immature Female Snow Crab')

# Plot the fits
females_imm_int <- gbm.interactions(females_imm_final)
females_imm_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(females_imm_final,
            7, 2,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_imm_final,
            6, 1,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_imm_final,
            2, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Training data map
nlat = 40
nlon = 60
latd = seq(min(imm_female_train$latitude), max(imm_female_train$latitude), length.out = nlat)
lond = seq(min(imm_female_train$longitude), max(imm_female_train$longitude), length.out = nlon)
spatial_grid_imm_female <- expand.grid(lond, latd)
names(spatial_grid_imm_female) <- c('longitude', 'latitude')
spatial_grid_imm_female$dist <- NA
for (k in 1:nrow(spatial_grid_imm_female)) {
  dist <-  distance_function(spatial_grid_imm_female$latitude[k],
                             spatial_grid_imm_female$longitude[k],
                             imm_female_train$latitude,
                             imm_female_train$longitude)
  spatial_grid_imm_female$dist[k] <- min(dist)
}
spatial_grid_imm_female$year_f <- as.factor('2010')
spatial_grid_imm_female$depth <- median(imm_female_train$depth, na.rm = T)
spatial_grid_imm_female$phi <- median(imm_female_train$phi, na.rm = T)
spatial_grid_imm_female$julian <- median(imm_female_train$julian, na.rm = T)
spatial_grid_imm_female$temperature <- median(imm_female_train$temperature, na.rm = T)
spatial_grid_imm_female$ice_mean <- median(imm_female_train$ice_mean, na.rm = T)
spatial_grid_imm_female$female_loading <- median(imm_female_train$female_loading, na.rm = T) 
spatial_grid_imm_female$log_pcod_cpue <- median(imm_female_train$log_pcod_cpue, na.rm = T)

preds_imm_female <- predict(females_imm_final,
                            spatial_grid_imm_female,
                            n.trees = females_imm_final$gbm.call$best.trees, 
                            type = "response",
                            const = add)
base_imm_female <- predict(brt_imm_females_base,
                           spatial_grid_imm_female,
                           n.trees = brt_imm_females_base$gbm.call$best.trees,
                           type = "response",
                           const = add)
spatial_grid_imm_female$pred <- preds_imm_female * base_imm_female

summary(preds_imm_female)

spatial_grid_imm_female$pred[spatial_grid_imm_female$dist > 30000] <- NA

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(spatial_grid_imm_female$pred, na.rm = T),
                               max(spatial_grid_imm_female$pred, na.rm = T))))
color_sequence = seq(max(spatial_grid_imm_female$pred, na.rm = T), 
                     min(spatial_grid_imm_female$pred, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(spatial_grid_imm_female$pred, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(spatial_grid_imm_female$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(imm_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(lond,
      latd,
      t(matrix(spatial_grid_imm_female$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-181, -156),
      ylim = range(imm_female_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Immature Females (BRT)",
      cex.main = 1.7,
      cex.lab = 1.7,
      cex.axis = 1.5)
# symbols(imm_female_train$longitude[imm_female_train$lncount_imm_female > 0],
#         imm_female_train$latitude[imm_female_train$lncount_imm_female > 0],
#         circles = log(imm_female_train$lncount_imm_female + 1)[imm_female_train$lncount_imm_female > 0],
#         inches = 0.1,
#         bg = alpha('grey', 0.3),
#         fg = alpha('black', 0.1),
#         add = T)
# points(imm_female_train$longitude[imm_female_train$lncount_imm_female == 0],
#        imm_female_train$latitude[imm_female_train$lncount_imm_female == 0],
#        pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_imm_female$pred, na.rm = T), 
                    max(spatial_grid_imm_female$pred, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))


## Legal Males ----
brt_leg_males_base <- gbm.step(data = leg_male_train,
                           gbm.x = c(1:10, 13),
                           gbm.y = 12,
                           family = 'bernoulli',
                           tree.complexity = 10,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_leg_males_base)

brt_leg_males1 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_leg_males1) 


brt_leg_males2 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males2)


brt_leg_males3 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 7,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_leg_males3)


brt_leg_males4 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0,],
                           gbm.x = c(1:10, 13),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 10,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_leg_males4)


brt_leg_males5 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.75)
summary(brt_leg_males5)


brt_leg_males6 <- gbm.step(data = leg_male_train[leg_male_train$lncount_leg_male > 0, ],
                             gbm.x = c(1:10),
                             gbm.y = 10,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.25)
summary(brt_leg_males6)

# Attempt dropping variable
leg_male_simp <- gbm.simplify(brt_leg_males4, n.drops = 8) # this takes forever
summary(leg_male_simp)

# Choose final model
leg_male_final <- brt_leg_males4 # Change this once decision made

# Plot the variables
windows()
gbm.plot(leg_male_final,
         n.plots = 8,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
leg_male_effects <- tibble::as_tibble(summary.gbm(leg_male_final, plotit = F))
leg_male_effects %>% arrange(desc(rel.inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var,
                                      .x = rel.inf),
             y = rel.inf,
             fill = rel.inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.title = element_text()) +
  xlab('Variable') +
  ylab('Relative Influence') +
  ggtitle('Variable Influence on Legal Male Snow Crab')

# Plot the fits
leg_male_int <- gbm.interactions(leg_male_final)
leg_male_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(leg_male_final,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(leg_male_final,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(leg_male_final,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Training data map
nlat = 40
nlon = 60
latd = seq(min(leg_male_train$latitude), max(leg_male_train$latitude), length.out = nlat)
lond = seq(min(leg_male_train$longitude), max(leg_male_train$longitude), length.out = nlon)
spatial_grid_leg_male <- expand.grid(lond, latd)
names(spatial_grid_leg_male) <- c('longitude', 'latitude')
spatial_grid_leg_male$dist <- NA
for (k in 1:nrow(spatial_grid_leg_male)) {
  dist <-  distance_function(spatial_grid_leg_male$latitude[k],
                             spatial_grid_leg_male$longitude[k],
                             leg_male_train$latitude,
                             leg_male_train$longitude)
  spatial_grid_leg_male$dist[k] <- min(dist)
}
spatial_grid_leg_male$year_f <- as.factor('2010')
spatial_grid_leg_male$depth <- median(leg_male_train$depth, na.rm = T)
spatial_grid_leg_male$phi <- median(leg_male_train$phi, na.rm = T)
spatial_grid_leg_male$julian <- median(leg_male_train$julian, na.rm = T)
spatial_grid_leg_male$temperature <- median(leg_male_train$temperature, na.rm = T)
spatial_grid_leg_male$ice_mean <- median(leg_male_train$ice_mean, na.rm = T)
spatial_grid_leg_male$legal_male_loading <- median(leg_male_train$legal_male_loading, na.rm = T) 
spatial_grid_leg_male$log_pcod_cpue <- median(leg_male_train$log_pcod_cpue, na.rm = T)

preds_leg_male <- predict.gbm(leg_male_final,
                              spatial_grid_leg_male,
                              n.trees = leg_male_final$gbm.call$best.trees,
                              type = "response",
                              const = add)
base_leg_male <- predict.gbm(brt_leg_males_base,
                             spatial_grid_leg_male,
                             n.trees = brt_leg_males_base$gbm.call$best.trees,
                             type = "response",
                             const = add)
spatial_grid_leg_male$pred <- preds_leg_male * base_leg_male

summary(preds_leg_male)

spatial_grid_leg_male$pred[spatial_grid_leg_male$dist > 30000] <- NA

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(spatial_grid_leg_male$pred, na.rm = T),
                               max(spatial_grid_leg_male$pred, na.rm = T))))
color_sequence = seq(max(spatial_grid_leg_male$pred, na.rm = T), 
                     min(spatial_grid_leg_male$pred, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(spatial_grid_leg_male$pred, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(spatial_grid_leg_male$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(leg_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(lond,
      latd,
      t(matrix(spatial_grid_leg_male$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-181, -156),
      ylim = range(leg_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Legal Males (BRT)",
      cex.main = 1.7,
      cex.lab = 1.7,
      cex.axis = 1.5)
# symbols(leg_male_train$longitude[leg_male_train$lncount_leg_male > 0],
#         leg_male_train$latitude[leg_male_train$lncount_leg_male > 0],
#         circles = log(leg_male_train$lncount_leg_male + 1)[leg_male_train$lncount_leg_male > 0],
#         inches = 0.1,
#         bg = alpha('grey', 0.3),
#         fg = alpha('black', 0.1),
#         add = T)
# points(leg_male_train$longitude[leg_male_train$lncount_leg_male == 0],
#        leg_male_train$latitude[leg_male_train$lncount_leg_male == 0],
#        pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_leg_male$pred, na.rm = T), 
                    max(spatial_grid_leg_male$pred, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))


## Sublegal Males ----
brt_sub_males_base <- gbm.step(data = sub_male_train,
                               gbm.x = c(1:10, 13),
                               gbm.y = 12,
                               family = 'bernoulli',
                               tree.complexity = 10,
                               learning.rate = 0.05,
                               bag.fraction = 0.5) 
summary(brt_sub_males_base)


brt_sub_males1 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_sub_males1) 


brt_sub_males2 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males2)


brt_sub_males3 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 7,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_sub_males3)


brt_sub_males4 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10, 13),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 10,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_sub_males4)


brt_sub_males5 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.05,
                           bag.fraction = 0.75)
summary(brt_sub_males5)


brt_sub_males6 <- gbm.step(data = sub_male_train[sub_male_train$lncount_sub_male > 0, ],
                           gbm.x = c(1:10),
                           gbm.y = 10,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.05,
                           bag.fraction = 0.25)
summary(brt_sub_males6)

# Attempt dropping variable
sub_male_simp <- gbm.simplify(brt_sub_males4, n.drops = 8) # this takes forever
summary(sub_male_simp)

# Choose final model
sub_male_final <- brt_sub_males4 # Change this once decision made

# Plot the variables
windows()
gbm.plot(sub_male_final,
         n.plots = 8,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
sub_male_effects <- tibble::as_tibble(summary.gbm(sub_male_final, plotit = F))
sub_male_effects %>% arrange(desc(rel.inf)) %>%
  ggplot(aes(x = forcats::fct_reorder(.f = var,
                                      .x = rel.inf),
             y = rel.inf,
             fill = rel.inf)) +
  geom_col() +
  coord_flip() +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.title = element_text()) +
  xlab('Variable') +
  ylab('Relative Influence') +
  ggtitle('Variable Influence on Sublegal Male Snow Crab')

# Plot the fits
sub_male_int <- gbm.interactions(sub_male_final)
sub_male_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(sub_male_final,
            2, 3,
            z.range = c(0, 6.25),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(sub_male_final,
            2, 7,
            z.range = c(-1.3, 4.8),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(sub_male_final,
            1, 3,
            z.range = c(-0.1, 5.7),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

# Training data map
nlat = 40
nlon = 60
latd = seq(min(sub_male_train$latitude), max(sub_male_train$latitude), length.out = nlat)
lond = seq(min(sub_male_train$longitude), max(sub_male_train$longitude), length.out = nlon)
spatial_grid_sub_male <- expand.grid(lond, latd)
names(spatial_grid_sub_male) <- c('longitude', 'latitude')
spatial_grid_sub_male$dist <- NA
for (k in 1:nrow(spatial_grid_sub_male)) {
  dist <-  distance_function(spatial_grid_sub_male$latitude[k],
                             spatial_grid_sub_male$longitude[k],
                             sub_male_train$latitude,
                             sub_male_train$longitude)
  spatial_grid_sub_male$dist[k] <- min(dist)
}
spatial_grid_sub_male$year_f <- as.factor('2010')
spatial_grid_sub_male$depth <- median(sub_male_train$depth, na.rm = T)
spatial_grid_sub_male$phi <- median(sub_male_train$phi, na.rm = T)
spatial_grid_sub_male$julian <- median(sub_male_train$julian, na.rm = T)
spatial_grid_sub_male$temperature <- median(sub_male_train$temperature, na.rm = T)
spatial_grid_sub_male$ice_mean <- median(sub_male_train$ice_mean, na.rm = T)
spatial_grid_sub_male$sublegal_male_loading <- median(sub_male_train$sublegal_male_loading, na.rm = T) 
spatial_grid_sub_male$log_pcod_cpue <- median(sub_male_train$log_pcod_cpue, na.rm = T)

preds_sub_male <- predict.gbm(sub_male_final,
                              spatial_grid_sub_male,
                              n.trees = sub_male_final$gbm.call$best.trees, 
                              type = "response",
                              const = add)
base_sub_male <- predict.gbm(brt_sub_males_base,
                             spatial_grid_sub_male,
                             n.trees = brt_sub_males_base$gbm.call$best.trees, 
                             type = "response",
                             const = add)
spatial_grid_sub_male$pred <- preds_sub_male * base_sub_male

summary(preds_sub_male)

spatial_grid_sub_male$pred[spatial_grid_sub_male$dist > 30000] <- NA

my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
color_levels = 100
max_absolute_value = max(abs(c(min(spatial_grid_sub_male$pred, na.rm = T),
                               max(spatial_grid_sub_male$pred, na.rm = T))))
color_sequence = seq(max(spatial_grid_sub_male$pred, na.rm = T), 
                     min(spatial_grid_sub_male$pred, na.rm = T),
                     length.out = color_levels + 1)
n_in_class = hist(spatial_grid_sub_male$pred, breaks = color_sequence, plot = F)$counts > 0
col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(spatial_grid_sub_male$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(sub_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(lond,
      latd,
      t(matrix(spatial_grid_sub_male$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(n = color_levels)[col_to_include],
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-181, -156),
      ylim = range(sub_male_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Sublegal Males (BRT)",
      cex.main = 1.7,
      cex.lab = 1.7,
      cex.axis = 1.5)
# symbols(sub_male_train$longitude[sub_male_train$lncount_sub_male > 0],
#         sub_male_train$latitude[sub_male_train$lncount_sub_male > 0],
#         circles = log(sub_male_train$lncount_sub_male + 1)[sub_male_train$lncount_sub_male > 0],
#         inches = 0.1,
#         bg = alpha('grey', 0.3),
#         fg = alpha('black', 0.1),
#         add = T)
# points(sub_male_train$longitude[sub_male_train$lncount_sub_male == 0],
#        sub_male_train$latitude[sub_male_train$lncount_sub_male == 0],
#        pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.28, .31, .27, .42),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.3,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_sub_male$pred, na.rm = T), 
                    max(spatial_grid_sub_male$pred, na.rm = T)),
           legend.args = list("log(count+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))