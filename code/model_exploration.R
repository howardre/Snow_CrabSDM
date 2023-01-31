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
                     lncount_sub_male = log(immature_male + 1),
                     pres_imm_female = ifelse(mature_female > 0, 1, 0),
                     pres_mat_female = ifelse(immature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(immature_male > 0, 1, 0))

# Look at correlation
pairs(crab_trans[, c(10, 11, 20, 21)], cex = 0.3) # potential issue with phi and depth

summary(lm(phi ~ depth, data = crab_trans)) # R2: 0.27

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
# Need to filter out stations without observer data in order to get accurate comparison
crab_train <- as.data.frame(crab_trans %>% 
                              filter(year < 2015))
crab_test <- as.data.frame(crab_trans %>% 
                             filter(year > 2014))

# GAMs ----
## Mature Female ----
# Gaussian
hist(crab_train$lncount_mat_female,
     main = "Mature Female Catch",
     xlab = "log(count+1)") # zero-inflated
hist(crab_train$lncount_mat_female[crab_train$lncount_mat_female > 0])

# Base model with presence/absence
mat_female_gam_base <- gam(pres_mat_female ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian),
                           data = crab_train,
                           family = "binomial")
summary(mat_female_gam_base) # 45.2% explained

# Add environmental data
mat_female_gam1 <- gam(lncount_mat_female ~ factor(year) + 
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index),
                       data = crab_train[crab_train$lncount_mat_female > 0, ])
summary(mat_female_gam1) # 37.1% explained

# Add obs data
mat_female_gam2 <- gam(lncount_mat_female ~ factor(year) +
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index) +
                         s(longitude, latitude, by = female_loading),
                       data = crab_train[crab_train$lncount_mat_female > 0, ])
summary(mat_female_gam2) # 38.8%

# Add pcod and bcs data
mat_female_gam3 <- gam(lncount_mat_female ~ factor(year) +
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index) +
                         s(longitude, latitude, by = female_loading) +
                         s(pcod_cpue),
                       data = crab_train[crab_train$lncount_mat_female > 0, ])
summary(mat_female_gam3) # 39.2%

par(mfrow = c(2, 2))
gam.check(mat_female_gam3)

par(mfrow = c(3, 3))
plot(mat_female_gam3)

# Tweedie
# Base model
mat_female_tweedie <- gam(mature_female + 1 ~ factor(year) +
                            te(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = crab_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(mat_female_tweedie) # 63.3% explained

# Add environmental data
mat_female_tweedie1 <- gam(mature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie1) # 64.5% explained

# Add obs data
mat_female_tweedie2 <- gam(mature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index) +
                             s(longitude, latitude, by = female_loading),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie2) # 67%

# Add pcod data
mat_female_tweedie3 <- gam(mature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index) +
                             s(longitude, latitude, by = female_loading) +
                             s(pcod_cpue),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie3) # 67.1%

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
          ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
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
## Immature Female ----
# Gaussian
hist(crab_train$lncount_imm_female,
     main = "Immature Female Catch",
     xlab = "log(count+1)") # zero-inflated
hist(crab_train$lncount_imm_female[crab_train$lncount_imm_female > 0])

imm_female_gam_base <- gam(pres_imm_female ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian),
                           data = crab_train,
                           family = "binomial")
summary(imm_female_gam_base) # 53.6% explained

# Add environmental data
imm_female_gam1 <- gam(lncount_imm_female ~ factor(year) + 
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index),
                       data = crab_train[crab_train$lncount_imm_female > 0, ])
summary(imm_female_gam1) # 44.8% explained

# Add obs data
imm_female_gam2 <- gam(lncount_imm_female ~ factor(year) +
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index) +
                         s(longitude, latitude, by = female_loading),
                       data = crab_train[crab_train$lncount_imm_female > 0, ])
summary(imm_female_gam2) # 47.1%

# Add pcod data
imm_female_gam3 <- gam(lncount_imm_female ~ factor(year) +
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index) +
                         s(longitude, latitude, by = female_loading) +
                         s(pcod_cpue),
                       data = crab_train[crab_train$lncount_imm_female > 0, ])
summary(imm_female_gam3) # 47.9%

par(mfrow = c(2, 2))
gam.check(imm_female_gam3)

par(mfrow = c(3, 3))
plot(imm_female_gam3)

# Tweedie
# Base model
imm_female_tweedie <- gam(immature_female + 1 ~ factor(year) +
                            te(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = crab_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(imm_female_tweedie) # 67.6% explained

# Add environmental data
imm_female_tweedie1 <- gam(immature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie1) # 68.7% explained

# Add obs data
imm_female_tweedie2 <- gam(immature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index) +
                             te(longitude, latitude, by = female_loading),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie2) # 69.9%

# Add pcod data
imm_female_tweedie3 <- gam(immature_female + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index) +
                             te(longitude, latitude, by = female_loading) +
                             s(pcod_cpue),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(imm_female_tweedie3) # 70.1%

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
          ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
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
hist(crab_train$lncount_leg_male,
     main = "Legal Male Catch",
     xlab = "log(count+1)") 
hist(crab_train$lncount_leg_male[crab_train$lncount_leg_male > 0])

leg_male_gam_base <- gam(pres_leg_male ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian),
                         data = crab_train,
                         family = "binomial")
summary(leg_male_gam_base) # 55.6% explained

# Add environmental data
leg_male_gam1 <- gam(lncount_leg_male ~ factor(year) + 
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index),
                       data = crab_train[crab_train$lncount_leg_male > 0, ])
summary(leg_male_gam1) # 48.1% explained

# Add obs data
leg_male_gam2 <- gam(lncount_leg_male ~ factor(year) +
                         te(longitude, latitude) +
                         s(julian) +
                         s(depth) +
                         s(phi) +
                         s(temperature) +
                         s(ice_index) +
                         s(longitude, latitude, by = legal_male_loading),
                       data = crab_train[crab_train$lncount_leg_male > 0, ])
summary(leg_male_gam2) # 50.9% explained

par(mfrow = c(2, 2))
gam.check(leg_male_gam2)

par(mfrow = c(3, 3))
plot(leg_male_gam2)

# Tweedie
# Base model
leg_male_tweedie <- gam(legal_male + 1 ~ factor(year) +
                            te(longitude, latitude) +
                            s(julian) +
                            s(depth),
                          data = crab_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(leg_male_tweedie) # 62.5% explained

# Add environmental data
leg_male_tweedie1 <- gam(legal_male + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(leg_male_tweedie1) # 65.2% explained

# Add obs data
leg_male_tweedie2 <- gam(legal_male + 1 ~ factor(year) +
                             te(longitude, latitude) +
                             s(julian) +
                             s(depth) +
                             s(phi) +
                             s(temperature) +
                             s(ice_index) +
                             s(longitude, latitude, by = legal_male_loading),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(leg_male_tweedie2) # 67.4%

# Add pcod data
leg_male_tweedie3 <- gam(legal_male + 1 ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_index) +
                           s(longitude, latitude, by = legal_male_loading) +
                           s(pcod_cpue),
                         data = crab_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(leg_male_tweedie3) # 67.4%

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
          ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
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
hist(crab_train$lncount_sub_male,
     main = "Sublegal Male Catch",
     xlab = "log(count+1)") 
hist(crab_train$lncount_sub_male[crab_train$lncount_sub_male > 0])

sub_male_gam_base <- gam(pres_sub_male ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian) +
                           s(depth),
                         data = crab_train,
                         family = "binomial")
summary(sub_male_gam_base) # 64.3% explained

# Add environmental data
sub_male_gam1 <- gam(lncount_sub_male ~ factor(year) + 
                       te(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_index),
                     data = crab_train[crab_train$lncount_sub_male > 0, ])
summary(sub_male_gam1) # 63.9% explained

# Add obs data
sub_male_gam2 <- gam(lncount_sub_male ~ factor(year) +
                       te(longitude, latitude) +
                       s(julian) +
                       s(depth) +
                       s(phi) +
                       s(temperature) +
                       s(ice_index) +
                       s(longitude, latitude, by = sublegal_male_loading),
                     data = crab_train[crab_train$lncount_sub_male > 0, ])
summary(sub_male_gam2) # 65.7% explained

par(mfrow = c(2, 2))
gam.check(sub_male_gam2)

par(mfrow = c(3, 3))
plot(sub_male_gam2)

# Tweedie
# Base model
sub_male_tweedie <- gam(immature_male + 1 ~ factor(year) +
                          te(longitude, latitude) +
                          s(julian) +
                          s(depth),
                        data = crab_train,
                        family = tw(link = "log"),
                        method = "REML")
summary(sub_male_tweedie) # 70.9% explained

# Add environmental data
sub_male_tweedie1 <- gam(immature_male + 1 ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_index),
                         data = crab_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie1) # 72.5% explained

# Add obs data
sub_male_tweedie2 <- gam(immature_male + 1 ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_index) +
                           s(longitude, latitude, by = sublegal_male_loading),
                         data = crab_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie2) # 74.7%

# Add pcod data
sub_male_tweedie3 <- gam(immature_male + 1 ~ factor(year) +
                           te(longitude, latitude) +
                           s(julian) +
                           s(depth) +
                           s(phi) +
                           s(temperature) +
                           s(ice_index) +
                           s(longitude, latitude, by = sublegal_male_loading) +
                           s(pcod_cpue),
                         data = crab_train,
                         family = tw(link = "log"),
                         method = "REML")
summary(sub_male_tweedie3) # 74.8%

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
          ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
          family = "serif",
          xlab = "Longitude",
          ylab = "Latitude",
          main = " ",
          cex.lab = 1.7,
          cex.axis =  1.7)
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

## Mature females ----
brt_mat_females1 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 27,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_mat_females1) 


brt_mat_females2 <- gbm.step(data = crab_train,
                         gbm.x = c(8:10, 19:23, 26),
                         gbm.y = 27,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females2)


brt_mat_females3 <- gbm.step(data = crab_train,
                         gbm.x = c(8:10, 19:23, 26),
                         gbm.y = 27,
                         family = 'gaussian',
                         tree.complexity = 3,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females3)


brt_mat_females4 <- gbm.step(data = crab_train,
                         gbm.x = c(8:10, 19:23, 26),
                         gbm.y = 27,
                         family = 'gaussian',
                         tree.complexity = 10,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females4)


brt_mat_females5 <- gbm.step(data = crab_train,
                         gbm.x = c(8:10, 19:23, 26),
                         gbm.y = 27,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.75)
summary(brt_mat_females5)


brt_mat_females6 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 27,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_mat_females6)

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

# par(mfrow = c(1, 3))
# gbm.perspec(females_final,
#             3,
#             z.range = c(0, 7.69),
#             theta = 60,
#             col = "light blue",
#             cex.axis = 0.8,
#             cex.lab = 1,
#             ticktype = "detailed")
# 
# gbm.perspec(females_final,
#             1, 3,
#             z.range = c(0, 10.75),
#             theta = 60,
#             col = "light blue",
#             cex.axis = 0.8,
#             cex.lab = 1,
#             ticktype = "detailed")
# 
# gbm.perspec(females_final,
#             1, 2,
#             z.range = c(0, 9.45),
#             theta = 60,
#             col = "light blue",
#             cex.axis = 0.8,
#             cex.lab = 1,
#             ticktype = "detailed")

# Training data map
nlat = 40
nlon = 60
latd = seq(min(crab_train$latitude), max(crab_train$latitude), length.out = nlat)
lond = seq(min(crab_train$longitude), max(crab_train$longitude), length.out = nlon)
spatial_grid_mat_female <- expand.grid(lond, latd)
names(spatial_grid_mat_female) <- c('longitude', 'latitude')
spatial_grid_mat_female$dist <- NA
for (k in 1:nrow(spatial_grid_mat_female)) {
  dist <-  distance_function(spatial_grid_mat_female$latitude[k],
                             spatial_grid_mat_female$longitude[k],
                             crab_train$latitude,
                             crab_train$longitude)
  spatial_grid_mat_female$dist[k] <- min(dist)
}
spatial_grid_mat_female$year <- 2010
spatial_grid_mat_female$depth <- median(crab_train$depth, na.rm = T)
spatial_grid_mat_female$phi <- median(crab_train$phi, na.rm = T)
spatial_grid_mat_female$julian <- median(crab_train$julian, na.rm = T)
spatial_grid_mat_female$temperature <- median(crab_train$temperature, na.rm = T)
spatial_grid_mat_female$ice_index <- median(crab_train$ice_index, na.rm = T)
spatial_grid_mat_female$female_loading <- median(crab_train$female_loading, na.rm = T) 

preds_mat_female <- predict.gbm(females_mat_final,
                                spatial_grid_mat_female,
                                n.trees = females_mat_final$gbm.call$best.trees, 
                                type = "response")
summary(preds_mat_female)

spatial_grid_mat_female$pred <- predict(females_mat_final,
                                        spatial_grid_mat_female,
                                        n.trees = females_mat_final$gbm.call$best.trees, 
                                        type = "response")



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

par(mfrow = c(1, 1),
    family = "serif")
image(lond,
      latd,
      t(matrix(spatial_grid_mat_female$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-181, -156),
      ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
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
      ylim = range(crab_train$latitude, na.rm = TRUE) + c(-.4, .5),
      main = "Distribution of Mature Females (BRT)",
      cex.main = 1.5,
      cex.lab = 1.5,
      cex.axis = 1.5)
symbols(crab_train$longitude[crab_train$lncount_mat_female > 0],
        crab_train$latitude[crab_train$lncount_mat_female > 0],
        circles = log(crab_train$lncount_mat_female + 1)[crab_train$lncount_mat_female > 0],
        inches = 0.1,
        bg = alpha('grey', 0.3),
        fg = alpha('black', 0.1),
        add = T)
points(crab_train$longitude[crab_train$lncount_mat_female == 0],
       crab_train$latitude[crab_train$lncount_mat_female == 0],
       pch =  '')
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = jet.colors(100),
           legend.shrink = 0.2,
           smallplot = c(.1, .13, .14, .29),
           legend.cex = 1,
           axis.args = list(cex.axis = 1.2),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(spatial_grid_mat_female$pred, na.rm = T), 
                    max(spatial_grid_mat_female$pred, na.rm = T)),
           legend.args = list("log(count + 1)",
                              side = 2, cex = 1.1))



## Immature females ----
brt_imm_females1 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_imm_females1) 


brt_imm_females2 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females2)


brt_imm_females3 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females3)


brt_imm_females4 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 10,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females4)


brt_imm_females5 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.75)
summary(brt_imm_females5)


brt_imm_females6 <- gbm.step(data = crab_train,
                             gbm.x = c(8:10, 19:23, 26),
                             gbm.y = 25,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_imm_females6)


# Attempt dropping variable
females_imm_simp <- gbm.simplify(brt_imm_females1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
females_imm_final <- brt_imm_females5 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_imm_final,
         n.plots = 9,
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

## Legal Males ----
base_brt_leg_males1 <- gbm.step(data = crab_train,
                                gbm.x = c(6, 19, 20),
                                gbm.y = 33,
                                grid)



brt_leg_males1 <- gbm.step(data = crab_train[crab_train$legal_male > 0, ],
                             gbm.x = c(8:10, 19:23, 21),
                             gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_leg_males1) 


brt_leg_males2 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 21),
                           gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males2)


brt_leg_males3 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 21),
                           gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males3)


brt_leg_males4 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 21),
                           gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 10,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males4)


brt_leg_males5 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 21),
                           gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.75)
summary(brt_leg_males5)


brt_leg_males6 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 21),
                           gbm.y = 26,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_leg_males6)


# Attempt dropping variable
males_leg_simp <- gbm.simplify(brt_leg_males1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
males_leg_final <- brt_leg_males5 # Change this once decision made

# Plot the variables
windows()
gbm.plot(males_leg_final,
         n.plots = 9,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
male_leg_effects <- tibble::as_tibble(summary.gbm(males_leg_final, plotit = F))
male_leg_effects %>% arrange(desc(rel.inf)) %>%
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
males_leg_int <- gbm.interactions(males_leg_final)
males_leg_int$interactions

## Sublegal Males ----
brt_sub_males1 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_sub_males1) 


brt_sub_males2 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males2)


brt_sub_males3 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 3,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males3)


brt_sub_males4 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 10,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males4)


brt_sub_males5 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.75)
summary(brt_sub_males5)


brt_sub_males6 <- gbm.step(data = crab_train,
                           gbm.x = c(8:10, 19:23, 22),
                           gbm.y = 27,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.25)
summary(brt_sub_males6)


# Attempt dropping variable
males_sub_simp <- gbm.simplify(brt_sub_males1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
males_sub_final <- brt_sub_males5 # Change this once decision made

# Plot the variables
windows()
gbm.plot(males_sub_final,
         n.plots = 9,
         plot.layout = c(3, 3),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

windows()
male_sub_effects <- tibble::as_tibble(summary.gbm(males_sub_final, plotit = F))
male_sub_effects %>% arrange(desc(rel.inf)) %>%
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
males_sub_int <- gbm.interactions(males_sub_final)
males_sub_int$interactions