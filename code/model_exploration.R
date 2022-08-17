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
source(here('code/functions', 'vis_gam_COLORS.R'))

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

crab_filtered <- crab_summary %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.)))

# Transform female and male data
crab_filtered$lncpue_female <- log(crab_filtered$female + 1)
crab_filtered$lncpue_male <- log(crab_filtered$male + 1)

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
crab_train <-as.data.frame(crab_filtered %>% 
                              filter(year < 2013))
crab_test <- as.data.frame(crab_filtered %>% 
                             filter(year > 2012))

# GAMs ----
# Female
# Gaussian
hist(crab_train$lncpue_female) # zero-inflated

female_gam_base <- gam(lncpue_female ~ s(latitude, longitude) +
                         s(doy) +
                         s(depth),
                       data = crab_train[crab_train$lncpue_female > 0, ])
summary(female_gam_base) # 28.6% explained

# Add environmental data
female_gam1 <- gam(lncpue_female ~ s(latitude, longitude) +
                     s(doy) +
                     s(phi) +
                     s(sst) +
                     s(ice) +
                     s(depth),
                   data = crab_train[crab_train$lncpue_female > 0, ])
summary(female_gam1) # 32.3% explained

# Add survey data
female_gam2 <- gam(lncpue_female ~ s(latitude, longitude) +
                     s(doy) +
                     s(phi) +
                     s(sst) +
                     s(ice) +
                     s(depth) +
                     s(female_immature) +
                     s(female_mature),
                   data = crab_train[crab_train$lncpue_female > 0, ])
summary(female_gam2) # 32.4% explained, added variables not significant
par(mfrow = c(2, 2))
gam.check(female_gam2)

# Tweedie
# Base model
female_tweedie <- gam(female + 1 ~ s(latitude, longitude) +
                        s(doy) +
                        s(depth),
                      data = crab_train,
                      family = tw(link = "log"),
                      method = "REML")
summary(female_tweedie) # 54.9% explained

# Add environmental data
female_tweedie1 <- gam(female + 1 ~ s(latitude, longitude) +
                         s(doy) +
                         s(phi) +
                         s(sst) +
                         s(ice) +
                         s(depth),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(female_tweedie1) # 56.6% explained

# Add survey data
female_tweedie2 <- gam(female + 1 ~ s(latitude, longitude) +
                         s(doy) +
                         s(phi) +
                         s(sst) +
                         s(ice) +
                         s(depth) +
                         s(female_immature) +
                         s(female_mature),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(female_tweedie2) # 56.9%

par(mfrow = c(2, 2))
gam.check(female_tweedie2)

par(mfrow = c(2, 4))
plot(female_tweedie2)


contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(female_tweedie2,
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
           zlim = c(min(female_tweedie2$linear.predictors),
                    max(female_tweedie2$linear.predictors)),
           legend.args = list("log(cpue+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

## ZIPLSS
female_ziplss <- gam(list(female ~ s(latitude, longitude) +
                            s(doy) +
                            s(phi) +
                            s(sst) +
                            s(ice) +
                            s(depth) +
                            s(female_immature) +
                            s(female_mature),
                          ~ s(latitude, longitude) +
                            s(doy)),
                     data = crab_train,
                     family = ziplss())
summary(female_ziplss) #33.8%

# Male
# Gaussian
hist(crab_train$lncpue_male) # left skewed 

male_gam_base <- gam(lncpue_male ~ s(latitude, longitude) +
                         s(doy) +
                         s(depth),
                       data = crab_train[crab_train$lncpue_male > 0, ])
summary(male_gam_base) # 13.9% explained

# Add environmental data
male_gam1 <- gam(lncpue_male ~ s(latitude, longitude) +
                     s(doy) +
                     s(phi) +
                     s(sst) +
                     s(ice) +
                     s(depth),
                   data = crab_train[crab_train$lncpue_male > 0, ])
summary(male_gam1) # 16.9% explained

# Add survey data
male_gam2 <- gam(lncpue_male ~ s(latitude, longitude) +
                     s(doy) +
                     s(phi) +
                     s(sst) +
                     s(ice) +
                     s(depth) +
                     s(male_immature) +
                     s(male_mature),
                   data = crab_train[crab_train$lncpue_male > 0, ])
summary(male_gam2) # 17.1% explained, added variables not significant
par(mfrow = c(2, 2))
gam.check(male_gam2)

# Tweedie
# Base model
male_tweedie <- gam(male + 1 ~ s(latitude, longitude) +
                        s(doy) +
                        s(depth),
                      data = crab_train,
                      family = tw(link = "log"),
                      method = "REML")
summary(male_tweedie) # 11.4% explained

# Add environmental data
male_tweedie1 <- gam(male + 1 ~ s(latitude, longitude) +
                         s(doy) +
                         s(phi) +
                         s(sst) +
                         s(ice) +
                         s(depth),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(male_tweedie1) # 14.5% explained

# Add survey data
male_tweedie2 <- gam(male + 1 ~ s(latitude, longitude) +
                         s(doy) +
                         s(phi) +
                         s(sst) +
                         s(ice) +
                         s(depth) +
                         s(male_immature) +
                         s(male_mature),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(male_tweedie2) # 14.5%

par(mfrow = c(2, 2))
gam.check(male_tweedie2)

par(mfrow = c(2, 4))
plot(male_tweedie2)


windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(male_tweedie2,
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
           zlim = c(min(male_tweedie2$linear.predictors),
                    max(male_tweedie2$linear.predictors)),
           legend.args = list("log(cpue+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

## ZIPLSS
male_ziplss <- gam(list(male ~ s(latitude, longitude) +
                            s(doy) +
                            s(phi) +
                            s(sst) +
                            s(ice) +
                            s(depth) +
                            s(male_immature) +
                            s(male_mature),
                          ~ s(latitude, longitude) +
                            s(doy)),
                     data = crab_train,
                     family = ziplss())
summary(male_ziplss) # 16.7%

# Random forests ----
# Females
set.seed(1993)
rf_females <- randomForest(lncpue_female ~ latitude +
                             longitude +
                             doy,
                            data = na.exclude(crab_train), # throws error if NAs included
                            ntree = 1000,
                            mtry = 2,
                            importance = T,
                            proximity = T)
print(rf_females) 
# 70.99% variance explained
# Mean squared residuals: 0.18

plot(rf_females$mse)



set.seed(1993)
rf_females1 <- randomForest(lncpue_female ~ depth +
                              latitude +
                              longitude +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature,
                            data = na.exclude(crab_train), 
                            ntree = 1000,
                            mtry = 2,
                            importance = T,
                            proximity = T)
print(rf_females1) 
# 70.94% variance explained
# Mean squared residuals: 0.18

plot(rf_females1$mse)

rf_females1$importance
varImpPlot(rf_females1)


set.seed(1993)
rf_females2 <- randomForest(lncpue_female ~ depth +
                              latitude +
                              longitude +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature +
                              bottom_temp,
                            data = na.exclude(crab_train), 
                            ntree = 1000,
                            mtry = 2,
                            importance = T,
                            proximity = T)
print(rf_females2) 
# 69.85% variance explained
# Mean squared residuals: 0.18

plot(rf_females2$mse) # adding bottom temp made it worse

rf_females2$importance
varImpPlot(rf_females2)


set.seed(1993)
rf_females3 <- randomForest(lncpue_female ~ depth +
                              latitude +
                              longitude +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature,
                            data = na.exclude(crab_train), 
                            ntree = 3000,
                            mtry = 2,
                            importance = T,
                            proximity = T)
print(rf_females3) 
# 71.12% variance explained
# Mean squared residuals: 0.18

plot(rf_females3$mse) 


set.seed(1993)
rf_females4 <- randomForest(lncpue_female ~ depth +
                              latitude +
                              longitude +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature,
                            data = na.exclude(crab_train), 
                            ntree = 3000,
                            mtry = 5,
                            importance = T,
                            proximity = T)
print(rf_females4) 
# 70.7% variance explained
# Mean squared residuals: 0.18

plot(rf_females4$mse) 



set.seed(1993)
rf_females5 <- randomForest(lncpue_female ~ depth +
                              latitude +
                              longitude +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature,
                            data = na.exclude(crab_train), 
                            ntree = 3000,
                            mtry = 3,
                            importance = T,
                            proximity = T)
print(rf_females5) 
# 71.2% variance explained
# Mean squared residuals: 0.18

plot(rf_females5$mse) 

# Males
set.seed(1993)
rf_males <- randomForest(lncpue_male ~ latitude +
                           longitude +
                           doy,
                          data = na.exclude(crab_train),
                          ntree = 1000,
                          mtry = 2,
                          importance = T,
                          proximity = T)
print(rf_males)
# 56.96% variance explained
# Mean squared residuals: 0.42

plot(rf_males$mse)


set.seed(1993)
rf_males1 <- randomForest(lncpue_male ~ depth +
                            latitude +
                            longitude +
                            phi +
                            ice +
                            sst +
                            doy +
                            male_mature +
                            male_immature,
                          data = na.exclude(crab_train),
                          ntree = 1000,
                          mtry = 5,
                          importance = T,
                          proximity = T)
print(rf_males1)
# 56.96% variance explained
# Mean squared residuals: 0.42

plot(rf_males1$mse)

rf_males1$importance
varImpPlot(rf_males1)


# Boosted regression trees ----
# Adjust the bag fraction to a value between 0.5-0.75 as suggested by Elith et al. (2008)
# The learning rate could range from 0.1-0.0001, higher value usually means less trees
# Depending on the number of samples, want tree complexity to be high enough (likely using 5)
# Want at least 1000 trees, but don't need to go way beyond it

# Females
brt_females1 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.05,
                         bag.fraction = 0.5) 
summary(brt_females1) 
# 4250 trees
# dev: 0.353
# res dev: 0.107
# corr: 0.838
# cv corr: 0.7

brt_females2 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females2)
# 8300 trees
# dev: 0.353
# res dev: 0.139
# corr: 0.783
# cv corr: 0.688

brt_females3 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 3,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females3)
# 10000 trees
# dev: 0.353
# res dev: 0.164
# corr: 0.735
# cv corr: 0671

brt_females4 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 10,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females4)
# 6800 trees
# dev: 0.353
# res dev: 0.11
# corr: 0.835
# cv corr: 0.703

brt_females5 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:31, 34:37),
                       gbm.y = 38,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.1,
                       bag.fraction = 0.5)
summary(brt_females5)

brt_females6 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:31, 34:37),
                       gbm.y = 38,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.1,
                       bag.fraction = 0.75)
summary(brt_females6)

brt_females7 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.05,
                         bag.fraction = 0.75)
summary(brt_females7)

# Attempt dropping variable
females_simp <- gbm.simplify(brt_females1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
females_final <- brt_females1 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_final,
         n.plots = 7,
         plot.layout = c(4, 2),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

# Plot the fits
females_int <- gbm.interactions(females_final)
females_int$interactions

# par(mfrow = c(1, 3))
# gbm.perspec(females_final,
#             2, 3,
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


# Males
brt_males1 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.05,
                       bag.fraction = 0.5)
summary(brt_males1)
# 10000 trees
# dev: 1.706
# res dev: 0.492
# corr: 0.851
# cv corr: 0.693

brt_males2 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.01,
                       bag.fraction = 0.5)
summary(brt_males2)

brt_males3 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 3,
                       learning.rate = 0.01,
                       bag.fraction = 0.5)
summary(brt_males3)

brt_males4 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 10,
                       learning.rate = 0.01,
                       bag.fraction = 0.5)
summary(brt_males4)

brt_males5 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 10,
                       learning.rate = 0.01,
                       bag.fraction = 0.75)
summary(brt_males5)

brt_males6 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.1,
                       bag.fraction = 0.5)
summary(brt_males6)

brt_males7 <- gbm.step(data = crab_train,
                       gbm.x = c(8, 25, 28:30, 32, 33, 35:37),
                       gbm.y = 39,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.1,
                       bag.fraction = 0.75)
summary(brt_males7)

# Attempt dropping variable
males_simp <- gbm.simplify(brt_males6, n.drops = 5) # this takes forever
summary(males_simp)

# Choose final model
males_final <- brt_males1 # Change this once decision made

# Plot the variables
windows()
gbm.plot(males_final,
         n.plots = 7,
         plot.layout = c(4, 2),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

# Plot the fits
males_int <- gbm.interactions(males_final)
males_int$interactions