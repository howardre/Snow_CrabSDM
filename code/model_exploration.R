# Title: Snow Crab Model Exploration
# Purpose: Investigate potential models
# Data created: 07/20/2022

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(randomForest)
library(here)

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

crab_filtered <- crab_summary %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.)))

# Transform female and male data
crab_filtered$lncpue_female <- log(crab_filtered$female + 1)
crab_filtered$lncpue_male <- log(crab_filtered$male + 1)

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
crab_train <- as.data.frame(crab_filtered %>% 
                              filter(year < 2013))
crab_test <- as.data.frame(crab_filtered %>% 
                             filter(year > 2012))

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

brt_females2 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females2)

brt_females3 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 3,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females3)

brt_females4 <- gbm.step(data = crab_train,
                         gbm.x = c(8, 25, 28:31, 34:37),
                         gbm.y = 38,
                         family = 'gaussian',
                         tree.complexity = 10,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_females4)

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
males_final <- brt_males6 # Change this once decision made

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