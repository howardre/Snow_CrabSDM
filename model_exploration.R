# Title: Snow Crab Model Exploration
# Purpose: Investigate potential models
# Data created: 07/20/2022

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(randomForest)

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

crab_filtered <- crab_summary %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.)))

# Response variable has large values, which leads to issues with gbm.step() so it was rescaled
crab_filtered$scaled_female <- scales::rescale(crab_filtered$female)
crab_filtered$scaled_male <- scales::rescale(crab_filtered$male)

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
crab_train <- as.data.frame(crab_filtered %>% 
                              filter(year < 2016))
crab_test <- as.data.frame(crab_filtered %>% 
                             filter(year > 2015))

# Random forests ----
# Females
set.seed(1993)
rf_females1 <- randomForest(scaled_female ~ depth +
                              phi +
                              ice +
                              sst +
                              doy +
                              female_mature +
                              female_immature,
                            data = na.exclude(crab_train), # throws error if NAs included
                            ntree = 1000,
                            mtry = 5,
                            importance = T,
                            proximity = T)
print(rf_females1)

# Males
set.seed(1993)
rf_males1 <- randomForest(scaled_male ~ depth +
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
print(rf_females1)


# Boosted regression trees ----
# Females
brt_females1 <- gbm.step(data = crab_train,
                         gbm.x = c(9, 26, 29:32, 35),
                         gbm.y = 36,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.05,
                         bag.fraction = 0.5) # probably producing too many trees
summary(brt_females1) 

brt_females2 <- gbm.step(data = crab_train,
                         gbm.x = c(9, 26, 29:32, 35),
                         gbm.y = 36,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) # 
summary(brt_females2)

brt_females3 <- gbm.step(data = crab_train,
                         gbm.x = c(9, 26, 29:32, 35),
                         gbm.y = 36,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) # 
summary(brt_females3)

# Attempt dropping variable
females_simp <- gbm.simplify(brt_females1, n.drops = 5)
summary(females_simp)

# Choose final model
females_final <- females_simp # Change this once decision made

# Plot the variables
gbm.plot(females_final,
         n.plots = 7,
         plot.layout = c(3, 4),
         write.title = F,
         smooth = T,
         common.scale = T,
         cex.axis = 1.7,
         cex.lab = 1.7,
         lwd = 1.5)

# Plot the fits
females_int <- gbm.interactions(females_final)
females_int$interactions

par(mfrow = c(1, 3))
gbm.perspec(females_final,
            2, 3,
            z.range = c(0, 7.69),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_final,
            1, 3,
            z.range = c(0, 10.75),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")

gbm.perspec(females_final,
            1, 2,
            z.range = c(0, 9.45),
            theta = 60,
            col = "light blue",
            cex.axis = 0.8,
            cex.lab = 1,
            ticktype = "detailed")


# Males
brt_males1 <- gbm.step(data = crab_train,
                       gbm.x = c(9, 26, 29:31, 33:34),
                       gbm.y = 37,
                       family = 'gaussian',
                       tree.complexity = 5,
                       learning.rate = 0.05,
                       bag.fraction = 0.5)
summary(brt_males1)
