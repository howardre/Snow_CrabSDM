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


# Boosted regression trees ----
# Females
brt_females1 <- gbm.step(data = crab_train,
                         gbm.x = c(9, 26, 29:32, 35),
                         gbm.y = 36,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5)
summary(brt_females1)
