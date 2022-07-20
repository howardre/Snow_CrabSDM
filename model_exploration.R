# Title: Snow Crab Model Exploration
# Purpose: Investigate potential models
# Data created: 07/20/2022

# Load libraries ----
library(gbm)
library(dismo)
library(scales)
library(randomForest)

# Load data ----
crab_filtered <- crab_summary %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.)))

# Response variable has large values, which leads to issues with gbm.step() so it was rescaled
crab_filtered$scaled_female <- rescale(crab_filtered$female)
crab_filtered$scaled_male <- rescale(crab_filtered$male)

# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
crab_train <- as.data.frame(crab_filtered %>% 
                              filter(year < 2016))
crab_test <- as.data.frame(crab_filtered %>% 
                             filter(year > 2015))

# Random forests ----
# Females


# Boosted regression trees ----
# Females
brt_females1 <- gbm.step(data = crab_train,
                         gbm.x = c(9, 10, 26),
                         gbm.y = 29,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5)
