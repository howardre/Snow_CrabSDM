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

# Load data ----
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncpue_mat_female = log(mature_female + 1),
                     lncpue_imm_female = log(immature_female + 1),
                     lncpue_leg_male = log(legal_male + 1),
                     lncpue_sub_male = log(immature_male + 1),
                     lncpue_obs_female = log(obs_female + 1),
                     lncpue_obs_sub_male = log(obs_male_sub + 1),
                     lncpue_obs_leg_male = log(obs_male_legal + 1))


# Create train and test datasets
# Considering using blocked approach but current discussion pointed toward using certain years
# Need to filter out stations without observer data in order to get accurate comparison
crab_final <- crab_trans[complete.cases(crab_trans), ]
crab_train <- as.data.frame(crab_final %>% 
                              filter(year < 2015))
crab_test <- as.data.frame(crab_final %>% 
                             filter(year > 2014))

# GAMs ----
# Mature Female
# Gaussian
hist(crab_train$lncpue_mat_female) # zero-inflated
hist(crab_train$lncpue_mat_female[crab_train$lncpue_mat_female > 0])

mat_female_gam_base <- gam(lncpue_mat_female ~ factor(year) +
                             s(latitude, longitude) +
                             s(julian) +
                             s(depth),
                           data = crab_train[crab_train$lncpue_mat_female > 0, ])
summary(mat_female_gam_base) # 34.1% explained

# Add environmental data
mat_female_gam1 <- gam(lncpue_mat_female ~ factor(year) + 
                         s(latitude, longitude) +
                         s(julian) +
                         s(phi) +
                         s(temperature) +
                         s(depth),
                       data = crab_train[crab_train$lncpue_mat_female > 0, ])
summary(mat_female_gam1) # 32.9% explained

# Add survey data
mat_female_gam2 <- gam(lncpue_mat_female ~ factor(year) +
                         s(latitude, longitude) +
                         s(julian) +
                         s(phi) +
                         s(temperature) +
                         s(depth) +
                         s(lncpue_obs_female),
                       data = crab_train[crab_train$lncpue_mat_female > 0, ])
summary(mat_female_gam2) # 46.4

par(mfrow = c(2, 2))
gam.check(mat_female_gam2)

# Tweedie
# Base model
mat_female_tweedie <- gam(mature_female + 1 ~ factor(year) +
                            s(latitude, longitude) +
                            s(julian) +
                            s(depth),
                          data = crab_train,
                          family = tw(link = "log"),
                          method = "REML")
summary(mat_female_tweedie) # 49.9% explained

# Add environmental data
mat_female_tweedie1 <- gam(mature_female + 1 ~ factor(year) +
                             s(latitude, longitude) +
                             s(julian) +
                             s(phi) +
                             s(temperature) +
                             s(depth),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie1) # 51.9% explained

# Add survey data
mat_female_tweedie2 <- gam(mature_female + 1 ~ factor(year) +
                             s(latitude, longitude) +
                             s(julian) +
                             s(phi) +
                             s(temperature) +
                             s(depth) +
                             s(lncpue_obs_female),
                           data = crab_train,
                           family = tw(link = "log"),
                           method = "REML")
summary(mat_female_tweedie2) # 

par(mfrow = c(2, 2))
gam.check(mat_female_tweedie2)

par(mfrow = c(3, 2))
plot(mat_female_tweedie2)


contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

windows(width = 12, height = 10)
par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
myvis_gam(mat_female_tweedie2,
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
           zlim = c(min(mat_female_tweedie2$linear.predictors),
                    max(mat_female_tweedie2$linear.predictors)),
           legend.args = list("log(cpue+1)",
                              side = 2,
                              cex = 1.5,
                              family =  "serif"))

## ZIPLSS
mat_female_ziplss <- gam(list(mature_female ~ factor(year) +
                                s(latitude, longitude) +
                                s(julian) +
                                s(phi) +
                                s(temperature) +
                                s(depth) +
                                s(lncpue_obs_female),
                              ~ factor(year) +
                                s(latitude, longitude) +
                                s(julian)),
                         data = crab_train,
                         family = ziplss())
summary(female_ziplss) # not working because counts are needed

# Male
# Gaussian
hist(crab_train$lncpue_leg_male) # left skewed 

leg_male_gam_base <- gam(lncpue_leg_male ~ s(latitude, longitude) +
                         s(julian) +
                         s(depth),
                       data = crab_train[crab_train$lncpue_leg_male > 0, ])
summary(leg_male_gam_base) # 39.6% explained

# Add environmental data
leg_male_gam1 <- gam(lncpue_leg_male ~ s(latitude, longitude) +
                     s(julian) +
                     s(phi) +
                     s(temperature) +
                     
                     s(depth),
                   data = crab_train[crab_train$lncpue_leg_male > 0, ])
summary(leg_male_gam1) # 41.3% explained

# Add survey data
leg_male_gam2 <- gam(lncpue_leg_male ~ s(latitude, longitude) +
                     s(julian) +
                     s(phi) +
                     s(temperature) +
                     
                     s(depth) +
                     s(obs_male_legal),
                   data = crab_train[crab_train$lncpue_leg_male > 0, ])
summary(leg_male_gam2) # 
par(mfrow = c(2, 2))
gam.check(leg_male_gam2)

# Tweedie
# Base model
leg_male_tweedie <- gam(legal_male + 1 ~ s(latitude, longitude) +
                        s(julian) +
                        s(depth),
                      data = crab_train,
                      family = tw(link = "log"),
                      method = "REML")
summary(leg_male_tweedie) # 55% explained

# Add environmental data
leg_male_tweedie1 <- gam(legal_male + 1 ~ s(latitude, longitude) +
                         s(julian) +
                         s(phi) +
                         s(temperature) +
                         
                         s(depth),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(leg_male_tweedie1) # 56% explained

# Add survey data
leg_male_tweedie2 <- gam(legal_male + 1 ~ s(latitude, longitude) +
                         s(julian) +
                         s(phi) +
                         s(temperature) +
                         
                         s(depth) +
                         s(obs_male_legal),
                       data = crab_train,
                       family = tw(link = "log"),
                       method = "REML")
summary(leg_male_tweedie2) # 14.5%

par(mfrow = c(2, 2))
gam.check(leg_male_tweedie2)

par(mfrow = c(2, 4))
plot(leg_male_tweedie2)


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
                            s(julian) +
                            s(phi) +
                            s(temperature) +
                            
                            s(depth) +
                            s(male_immature) +
                            s(male_mature),
                          ~ s(latitude, longitude) +
                            s(julian)),
                     data = crab_train,
                     family = ziplss())
summary(male_ziplss) # 16.7%

# Random forests ----
# Females
set.seed(1993)
rf_females <- randomForest(lncpue_female ~ latitude +
                             longitude +
                             julian,
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
                              temperature +
                              julian +
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
                              temperature +
                              julian +
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
                              temperature +
                              julian +
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
                              temperature +
                              julian +
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
                              temperature +
                              julian +
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
                           julian,
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
                            temperature +
                            julian +
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

# Mature females
brt_mat_females1 <- gbm.step(data = crab_train,
                         gbm.x = c(1, 3:5, 14, 16, 17, 23),
                         gbm.y = 19,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.05,
                         bag.fraction = 0.5) 
summary(brt_mat_females1) 
# 550 trees
# dev: 17.746
# res dev: 3.637
# corr: 0.897
# cv corr: 0.719

brt_mat_females2 <- gbm.step(data = crab_train,
                         gbm.x = c(1, 3:5, 14, 16, 17, 23),
                         gbm.y = 19,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females2)
# 1900 trees
# dev: 17.746
# res dev: 4.342
# corr: 0.875
# cv corr: 0.715

brt_mat_females3 <- gbm.step(data = crab_train,
                         gbm.x = c(1, 3:5, 14, 16, 17, 23),
                         gbm.y = 19,
                         family = 'gaussian',
                         tree.complexity = 3,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females3)
# 3250 trees
# dev: 17.746
# res dev: 4.724
# corr: 0.862
# cv corr: 0.718

brt_mat_females4 <- gbm.step(data = crab_train,
                         gbm.x = c(1, 3:5, 14, 16, 17, 23),
                         gbm.y = 19,
                         family = 'gaussian',
                         tree.complexity = 10,
                         learning.rate = 0.01,
                         bag.fraction = 0.5) 
summary(brt_mat_females4)
# 700 trees
# dev: 17.746
# res dev: 4.612
# corr: 0.867
# cv corr: 0.722

brt_mat_females5 <- gbm.step(data = crab_train,
                         gbm.x = c(1, 3:5, 14, 16, 17, 23),
                         gbm.y = 19,
                         family = 'gaussian',
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         bag.fraction = 0.75)
summary(brt_mat_females5)
# 1300 trees
# dev: 17.746
# res dev: 5.041
# corr: 0.852
# cv corr: 0.713

brt_mat_females6 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 19,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_mat_females5)
# 1300 trees
# dev: 17.746
# res dev: 4.905
# corr: 0.855
# cv corr: 0.719

# Attempt dropping variable
females_mat_simp <- gbm.simplify(brt_mat_females1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
females_mat_final <- brt_mat_females2 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_mat_final,
         n.plots = 7,
         plot.layout = c(4, 2),
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

# Immature females
brt_imm_females1 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_imm_females1) 
# 200 trees
# dev: 
# res dev: 
# corr: 
# cv corr: 

brt_imm_females2 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females2)
# 750 trees
# dev: 5.505
# res dev: 3.102
# corr: 0.693
# cv corr: 0.425

brt_imm_females3 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females3)
# 1000 trees
# dev: 5.505
# res dev: 3.46
# corr: 0.637
# cv corr: 0.418

brt_imm_females4 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 10,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_imm_females4)
# 450 trees
# dev: 5.505
# res dev: 2.856
# corr: 0.734
# cv corr: 0.44

brt_imm_females5 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.75)
summary(brt_imm_females5)
# 850 trees
# dev: 5.505
# res dev: 2.992
# corr: 0.713
# cv corr: 0.421

brt_imm_females6 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 23),
                             gbm.y = 20,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_imm_females5)
# 900 trees
# dev: 5.505
# res dev: 3.184
# corr: 0.674
# cv corr: 0.417

# Attempt dropping variable
females_imm_simp <- gbm.simplify(brt_imm_females1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
females_imm_final <- brt_imm_females5 # Change this once decision made

# Plot the variables
windows()
gbm.plot(females_imm_final,
         n.plots = 7,
         plot.layout = c(4, 2),
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

# Legal Males
brt_leg_males1 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.05,
                             bag.fraction = 0.5) 
summary(brt_leg_males1) 
# 200 trees
# dev: 
# res dev: 
# corr: 
# cv corr: 

brt_leg_males2 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males2)
# 1000 trees
# dev: 9.102
# res dev: 2.747
# corr: 0.841
# cv corr: 0.716

brt_leg_males3 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males3)
# 1350 trees
# dev: 9.102
# res dev: 3.147
# corr: 0.813
# cv corr: 0.708

brt_leg_males4 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 10,
                             learning.rate = 0.01,
                             bag.fraction = 0.5) 
summary(brt_leg_males4)
# 700 trees
# dev: 9.102
# res dev: 2.205
# corr: 0.877
# cv corr: 0.724

brt_leg_males5 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.75)
summary(brt_leg_males5)
# 950 trees
# dev: 9.102
# res dev: 2.744
# corr: 0.842
# cv corr: 0.719

brt_leg_males6 <- gbm.step(data = crab_train,
                             gbm.x = c(1, 3:5, 14, 16, 17, 25),
                             gbm.y = 21,
                             family = 'gaussian',
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             bag.fraction = 0.25)
summary(brt_leg_males5)
# 1650 trees
# dev: 9.102
# res dev: 2.5
# corr: 0.856
# cv corr: 0.719

# Attempt dropping variable
males_leg_simp <- gbm.simplify(brt_leg_males1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
males_leg_final <- brt_leg_males6 # Change this once decision made

# Plot the variables
windows()
gbm.plot(males_leg_final,
         n.plots = 7,
         plot.layout = c(4, 2),
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

# Sublegal Males
brt_sub_males1 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.05,
                           bag.fraction = 0.5) 
summary(brt_sub_males1) 
# 300 trees
# dev: 
# res dev: 
# corr: 
# cv corr: 

brt_sub_males2 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males2)
# 1050 trees
# dev: 9.089
# res dev: 3.029
# corr: 0.822
# cv corr: 0.686

brt_sub_males3 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 3,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males3)
# 1500 trees
# dev: 9.089
# res dev: 3.332
# corr: 0.8
# cv corr: 0.686

brt_sub_males4 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 10,
                           learning.rate = 0.01,
                           bag.fraction = 0.5) 
summary(brt_sub_males4)
# 750 trees
# dev: 9.089
# res dev: 2.448
# corr: 0.862
# cv corr: 0.69

brt_sub_males5 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.75)
summary(brt_sub_males5)
# 1100 trees
# dev: 9.089
# res dev: 2.921
# corr: 0.83
# cv corr: 0.686

brt_sub_males6 <- gbm.step(data = crab_train,
                           gbm.x = c(1, 3:5, 14, 16, 17, 24),
                           gbm.y = 22,
                           family = 'gaussian',
                           tree.complexity = 5,
                           learning.rate = 0.01,
                           bag.fraction = 0.25)
summary(brt_sub_males5)
# 1400 trees
# dev: 9.089
# res dev: 2.909
# corr: 0.83
# cv corr: 0.693

# Attempt dropping variable
males_sub_simp <- gbm.simplify(brt_sub_males1, n.drops = 5) # this takes forever
summary(females_simp)

# Choose final model
males_sub_final <- brt_sub_males4 # Change this once decision made

# Plot the variables
windows()
gbm.plot(males_sub_final,
         n.plots = 7,
         plot.layout = c(4, 2),
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