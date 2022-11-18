# Title: Observer Data PCA
# Purpose: Determine if PCA can be used for observer data
# Author:  Rebecca Howard

# Load libraries ----
library(here)
library(vegan)
library(dplyr)

# Functions ----

# Load data ----
observer_df <- readRDS(here('data/Snow_CrabData', 'observer_summarized.rds'))
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

# Filter to groups
observer_df$lncpue_male_legal <- log(observer_df$obs_male_legal + 1)
observer_df$lncpue_male_sublegal <- log(observer_df$obs_male_sub + 1)
observer_df$lncpue_female <- log(observer_df$obs_female + 1)

# Average values for the stations and years

# Make matrices with year and station ID for each group
# Sites as rows, years as columns?
legal_male_mat <- fossil::create.matrix(as.data.frame(observer_df),
                                        tax.name = "STATIONID",
                                        locality = "year_lag",
                                        time.col = NULL,
                                        time = NULL,
                                        abund = T,
                                        abund.col = "lncpue_male_legal")

sublegal_male_mat <- fossil::create.matrix(as.data.frame(observer_df),
                                           tax.name = "STATIONID",
                                           locality = "year_lag",
                                           time.col = NULL,
                                           time = NULL,
                                           abund = T,
                                           abund.col = "lncpue_male_sublegal")

female_mat <- fossil::create.matrix(as.data.frame(observer_df),
                                        tax.name = "STATIONID",
                                        locality = "year_lag",
                                        time.col = NULL,
                                        time = NULL,
                                        abund = T,
                                        abund.col = "lncpue_female")

# Create matrices
corrplot(cor(legal_male_mat), is.corr = FALSE)
corrplot(cor(sublegal_male_mat), is.corr = FALSE)
corrplot(cor(female_mat), is.corr = FALSE)

# Scale values
legal_male_scaled <- scale(legal_male_mat)
sublegal_male_scaled <- scale(sublegal_male_mat)
female_scaled <- scale(female_mat)

# PCA
# prcomp version, using vegan below
# legal_male_pca <- prcomp(legal_male_scaled,
#                          scale = F,
#                          center = F)
# summary(legal_male_pca)
# 
# legal_male_pca$sdev # eigenvalues
# legal_male_pca$x[, 1][1:10] # first 10 values of pca 1
# legal_male_pca$rotation[, 1] # loadings for first pca
# 
# # Checks
# legal_male_pca1_1 <- legal_male_pca$rotation[, 1]%*%legal_male_scaled[1, ]
# legal_male_pca1_1;legal_male_pca$x[1, ][1]
# 
# # Plots
# plot(legal_male_pca$sdev, 
#      type = 'h',
#      ylab = 'Standard Deviation',
#      xlab = 'PC')
# plot(legal_male_pca$x[, 1],
#      xlab = 'Values',
#      ylab = 'PC1',
#      type = 'b')
# plot(legal_male_pca$rotation[, 1],
#      ylab = 'Loadings',
#      xlab = 'Values')

# Vegan version
legal_male_pca <- princomp(legal_male_scaled, cor = TRUE)
summary(legal_male_pca)
head(summary(legal_male_pca))

legal_male_pca$sdev
legal_male_pca$loadings[, 1] # Use these in the model
legal_male_pca$scores[, 1]
legal_male_pca$scale
legal_male_pca$center

sublegal_male_pca <- princomp(sublegal_male_scaled, cor = TRUE)
summary(sublegal_male_pca)
head(summary(sublegal_male_pca))

sublegal_male_pca$sdev
sublegal_male_pca$loadings[, 1] # Use these in the model
sublegal_male_pca$scores[, 1]
sublegal_male_pca$scale
sublegal_male_pca$center

female_pca <- princomp(female_scaled, cor = TRUE)
summary(female_pca)
head(summary(female_pca))

female_pca$sdev
female_pca$loadings[, 1] # Use these in the model
female_pca$scores[, 1]
female_pca$scale
female_pca$center

# Plot
biplot(legal_male_pca)
biplot(sublegal_male_pca)
biplot(female_pca)

## Match to data frame
legal_male_loadings <- as.data.frame(legal_male_pca$loadings[, 1])
legal_male_loadings <- tibble::rownames_to_column(legal_male_loadings, "year")
colnames(legal_male_loadings)[2] <- "legal_male_loadings"

plot(legal_male_loadings,
     main = "Legal Male Loadings",
     ylab = "value")

sublegal_male_loadings <- as.data.frame(sublegal_male_pca$loadings[, 1])
sublegal_male_loadings <- tibble::rownames_to_column(sublegal_male_loadings, "year")
colnames(sublegal_male_loadings)[2] <- "sublegal_male_loadings"

plot(sublegal_male_loadings,
     main = "Sublegal Male Loadings",
     ylab = "value")

female_loadings <- as.data.frame(female_pca$loadings[, 1])
female_loadings <- tibble::rownames_to_column(female_loadings, "year")
colnames(female_loadings)[2] <- "female_loadings"

plot(female_loadings,
     main = "Female Loadings",
     ylab = "value")

crab_summary$legal_male_loading <- legal_male_loadings$legal_male_loadings[match(crab_summary$year, legal_male_loadings$year)]
crab_summary$sublegal_male_loading <- sublegal_male_loadings$sublegal_male_loadings[match(crab_summary$year, sublegal_male_loadings$year)]
crab_summary$female_loading <- female_loadings$female_loadings[match(crab_summary$year, female_loadings$year)]

# Save data
saveRDS(crab_summary, file = here('data/Snow_CrabData', 'crab_pca.rds'))
