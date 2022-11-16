# Title: Observer Data PCA
# Purpose: Determine if PCA can be used for observer data
# Author:  Rebecca Howard

# Load libraries ----
library(here)
library(vegan)
library(dplyr)

# Functions ----

# Load data ----
observer_df <- readRDS(here('data/Snow_CrabData', 'observer_df.rds'))
bycatch_df <- readRDS(here('data/Snow_CrabData', 'bycatch_df.rds'))

# Filter to groups
observer_df$index <- 1:nrow(observer_df)
observer_df$lncpue_male_legal <- log(observer_df$tot_legal + 1)
observer_df$lncpue_male_sublegal <- log(observer_df$sublegal + 1)
observer_df$lncpue_female <- log(observer_df$female + 1)

# Make matrices with year and station ID for each group
# Sites as rows, years as columns?
legal_male_mat <- fossil::create.matrix(as.data.frame(observer_df),
                                        tax.name = "STATIONID",
                                        locality = "year",
                                        time.col = NULL,
                                        time = NULL,
                                        abund = T,
                                        abund.col = "lncpue_male_legal")

# Create matrices
corrplot(cor(legal_male_mat), is.corr = FALSE)
legal_male_scaled <- scale(legal_male_mat)



# PCA
legal_male_pca <- prcomp(legal_male_scaled,
                         scale = F,
                         center = F)
summary(legal_male_pca)

legal_male_pca$sdev # eigenvalues
legal_male_pca$x[, 1][1:10] # first 10 values of pca 1
legal_male_pca$rotation[, 1] # loadings for first pca

# Checks
legal_male_pca1_1 <- legal_male_pca$rotation[, 1]%*%legal_male_scaled[1, ]
legal_male_pca1_1;legal_male_pca$x[1, ][1]

# Plots
plot(legal_male_pca$sdev, 
     type = 'h',
     ylab = 'Standard Deviation',
     xlab = 'PC')
plot(legal_male_pca$x[, 1],
     xlab = 'Values',
     ylab = 'PC1',
     type = 'b')
plot(legal_male_pca$rotation[, 1],
     ylab = 'Loadings',
     xlab = 'Values')

# Try with vegan
legal_male_pca <- princomp(legal_male_scaled, cor = TRUE)
summary(legal_male_pca)
head(summary(legal_male_pca))

legal_male_pca$sdev
legal_male_pca$loadings[, 1] # Use these in the model
legal_male_pca$scores[, 1]
legal_male_pca$scale
legal_male_pca$center

# Plot
biplot(legal_male_pca)

## Match to data frame

