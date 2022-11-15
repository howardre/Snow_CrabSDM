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
observer_filtered <- select(observer_df, lncpue_male_legal, lncpue_male_sublegal, lncpue_female)

# Create matrices
corrplot(cor(observer_filtered), is.corr = FALSE)
observer_scaled <- scale(observer_filtered)

# legal_male <- as.data.frame(legal_male_df[order(legal_male_df$index), ])
# legal_male_mat <- fossil::create.matrix(legal_male,
#                                         tax.name = "index",
#                                         locality = "doy",
#                                         time.col = NULL,
#                                         time = NULL,
#                                         abund = T,
#                                         abund.col = "anomalies")

# PCA
pca <- prcomp(observer_scaled,
              scale = F,
              center = F)
summary(pca)

pca$sdev # eigenvalues
pca$x[, 1][1:10] # first 10 values of pca 1
pca$rotation[, 1] # loadings for first pca

# Checks
pca1_1 <- pca$rotation[, 1]%*%observer_scaled[1, ]
pca1_1;pca$x[1, ][1]

# Plots
plot(pca$sdev, 
     type = 'h',
     ylab = 'Standard Deviation',
     xlab = 'PC')
plot(pca$x[, 1],
     xlab = 'Values',
     ylab = 'PC1',
     type = 'b')
plot(pca$rotation[, 1],
     ylab = 'Loadings',
     xlab = 'Life Stages')

# Try with vegan
vegan_pca <- princomp(observer_scaled, cor = TRUE)
summary(vegan_pca)
head(summary(vegan_pca))

vegan_pca$sdev
vegan_pca$loadings[, 1]
vegan_pca$scores[, 1]
vegan_pca$scale
vegan_pca$center

# Plot
biplot(vegan_pca)
