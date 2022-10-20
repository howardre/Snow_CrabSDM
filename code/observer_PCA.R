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
legal_male_df <- na.omit(select(observer_df, year, doy, tot_legal, latitude, longitude, depth))
sublegal_male_df <- na.omit(select(observer_df, year, doy, sublegal, latitude, longitude, depth))
female_df <- na.omit(select(observer_df, year, doy, female, latitude, longitude, depth))

# Add index for haul ID
legal_male_df$index <- 1:nrow(legal_male_df)
legal_male_df$anomalies <- legal_male_df$tot_legal - mean(legal_male_df$tot_legal)

# Create matrices
legal_male <- as.data.frame(legal_male_df[order(legal_male_df$index), ])
legal_male_mat <- fossil::create.matrix(legal_male,
                                        tax.name = "index",
                                        locality = "doy",
                                        time.col = NULL,
                                        time = NULL,
                                        abund = T,
                                        abund.col = "anomalies")

# PCA
legal_male_pca <- prcomp(legal_male_mat,
                         scale = F,
                         center = F)
summary(legal_male_pca)
plot(summary(legal_male_pca)$importance[2, ], 
     type = "h", 
     ylab = "Proportion of variance", 
     xlab = "PC")
barplot(-legal_male_pca$rotation[, 1],
        col = "red",
        ylim= range(-legal_male_pca$rotation[, 1]),
        main = "Legal Males",
        xlab = "doy",
        ylab = "Magnitude",
        las = 2)

# Try with vegan
legal_male_pca <- princomp(legal_male_mat, cor = F)
summary(legal_male_pca, loadings = T)
head(summary(legal_male_pca))

biplot(legal_male_mat)

# Outlier evaluation
legal_male_df <- filter(legal_male_df, depth < 400, tot_legal > 0)
sublegal_male_df <- filter(sublegal_male_df, depth < 400, sublegal > 0)
female_df <- filter(female_df, depth < 400, female > 0)


# PCA
variables <- c("year", "doy", "latitude", "depth")

legal_male_pca <- rda(legal_male_df[, variables], scale = T)
biplot(legal_male_pca, display = "species", scaling = "species")
head(summary(legal_male_pca))
legal_male_loadings <- scores(legal_male_pca, display = "species", scaling = 0)
sort(abs(legal_male_loadings[, 1]), decreasing = T)

sublegal_male_pca <- rda(sublegal_male_df[, variables], scale = T)
biplot(sublegal_male_pca)
head(summary(sublegal_male_pca))
sublegal_male_loadings <- scores(sublegal_male_pca, display = "species", scaling = 0)
sort(abs(sublegal_male_loadings[, 1]), decreasing = T)

female_pca <- rda(female_df[, variables], scale = T)
biplot(female_pca)
head(summary(female_pca))
female_loadings <- scores(female_pca, display = "species", scaling = 0)
sort(abs(female_loadings[, 1]), decreasing = T)
