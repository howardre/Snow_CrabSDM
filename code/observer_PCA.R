# Title: Observer Data PCA
# Purpose: Determine if PCA can be used for observer data
# Author:  Rebecca Howard

# Load libraries ----
library(here)
library(vegan)
library(dplyr)
library(corrplot)
library(ggplot2)

# Load data ----
observer_all <- readRDS(here('data/Snow_CrabData', 'observer_all.rds'))
crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_summary.rds'))

# Convert to lat lon, group
observer_locs <- observer_all %>% tidyr::extract(geometry, c('lon', 'lat'), '\\((.*), (.*)\\)', convert = TRUE)

observer_df <- observer_locs %>%
  group_by(year_lag, STATIONID) %>%
  summarise(obs_male_legal = mean(tot_legal),
            obs_male_sub = mean(sublegal),
            obs_female = mean(female),
            latitude = mean(lat),
            longitude = mean(lon))

# Filter to groups
observer_df$lncpue_male_legal <- log(observer_df$obs_male_legal + 1)
observer_df$lncpue_male_sublegal <- log(observer_df$obs_male_sub + 1)
observer_df$lncpue_female <- log(observer_df$obs_female + 1)

# Histograms
hist(observer_df$lncpue_male_legal)
hist(observer_df$lncpue_male_sublegal)
hist(observer_df$lncpue_female)

# Make matrices ----
# Year and station ID for each group
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

# Flip the matrix to get loadings for the stations instead
# Year and station ID for each group
legal_male_mat_st <- fossil::create.matrix(as.data.frame(observer_df),
                                           tax.name = "year_lag",
                                           locality = "STATIONID",
                                           time.col = NULL,
                                           time = NULL,
                                           abund = T,
                                           abund.col = "lncpue_male_legal")
legal_male_stations <- legal_male_mat_st[, colSums(abs(legal_male_mat_st)) != 0]

sublegal_male_mat_st <- fossil::create.matrix(as.data.frame(observer_df),
                                              tax.name = "year_lag",
                                              locality = "STATIONID",
                                              time.col = NULL,
                                              time = NULL,
                                              abund = T,
                                              abund.col = "lncpue_male_sublegal")
sublegal_male_stations <- sublegal_male_mat_st[, colSums(abs(sublegal_male_mat_st)) != 0]

female_mat_st <- fossil::create.matrix(as.data.frame(observer_df),
                                       tax.name = "year_lag",
                                       locality = "STATIONID",
                                       time.col = NULL,
                                       time = NULL,
                                       abund = T,
                                       abund.col = "lncpue_female")
female_stations <- female_mat_st[, colSums(abs(female_mat_st)) != 0]

# Correlation plots
corrplot(cor(legal_male_mat), is.corr = FALSE)
corrplot(cor(sublegal_male_mat), is.corr = FALSE)
corrplot(cor(female_mat), is.corr = FALSE)

# Scale values
legal_male_scaled <- scale(legal_male_mat)
sublegal_male_scaled <- scale(sublegal_male_mat)
female_scaled <- scale(female_mat)

legal_male_scaled_st <- scale(legal_male_stations)
sublegal_male_scaled_st <- scale(sublegal_male_stations)
female_scaled_st <- scale(female_stations)

# PCA ----
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
# For years
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

# Components
plot(legal_male_pca, type = "l")
plot(sublegal_male_pca, type = "l")
plot(female_pca, type = "l")

# For stations
legal_male_pca_st <- prcomp(legal_male_scaled_st,
                            scale = F,
                            center = F)
summary(legal_male_pca_st)
head(summary(legal_male_pca_st))

legal_male_pca$sdev # eigenvalues
legal_male_pca$x[, 1][1:10] # first 10 values of pca 1
legal_male_pca$rotation[, 1] # loadings for first pca

legal_male_pca_st$sdev
legal_male_pca_st$loadings[, 1] # Use these in the model
legal_male_pca_st$scores[, 1] 
legal_male_pca_st$scale
legal_male_pca_st$center

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

# Match to data ----
legal_male_loadings <- as.data.frame(legal_male_pca$loadings[, 1])
legal_male_loadings <- tibble::rownames_to_column(legal_male_loadings, "year")
colnames(legal_male_loadings)[2] <- "legal_male_loadings"

legal_male_scores <- as.data.frame(legal_male_pca$scores[, 1])
legal_male_scores <- tibble::rownames_to_column(legal_male_scores, "station")
colnames(legal_male_scores)[2] <- "legal_male_scores"

plot(legal_male_loadings,
     main = "Legal Male Loadings",
     ylab = "value")

sublegal_male_loadings <- as.data.frame(sublegal_male_pca$loadings[, 1])
sublegal_male_loadings <- tibble::rownames_to_column(sublegal_male_loadings, "year")
colnames(sublegal_male_loadings)[2] <- "sublegal_male_loadings"

sublegal_male_scores <- as.data.frame(sublegal_male_pca$scores[, 1])
sublegal_male_scores <- tibble::rownames_to_column(sublegal_male_scores, "station")
colnames(sublegal_male_scores)[2] <- "sublegal_male_scores"

plot(sublegal_male_loadings,
     main = "Sublegal Male Loadings",
     ylab = "value")

female_loadings <- as.data.frame(female_pca$loadings[, 1])
female_loadings <- tibble::rownames_to_column(female_loadings, "year")
colnames(female_loadings)[2] <- "female_loadings"

female_scores <- as.data.frame(female_pca$scores[, 1])
female_scores <- tibble::rownames_to_column(female_scores, "station")
colnames(female_scores)[2] <- "female_scores"

plot(female_loadings,
     main = "Female Loadings",
     ylab = "value")

crab_summary$legal_male_loading <- legal_male_loadings$legal_male_loadings[match(crab_summary$year, legal_male_loadings$year)]
crab_summary$sublegal_male_loading <- sublegal_male_loadings$sublegal_male_loadings[match(crab_summary$year, sublegal_male_loadings$year)]
crab_summary$female_loading <- female_loadings$female_loadings[match(crab_summary$year, female_loadings$year)]



# Save data
saveRDS(crab_summary, file = here('data/Snow_CrabData', 'crab_pca.rds'))

# Plot ----
legal_male_scores$latitude <- observer_df$latitude[match(legal_male_scores$station, observer_df$STATIONID)]
legal_male_scores$longitude <- observer_df$longitude[match(legal_male_scores$station, observer_df$STATIONID)]
sublegal_male_scores$latitude <- observer_df$latitude[match(sublegal_male_scores$station, observer_df$STATIONID)]
sublegal_male_scores$longitude <- observer_df$longitude[match(sublegal_male_scores$station, observer_df$STATIONID)]
female_scores$latitude <- observer_df$latitude[match(female_scores$station, observer_df$STATIONID)]
female_scores$longitude <- observer_df$longitude[match(female_scores$station, observer_df$STATIONID)]

legal_male_sf <- st_as_sf(legal_male_scores,
                          coords = c("longitude", "latitude"), 
                          crs = 4269)
sublegal_male_sf <- st_as_sf(sublegal_male_scores,
                             coords = c("longitude", "latitude"), 
                             crs = 4269)
female_sf <- st_as_sf(female_scores,
                      coords = c("longitude", "latitude"), 
                      crs = 4269)

ggplot() +
  geom_sf(data = legal_male_sf,
          aes(color = legal_male_scores),
          size = 4) 

ggplot() +
  geom_sf(data = sublegal_male_sf,
          aes(color = sublegal_male_scores),
          size = 4) 

ggplot() +
  geom_sf(data = female_sf,
          aes(color = female_scores),
          size = 4) 

# NMDS ----
# Ultimately inappropriate for replacing PCA loadings due to non-orthogonal 
# New matrices
# Year and station ID for each group
legal_male_mat2 <- fossil::create.matrix(as.data.frame(observer_df),
                                         tax.name = "STATIONID",
                                         locality = "year_lag",
                                         time.col = NULL,
                                         time = NULL,
                                         abund = T,
                                         abund.col = "obs_male_legal")

sublegal_male_mat2 <- fossil::create.matrix(as.data.frame(observer_df),
                                            tax.name = "STATIONID",
                                            locality = "year_lag",
                                            time.col = NULL,
                                            time = NULL,
                                            abund = T,
                                            abund.col = "obs_male_sub")

female_mat2 <- fossil::create.matrix(as.data.frame(observer_df),
                                     tax.name = "STATIONID",
                                     locality = "year_lag",
                                     time.col = NULL,
                                     time = NULL,
                                     abund = T,
                                     abund.col = "obs_female")

# Distance matrices
legal_male_dist <- as.matrix(vegdist(legal_male_mat2,
                                     method = "bray"), labels = TRUE)

set.seed(1993)
legal_male_nmds <- metaMDS(legal_male_dist, trace = FALSE)

plot(legal_male_nmds)
stressplot(legal_male_nmds)

legal_male_vec <- 1:10
legal_male_stress <- numeric(length(legal_male_vec))
legal_male_distance <- metaMDSdist(legal_male_mat2, trace = FALSE)
set.seed(2)
for(i in seq_along(legal_male_vec)) {
  legal_male_nmds <- metaMDSiter(legal_male_distance, 
                                 k = i,
                                 trace = FALSE)
  legal_male_stress[i] <- legal_male_nmds$stress
}

plot(k_vec, 
     legal_male_stress, 
     type = "b",
     ylab = "Stress",
     xlab = "Dimensions")

legal_male_scores <- scores(legal_male_nmds)

