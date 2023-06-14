grid_development <- function(train_data){
  nlat = 40
  nlon = 60
  latd = seq(min(train_data$latitude), max(train_data$latitude), length.out = nlat)
  lond = seq(min(train_data$longitude), max(train_data$longitude), length.out = nlon)
  spatial_grid <- expand.grid(lond, latd) # create grid
  names(spatial_grid) <- c('longitude', 'latitude')
  spatial_grid$dist <- NA # calculate distance from nearest station
  for (k in 1:nrow(spatial_grid)) {
    dist <-  distance_function(spatial_grid$latitude[k],
                               spatial_grid$longitude[k],
                               train_data$latitude,
                               train_data$longitude)
    spatial_grid$dist[k] <- min(dist)
  }
  
  # Add in LOESS interpolated values for phi and depth
  spatial_grid$depth <- as.vector(predict(depth_loess,
                                          newdata = spatial_grid))
  spatial_grid$phi <- as.vector(predict(phi_loess,
                                        newdata = spatial_grid))
  
  # Add median values for all other variables
  spatial_grid$year_f <- as.factor('2010')
  spatial_grid$julian <- median(train_data$julian, na.rm = T)
  spatial_grid$temperature <- median(train_data$temperature, na.rm = T)
  spatial_grid$ice_mean <- median(train_data$ice_mean, na.rm = T)
  spatial_grid$log_pcod_cpue <- median(train_data$log_pcod_cpue, na.rm = T)
  return(spatial_grid)
}