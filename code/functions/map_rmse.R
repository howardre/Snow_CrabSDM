map_rmse <- function(rmse_avg){
  # Get average RMSE per station
  df <- merge(rmse_avg,
              EBS_trans,
              by.x = "station",
              by.y =  "STATIONID") # make spatial
  sf <- st_as_sf(df, # turn into sf object to plot
                 crs = 4269)
  # Plot RMSE
  ggplot() +
    geom_sf(data = sf,
            aes(fill = rmse)) +
    scale_x_continuous(name = "Longitude", 
                       breaks = EBS$lon.breaks) + 
    scale_y_continuous(name = "Latitude", 
                       breaks = EBS$lat.breaks)
}