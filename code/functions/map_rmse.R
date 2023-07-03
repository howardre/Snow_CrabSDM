map_rmse <- function(rmse_avg, area){
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
    geom_polygon(aes(long, lat, group = group), data = bering_sea,
                 fill = "lightyellow4", 
                 colour = "black",
                 inherit.aes = TRUE) +
    scale_size_area() +
    coord_sf(xlim = c(-180, -156), 
             ylim = c(54, 66), 
             expand = FALSE) +
    scale_x_continuous(breaks = EBS$lon.breaks) + 
    scale_y_continuous(breaks = EBS$lat.breaks) +
    scale_fill_viridis(option = "inferno") +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 14),
          axis.title = element_text(family = "serif", size = 18),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          strip.text = element_text(family = "serif", size = 18),
          legend.title = element_text(family = "serif", size = 16),
          legend.text = element_text(family = "serif", size = 14)) +
    labs(x = "Longitude",
         y = "Latitude") 
}
