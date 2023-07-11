map_rmse <- function(rmse_avg, plot_title){
  # Get average RMSE per station
  df <- merge(rmse_avg,
              EBS_trans,
              by.x = "station",
              by.y =  "STATIONID") # make spatial
  sf <- sf::st_as_sf(df, # turn into sf object to plot
                     crs = st_crs(4269))
  # Plot RMSE
  windows(width = 12, height = 10)
  par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0),
      family = "serif")
  ggplot() +
    geom_polygon(aes(long, lat, group = group),
                 data = bering_sea,
                 fill = "lightyellow4",
                 colour = "black") +
    geom_sf(data = sf,
            aes(fill = rmse),
            inherit.aes = FALSE) +
    coord_sf(xlim = c(-179.5, -156),  # cannot extend out to -180, cuts off latitude labels
             ylim = c(54, 66), 
             expand = FALSE) +
    scale_size_area() +
    scale_fill_viridis(option = "inferno") +
    theme_classic() +
    theme(panel.background = element_rect(fill = "gray91", colour = "gray91"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 24, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 16),
          axis.title = element_text(family = "serif", size = 20),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          strip.text = element_text(family = "serif", size = 20),
          legend.title = element_text(family = "serif", size = 18),
          legend.text = element_text(family = "serif", size = 16)) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "RMSE",
         title = plot_title) 
}
