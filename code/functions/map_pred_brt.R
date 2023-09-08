map_pred_brt <- function(spatial_grid, train_data, title){
  # Create palette dependent on scale
  my_color = colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  color_levels = 100
  max_absolute_value = max(abs(c(min(spatial_grid$pred_brt, na.rm = T),
                                 max(spatial_grid$pred_brt, na.rm = T))))
  color_sequence = seq(max(spatial_grid$pred_brt, na.rm = T), 
                       min(spatial_grid$pred_brt, na.rm = T),
                       length.out = color_levels + 1)
  n_in_class = hist(spatial_grid$pred_brt, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  
  # Use for matrix
  nlat = 40
  nlon = 60
  latd = seq(min(train_data$latitude), max(train_data$latitude), length.out = nlat)
  lond = seq(min(train_data$longitude), max(train_data$longitude), length.out = nlon)
  
  # Make map
  windows(width = 12, height = 10)
  par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0),
      family = "serif")
  image(lond,
        latd,
        t(matrix(spatial_grid$pred_brt,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-181, -156),
        ylim = range(train_data$latitude, na.rm = TRUE) + c(-.4, .5),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(spatial_grid$pred_brt,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include],
        ylab = "Latitude \u00B0N",
        xlab = "Longitude \u00B0W",
        xlim = c(-181, -156),
        ylim = range(train_data$latitude, na.rm = TRUE) + c(-.4, .5),
        main = title,
        cex.main = 2.3,
        cex.lab = 2.5,
        cex.axis = 2.3,
        las = 1)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.2, .24, .17, .38),
             legend.cex = 1.8,
             axis.args = list(cex.axis = 2.1,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(spatial_grid$pred_brt, na.rm = T), 
                      max(spatial_grid$pred_brt, na.rm = T)),
             legend.args = list("log(count+1)",
                                side = 2,
                                cex = 2.3,
                                line = 1.3,
                                family =  "serif"))
}
