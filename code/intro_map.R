# Libraries
library(marmap)
library(ggOceanMaps)
library(ggplot2)

# ggOceanMaps method
bathy_palette <- colorRampPalette(c("lightsteelblue1", "lightsteelblue4"))(9)
BS_bathy <- data.frame(lon = c(-176.6, -156.5),
                       lat = c(51, 65.5))
bering_map <- basemap(data = BS_bathy,
                      bathymetry = TRUE,
                      rotate = TRUE,
                      legends = FALSE,
                      land.col = "wheat4",
                      grid.col = NA,
                      lon.interval = 5,
                      lat.interval = 2) +
  geom_polygon(data = transform_coord(BS_bathy),
               aes(x = lon, y = lat),
               fill = NA) +
  scale_fill_manual(values = bathy_palette) +
  labs(x = "Longitude",
       y = "Latitude") 

# marmap method
BS_bathy <- getNOAA.bathy(lon1 = -176.5, lon2 = -156.5, 
                          lat1 = 51, lat2 = 65.5, 
                          resolution = 1)
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

windows(width = 17,
        height = 14,)
par(mfrow = c(1, 1),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
plot.bathy(BS_bathy,
           image = T,
           axes = T,
           lwd = 0.03,
           land = T,
           n = 0, 
           bpal = list(c(0,
                         max(BS_bathy),
                         greys),
                       c(min(BS_bathy),
                         0,
                         blues)),
           xlim = c(-176.5, -156.8),
           ylim = c(51, 65.5),
           ylab = "Latitude °N",
           xlab = "Longitude °W",
           main = "",
           cex.lab = 2,
           cex.main = 2.5,
           cex.axis = 1.8)