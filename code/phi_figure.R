library(raster)
library(here)
library(ggplot2)
library(viridis)
library(tidyverse)
library(terra)
library(maps)
library(mapdata)
library(fields)

# Load phi data, specify CRS, convert to lat lon
phi_raster <- raster(here('data', 'EBS_phi_1km.gri'))
phi_cropped <- crop(phi_raster, c(-1400000, -200000, 500000, 1850000))
phi_pts <- rasterToPoints(phi_raster, spatial = T)
proj4string(phi_pts) # warning here related to the retirement of rgdal package
phi_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
phi_trans <- spTransform(phi_pts, CRS(phi_prj))
proj4string(phi_trans)
phi_data <- as.data.frame(phi_trans)
phi_matrix <- as.matrix(phi_data)

# Convert to plottable data with terra
phi_terra <- ext(apply(phi_matrix[, 2:3], 2, range))
phi_rast <- rast(phi_terra, ncol = 300, nrow = 300)
phi_x <- terra::rasterize(phi_matrix[, 2:3], phi_rast, phi_matrix[, 1], fun = mean)

# Plot
windows(width = 10, height = 10)
par(mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(phi_x,
      xlim = c(-179.5, -156),
      ylim = c(54, 66),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
image(phi_x,
      col = hcl.colors(18, "YlOrRd", rev = TRUE),
      xlim = c(-179.5, -156),
      ylim = c(54, 66),     
      ylab = "Latitude \u00B0N",
      xlab = "Longitude \u00B0W",
      cex.main = 2.3,
      cex.lab = 2.5,
      cex.axis = 2.3,
      las = 1)
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = hcl.colors(18, "YlOrRd", rev = TRUE),
           legend.shrink = 0.2,
           smallplot = c(.24, .26, .17, .35),
           legend.cex = 1.8,
           axis.args = list(cex.axis = 2.1,
                            family = "serif"),
           legend.width = 0.8,
           legend.mar = 6,
           zlim = c(min(phi_data$X3.pred, na.rm = T), 
                    max(phi_data$X3.pred, na.rm = T)),
           legend.args = list(expression(paste("Grain Size (", phi, ")")),
                              side = 2,
                              cex = 2.3,
                              line = 1.3,
                              family =  "serif"))
dev.copy(jpeg,
         here('results',
              'grain_size.jpg'),
         height = 10,
         width = 10,
         res = 300,
         units = 'in')
dev.off()
