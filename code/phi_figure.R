library(raster)
library(here)
library(ggplot2)
library(viridis)
library(tidyverse)
library(terra)

# Load phi data
phi_raster <- raster(here('data', 'EBS_phi_1km.gri'))
phi_cropped <- crop(phi_raster, c(-1400000, -200000, 500000, 1850000))
phi_pts <- rasterToPoints(phi_raster, spatial = T)
proj4string(phi_pts) # warning here related to the retirement of rgdal package
phi_prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
phi_trans <- spTransform(phi_pts, CRS(phi_prj))
gridded(phi_trans) = TRUE
phi_back <- raster(phi_trans)
projection(phi_back) <- CRS(phi_prj)

proj4string(phi_trans)
phi_data <- as.data.frame(phi_trans)

image(phi_back)

ggplot() +
  geom_tile(data = phi_data, 
            aes(x = coords.x1, y = coords.x2, fill = X3.pred)) +
  scale_fill_viridis()

ggplot() +
  geom_raster(data = phi_pts, 
              aes(x = coords.x1, y = coords.x2, fill = X3.pred))
