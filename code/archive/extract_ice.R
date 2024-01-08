# Modified code from Emily

# Libraries
library(ncdf4)
library(here)
library(tidyverse)
library(chron)
library(fields)
library(oce)

# Ice data
nc <- nc_open(here("data/ice_data", "ERA5_ice_1979-2022.nc"))

ice_data <- ncvar_get(nc, "siconc", verbose = F)

ice_data <- ice_data[, , 1 ,]

ice_data <- ifelse(ice_data < 0.0005, 0, ice_data)
dim(ice_data) # 81 lon, 37 lat, 304 months

#Process
h <- (ncvar_get(nc, "time") / 24)
d <- dates(h, origin = c(1, 1, 1900))  
m <- months(d)
yr <- chron::years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))

#Create ice data matrix
ice_data <- aperm(ice_data, 3:1)
mat_ice <- t(matrix(ice_data, 
                    nrow = dim(ice_data)[1], 
                    ncol = prod(dim(ice_data)[2:3])))

data.frame(lon = lon, lat = lat,  mat_ice) %>%
  pivot_longer(cols = c(3:306), names_to = "month", values_to = "ice") %>%
  mutate(month = rep(m, 2997), year = rep(yr, 2997), ice = ice) %>%
  na.omit()-> ice_latlon

write.csv(ice_latlon, "./data/ice_latlon.csv")
