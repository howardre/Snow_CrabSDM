library(here)
library(ncdf4)
library(tidync)
require(tidyverse)

# Currently unable to download directly from THREDDS
# create objects for known lats and longs and xi and eta axes
nc_grid <- nc_open(here('data', '/Bering10K_extended_grid.nc'))

lats <- ncvar_get(nc_grid, "lat_rho")
lons <- ncvar_get(nc_grid, "lon_rho")

nc_close(nc_grid)

# read in files
hindcast_temp_file_list <- list.files(path = here("data", "ROMS/"))

prestring <- here("data", "ROMS", "/")

hindcast_temp_dat_list <- list()
for(i in hindcast_temp_file_list){
  hindcast_temp_dat_list[[i]] <- paste0(prestring, i)
  hindcast_temp_dat_list
}

hindcast_temp_df_list <- list()
for(i in hindcast_temp_dat_list){
  hindcast_temp_df_list[[i]] <- tidync(i) %>%
    hyper_tibble(select_var = "temp")
  hindcast_temp_df_list
}

hindcast_temp_dfs <- bind_rows(hindcast_temp_df_list)

# add in lat/longs matched to xi/eta 
hindcast_temp_dfs$lon <- lons[cbind(hindcast_temp_dfs$xi_rho, hindcast_temp_dfs$eta_rho)]
hindcast_temp_dfs$lat <- lats[cbind(hindcast_temp_dfs$xi_rho, hindcast_temp_dfs$eta_rho)]

# create object for time axis
hindcast_temp_dfs$DateTime <- as.POSIXct(hindcast_temp_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")

hindcast_temp_dfs$date <- as.Date(hindcast_temp_dfs$DateTime) # date in date format
hindcast_temp_dfs$month <- month(hindcast_temp_dfs$date)
hindcast_temp_dfs$year <- year(hindcast_temp_dfs$date)

# remove all months aside from Feb - Aug
months <- c(5:8)

hindcast_temp_dfs_trim <- hindcast_temp_dfs %>%
  filter(., month %in% months) 

hindcast_temp_dfs_trim <- select(hindcast_temp_dfs_trim, -DateTime)

hindcast_temps <- hindcast_temp_dfs_trim %>% 
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon)) %>%
  filter(between(lat, 54, 66) &
           between(lon, -178, -159)) 

saveRDS(hindcast_temps, file = here('data', 'hindcast_temps.rds'))