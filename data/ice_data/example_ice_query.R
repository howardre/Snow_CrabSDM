# Code from Emily    
# produce ice cover time series for Bering (55-64N, 180-160W)
    
    library(tidyverse)
    library(ncdf4)
    library(zoo)
    library(maps)
    library(mapdata)
    library(chron)
    library(fields)
    library(oce)
    
    # set palettes
    new.col <- oceColorsPalette(64)
    cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    # set theme
    theme_set(theme_bw())
    
    
    ## load and process ------------------------
    
    # note that there are separate time series for 1950-1978 and 1979-present
    
    nc1 <- nc_open(here("data/ice_data", "ERA5_ice_1950-1978.nc"))
    
    # process
    
    ncvar_get(nc1, "time")   # hours since 1-1-1900
    raw <- ncvar_get(nc1, "time")
    h <- raw/24
    d1 <- dates(h, origin = c(1,1,1900))
    m1 <- months(d1)
    yr1 <- years(d1)
    
    x1 <- ncvar_get(nc1, "longitude")
    y1 <- ncvar_get(nc1, "latitude")
    
    ice1 <- ncvar_get(nc1, "siconc", verbose = F)
    dim(ice1) # 87 long, 37 lat, 203 months
    
    # reverse lat for plotting
    ice1 <- ice1[,37:1,]
    
    # reverse y too
    y1 <- rev(y1)
    
    ice1 <- aperm(ice1, 3:1)
    
    ice1 <- matrix(ice1, nrow=dim(ice1)[1], ncol=prod(dim(ice1)[2:3]))
    
    # plot to check
    
    ice.mean <- colMeans(ice1, na.rm=T)
    z <- t(matrix(ice.mean,length(y1)))
    image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x1, y1, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  
    
    # now the second time series
    
    nc2 <- nc_open(here("data/ice_data", "ERA5_ice_1979-2022.nc"))
    
    # process
    
    ncvar_get(nc2, "time")   # hours since 1-1-1900
    raw <- ncvar_get(nc2, "time")
    h <- raw/24
    d2 <- dates(h, origin = c(1,1,1900))
    m2 <- months(d2)
    yr2 <- years(d2)
    
    x2 <- ncvar_get(nc2, "longitude")
    y2 <- ncvar_get(nc2, "latitude")
    
    expver <-  ncvar_get(nc2, "expver", verbose = F)
    expver # 1 and 5??
    
    lat <- rep(y2, length(x2))
    lon <- rep(x2, each = length(y2))
   
    ice2 <- ncvar_get(nc2, "siconc", verbose = F)
    dim(ice2) # 87 long, 37 lat, 2 expver, 203 months
    
    # expver1 - this is ERA5
    
    ice2 <- ice2[,,1,]
    
    # reverse lat for plotting
 #   ice2 <- ice2[,37:1,]
    
    # reverse y too
#    y2 <- rev(y2)
    
    ice2 <- aperm(ice2, 3:1)
    
    ice2 <- matrix(ice2, nrow = dim(ice2)[1], ncol = prod(dim(ice2)[2:3]))
    

    data.frame(lon = lon, lat = lat, ice2) %>%
      pivot_longer(cols = c(3:304), names_to = "month", values_to = "ice") %>%
      mutate(month = rep(m2, ), year = rep(yr2, ), ice = ice2) -> ice_latlon2
    #dplyr::filter(lon > -166, lat < 63) 
    
    write.csv(ice_latlon2, "./Data/ice_latlon2.csv")
    
    
    # plot to check
    
    ice.mean <- colMeans(ice2, na.rm=T)
    z <- t(matrix(ice.mean,length(y2)))
    image.plot(x2,y2,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x2, y2, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  
    
    
    # check dimensions
    identical(x1, x2)
    identical(y1, y2)
    

    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y1, length(x1))
    lon <- rep(x1, each = length(y1))
    
  
    ice <- rbind(ice1, ice2)
    
    # drop E of 165 and N of 63
    drop <- lon > -165 | lat > 63
    ice[,drop] <- NA
    
    # plot to check
    ice.mean <- colMeans(ice, na.rm=T)
    z <- t(matrix(ice.mean,length(y1)))
    image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
    contour(x1, y1, z, add=T)
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3") # perfecto
    
    dimnames(ice) <- list(as.character(c(d1, d2)), paste("N", lat, "E", lon, sep=""))
    
    f <- function(x) colMeans(x, na.rm = T)
    
    m <- c(as.character(m1), as.character(m2))
    
    yr <- c(as.numeric(as.character(yr1)), as.numeric(as.character(yr2)))
    
    means <- data.frame(month = m,
                        year = as.numeric(as.character(yr)),
                        ice = rowMeans(ice, na.rm = T)) 
    
    
    ggplot(means, aes(year, ice, color = month)) +
      geom_line()
    
    # drop Oct - Dec
    means <- means %>%
      filter(!month %in% c("Oct", "Nov", "Dec"))
    
    
    # pivot wider
    means <- means %>% 
      pivot_wider(values_from = ice, names_from = month) %>%
      filter(year %in% 1953:2021) # fixed values through 1952, 2022 incomplete!
    
    means[,2:5] <- apply(means[,2:5], 2, scale)
    
    plot <- means %>%
      pivot_longer(cols = -year)
    
    ggplot(plot, aes(year, value, color = name)) +
      geom_line()
    
    
    # generate Jan-Feb and Mar-Apr means
    means$JanFeb_ice <- apply(means[,2:3], 1, mean)
    means$MarApr_ice <- apply(means[,4:5], 1, mean)
    
    # clean up
    means <- means %>%
      select(year, JanFeb_ice, MarApr_ice)
    
    plot <- means %>%
      pivot_longer(cols = -year)
    
    ggplot(plot, aes(year, value, color = name)) +
      geom_line()
      
      
# save 
    write.csv(means, "./data/ice.csv", row.names = F)
    

# Third timeseries
    # Load third dataset
    
    nc3 <- nc_open("./Data/ERA5_ice_sst_1979-2021.nc")
    
    #Process
    
    h <- (ncvar_get(nc3, "time")/24)
    d3 <- dates(h, origin = c(1, 1, 1900))  
    m3 <- months(d3)
    yr3 <- years(d3)
    
    x3 <- ncvar_get(nc3, "longitude")
    y3 <- ncvar_get(nc3, "latitude")
    
    #Ice cover
    ice3 <- ncvar_get(nc3, "siconc", verbose = F)
    dim(ice3) #101 lon, 45 lat, 516 months
    
    #Reverse lat for plotting
    ice3 <- ice3[,45:1,]
    
    #Reverse y too
    y3 <- rev(y3)
    
    #Create ice data matrix??
    ice3 <- aperm(ice3, 3:1)
    
    ice3 <- matrix(ice3, nrow = dim(ice3)[1], ncol = prod(dim(ice3)[2:3]))
    
    ice.mean <- colMeans(ice3, na.rm = T)
    
    z3 <- t(matrix(ice.mean, length(y3))) #mean across all years
    
    #Plot
    image.plot(x3, y3, z3, col=oceColorsPalette(64), xlab = "", ylab = "")
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")
    contour(x3, y3, z, add=T)
    
    #SST
    sst <- ncvar_get(nc3, "sst", verbose = F)-273.15 #change to celcius?
    
    #Reverse lat for plotting
    sst <- sst[,45:1,]
    
    #Create sst data matrix??
    sst <- aperm(sst, 3:1)
    
    sst <- matrix(sst, nrow = dim(sst)[1], ncol = prod(dim(sst)[2:3]))
    
    sst.mean <- colMeans(sst, na.rm = T)
    
    z <- t(matrix(sst.mean, length(y3)))
    
    #Plot
    image.plot(x3, y3, z, col=oceColorsPalette(64), xlab = "", ylab = "")
    map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")
    contour(x3, y3, z, add=T)
    
    #Calculate monthly means for ice and SST
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y3, length(x3))
    lon <- rep(x3, each = length(y3))
    
    
    # drop E of 165 and N of 63
    drop <- lon > -165 | lat > 63
    ice3[,drop] <- NA
    sst[,drop] <- NA
    
    #Calculate means
    ice3.mean <- colMeans(ice3, na.rm=T)
    sst.mean <- colMeans(sst, na.rm=T)
    
    #Set dimnames
    dimnames(ice3)<-list(as.character(d3), paste("N", lat, "E", lon, sep=""))
    
    m <- as.character(m3)
    yr <- as.numeric(as.character(yr3))
    
    #Means
    data.frame(month = m,
               year = yr,
               ice = rowMeans(ice3, na.rm=T),
               sst = rowMeans(sst, na.rm=T)) -> nc3varmeans
    
    #plot
    ice_plot <- ggplot(nc3varmeans, aes(year, ice, color = month)) +
      geom_line()
    
    sst_plot <- ggplot(nc3varmeans, aes(year, sst, color = month)) +
      geom_line()

    
#creating a datasets with lat/lon, month, and year for ice and sst----------------------------------------
  #ICE------------------
  # Load third dataset
  nc3 <- nc_open(here("data/ice_data", "ERA5_ice_1979-2022.nc"))
    
  #Process
  h <- (ncvar_get(nc3, "time")/24)
  d3 <- dates(h, origin = c(1, 1, 1900))  
  m3 <- months(d3)
  yr3 <- chron::years(d3)
    
  x3 <- ncvar_get(nc3, "longitude")
  y3 <- ncvar_get(nc3, "latitude")
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y3, length(x3))
  lon <- rep(x3, each = length(y3))
    
  #Ice cover
  ice3 <- ncvar_get(nc3, "siconc", verbose = F)
  dim(ice3) #101 lon, 45 lat, 516 months
    
  #Reverse lat for plotting
  #ice3 <- ice3[,45:1,]
    
  #Reverse y too
  #y3 <- rev(y3)
    
  # Create ice data matrix??
  ice3 <- aperm(ice3, c(4:1))
  mat_ice3 <- t(matrix(ice3, nrow = dim(ice3)[1], ncol = prod(dim(ice3)[c(2:3)])))

  data.frame(lon = lon3, lat = lat3,  mat_ice3) %>%
    pivot_longer(cols = c(3:518), names_to = "month", values_to = "ice") %>%
    mutate(month = rep(m3, 4545), year = rep(yr3, 4545), ice = ice) -> ice_latlon2
    #dplyr::filter(lon > -166, lat < 63) 

  write.csv(ice_latlon2, "./Data/ice_latlon2.csv")
  
  #SST-----------------------
  #SST
  sst <- ncvar_get(nc3, "sst", verbose = F)-273.15 #change to celcius?
  
  #Reverse lat for plotting
  sst <- sst[,45:1,]
  
  #Create sst data matrix??
  sst <- aperm(sst, 3:1)
  
  mat_sst <- t(matrix(sst, nrow = dim(sst)[1], ncol = prod(dim(sst)[2:3])))
  
  data.frame(lat = lat3, lon = lon3, mat_sst) %>%
    pivot_longer(cols = c(3:518), names_to = "month", values_to = "sst") %>%
    mutate(month = rep(m3, 4545), year = rep(yr3, 4545), sst = sst) -> sst_latlon
    #dplyr::filter(lon > -166, lat < 63) -> sst_latlon
  
  write.csv(sst_latlon, "./Data/sst_latlon.csv")
  
  
#############################################################################################################################
#############################################################################################################################
# NEW DATA
  
  # 1975-2021 data
    # ICE ------------------------------------------------------------------------
    ice_nc4 <- nc_open(here("data/ice_data", "ERA5_ice_1979-2022.nc"))
    
    ice4 <- ncvar_get(ice_nc4, "siconc", verbose = F)
    ice4 <- ifelse(ice4<0.0005, 0, ice4)
    dim(ice4) #101 lon, 33 lat, 564 months
  
    #Process
    h <- (ncvar_get(ice_nc4, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(ice_nc4, "longitude")
    y4 <- ncvar_get(ice_nc4, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
    
    #Reverse lat for plotting
    #ice4 <- ice4[,33:1,]
    
    #Reverse y too
    #y4 <- rev(y4)
    
    #Create ice data matrix??
    ice4 <- aperm(ice4, 4:1)
    mat_ice4 <- t(matrix(ice4, nrow = dim(ice4)[1], ncol = prod(dim(ice4)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_ice4) %>%
      pivot_longer(cols = c(3:566), names_to = "month", values_to = "ice") %>%
      mutate(month = rep(m4, 3333), year = rep(yr4, 3333), ice = ice) %>%
      na.omit()-> ice_latlon
    
    #2022 data
    nc4_22 <- nc_open(here("data/ice_data", "ERA5_ice_1979-2022.nc"))
    
    ice22 <- ncvar_get(nc4_22, "siconc", verbose = F)
    
    ice22 <- ice22[,,1,]
    
    ice22 <- ifelse(ice22<0.0005, 0, ice22)
    dim(ice22) #101 lon, 33 lat, 2 months
    
    #Process
    h <- (ncvar_get(nc4_22, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(nc4_22, "longitude")
    y4 <- ncvar_get(nc4_22, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
    
    #Reverse lat for plotting
    #ice4 <- ice4[,33:1,]
    
    #Reverse y too
    #y4 <- rev(y4)
    
    #Create ice data matrix??
    ice22 <- aperm(ice22, 3:1)
    mat_ice22 <- t(matrix(ice22, nrow = dim(ice22)[1], ncol = prod(dim(ice22)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_ice22) %>%
      pivot_longer(cols = c(3:306), names_to = "month", values_to = "ice") %>%
      mutate(month = rep(m4, 2997), year = rep(yr4, 2997), ice = ice) %>%
      na.omit()-> ice_latlon22
    
    #bind with earlier ice data
    rbind(ice_latlon, ice_latlon22) -> ice_latlon
    
    #write csv
    write.csv(ice_latlon, "./Data/ice_latlon.csv")
    
    # SST --------------------------------------------------
   
    sst_nc4 <- nc_open("./Data/ERA5_sst2_1975-2021.nc")
    
    sst4 <- ncvar_get(sst_nc4, "sst", verbose = F)-273.15
    dim(sst4) #101 lon, 33 lat, 564 months
    
      #Process
      h <- (ncvar_get(sst_nc4, "time")/24)
      d4 <- dates(h, origin = c(1, 1, 1900))  
      m4 <- months(d4)
      yr4 <- chron::years(d4)
      
      x4 <- ncvar_get(sst_nc4, "longitude")
      y4 <- ncvar_get(sst_nc4, "latitude")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y4, length(x4))
      lon <- rep(x4, each = length(y4))
      
      #Reverse lat for plotting
      #sst4 <- sst4[,33:1,]
      
      #Reverse y too
      #y4 <- rev(y4)
      
      #Create ice data matrix??
      sst4 <- aperm(sst4, 3:1)
      mat_sst4 <- t(matrix(sst4, nrow = dim(sst4)[1], ncol = prod(dim(sst4)[2:3])))
      
      data.frame(lon = lon, lat = lat,  mat_sst4) %>%
        pivot_longer(cols = c(3:566), names_to = "month", values_to = "sst") %>%
        mutate(month = rep(m4, 3333), year = rep(yr4, 3333), sst = sst) %>%
        na.omit()-> sst_latlon
      
      #2022 data
      sst22 <- ncvar_get(nc4_22, "sst", verbose = F) -273.15
  
      sst22 <- sst22[,,1,]
      
      dim(sst22) #101 lon, 33 lat, 2 months
      
      #Process
      h <- (ncvar_get(nc4_22, "time")/24)
      d4 <- dates(h, origin = c(1, 1, 1900))  
      m4 <- months(d4)
      yr4 <- chron::years(d4)
      
      x4 <- ncvar_get(nc4_22, "longitude")
      y4 <- ncvar_get(nc4_22, "latitude")
      
      # Keep track of corresponding latitudes and longitudes of each column:
      lat <- rep(y4, length(x4))
      lon <- rep(x4, each = length(y4))
      
      #Reverse lat for plotting
      #ice4 <- ice4[,33:1,]
      
      #Reverse y too
      #y4 <- rev(y4)
      
      #Create ice data matrix??
      sst22 <- aperm(sst22, 3:1)
      mat_sst22 <- t(matrix(sst22, nrow = dim(sst22)[1], ncol = prod(dim(sst22)[2:3])))
      
      data.frame(lon = lon, lat = lat,  mat_sst22) %>%
        pivot_longer(cols = c(3:12), names_to = "month", values_to = "sst") %>%
        mutate(month = rep(m4, 3333), year = rep(yr4, 3333), sst = sst) %>%
        na.omit()-> sst_latlon22
      
      #bind with earlier ice data
      rbind(sst_latlon, sst_latlon22) -> sst_latlon
      
      write.csv(sst_latlon, "./Data/sst_latlon.csv")
    
    
    
    
   