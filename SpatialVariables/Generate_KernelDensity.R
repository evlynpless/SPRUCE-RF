#Instruction to make kernel density map from sample site coordinates, with the dimensions below
#extent      : 21.5 42.5 -5 22.5  (xmin, xmax, ymin, ymax) #MAPS1: 29.5 44.5 -7.5 13.0
#dimensions  : 3300, 2520, 8316000  (nrow, ncol, ncell) #MAPS1: 1800, 2460, 4428000
#projection : stored in crs.geo

#Source: #https://www.samuelbosch.com/2014/02/creating-kernel-density-estimate-map-in.html

library("KernSmooth")
library("raster")
library("dplyr")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

points = read.table("Pagani_Scheinfeldt_ID_list_excludeHadBan.txt", header = T)

# compute the 2D binned kernel density estimate
coordinates <- points[,3:4]
coordinates2 <- distinct(coordinates)

est <- bkde2D(coordinates, 
              bandwidth=c(2,2), 
              gridsize=c(3300, 2520),
              range.x=list(c(21.5, 42.5),c(-5, 22.5)))

est <- bkde2D(coordinates2, 
              bandwidth=c(2,2), 
              gridsize=c(1800, 2460),
              range.x=list(c(29.5, 44.5),c(-7.5, 13.0)))

# create raster
est.raster = raster(list(x=est$x1,y=est$x2, z=est$fhat))
projection(est.raster) <- crs.geo
xmin(est.raster) <- 21.5
xmax(est.raster) <- 42.5
ymin(est.raster) <- -5
ymax(est.raster) <- 22.5

xmin(est.raster) <- 29.5
xmax(est.raster) <- 44.5
ymin(est.raster) <- -7.5
ymax(est.raster) <- 13.0

# visually inspect and save the raster output
plot(est.raster)

writeRaster(est.raster, "kernel", format = "GTiff")

