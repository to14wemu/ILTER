## Analysis of sites accreted on 28th of September 2020 ##
rm(list = ls())

#####
# 1. Load Packages
####
library(raster)
library(rgdal)
library(ggplot2)
library(rgeos)
library(epiDisplay)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyr)
library(purrr)
library(ncdf4)
library(sf)
library(spatialEco)
library(plyr)

# 2. Load Data
ILTER_acc_single  <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                            layer = "ilter_boundaries")
mollproj          <- '+proj=moll +ellps=WGS84'
ILTER_acc_points  <-readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                            layer = "ilter_all_formal")
plot(ILTER_acc_points);plot(ILTER_acc_single, add = T)

# 3. Analyses
# 3.1 Identify point data without associated Polygon
colnames(ILTER_acc_points@data)[1:2] <- c("title", "id")
Points_Without_Poly                  <- merge(ILTER_acc_points, ILTER_acc_single@data, by = "id")
Points_Without_Poly                  <- Points_Without_Poly[is.na(Points_Without_Poly@data$title.y),]
length(Points_Without_Poly)
#writeOGR(Points_Without_Poly, dsn = path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
#         layer = "Points_Without_Poly", driver = "ESRI Shapefile")

# 3.2 Assign reported site sizes to point data without polygon representation
ilter_size               <- read.csv("S:/GIS_data_original/ILTER/2020_09_28/ilter_sizes.csv", row.names = NULL)
colnames(ilter_size)[2]  <- "id"
Points_Without_Poly@data <- left_join(Points_Without_Poly@data, ilter_size, by = "id")
Points_with_size         <- Points_Without_Poly[Points_Without_Poly$size_ha > 0,]
length(Points_with_size) # count of point sites that have a size reported
length(Points_with_size[Points_with_size$size_ha > 8600,]) # count of sites that are bigger than one raster cell

crs(Points_with_size) # needs to be long/lat
Points_with_size$radius  <- sqrt((Points_with_size$size_ha*10000/3.14)) # calculate radius for buffer
Points_with_buffer       <- buffer(Points_with_size, width = Points_with_size$radius, dissolve = F)
plot(Points_with_buffer)
Points_with_buffer       <- spTransform(Points_with_buffer, CRSobj = mollproj)
plot(Points_with_buffer)
Points_with_buffer       <- as(Points_with_buffer, "SpatialPolygonsDataFrame")
writeOGR(Points_with_buffer, dsn = path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
         layer = "Points_with_buffer", driver = "ESRI Shapefile", overwrite_layer = T)

Points_without_size      <- Points_Without_Poly[Points_Without_Poly$size_ha == 0,]
length(Points_without_size)
writeOGR(Points_without_size, dsn = path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
         layer = "Points_without_size", driver = "ESRI Shapefile", overwrite_layer = T)

# 3.3 Identify sites > 100,000 km²
ILTER_acc_single         <- spTransform(ILTER_acc_single, CRSobj = mollproj)
ILTER_acc_single$area    <- raster::area(ILTER_acc_single)/1000000 # area in km²
ILTER_acc_single         <- ILTER_acc_single[ILTER_acc_single$area > 100000,]
length(ILTER_acc_single)
Points_100000            <- Points_with_size[Points_with_size$size_ha > 10000000,]
length(Points_100000)
