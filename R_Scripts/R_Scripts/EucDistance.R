# Producing Euclidean Distance Raster from the slightly flawed polygon-Euclidean Distance raster and the
# point-Euclidean Distance raster, both computed in ArcGIS v. 10.8
rm(list = ls())

library(raster)
library(rgdal)

euc_point  <- raster("S:/GIS_data_processed/EuclideanDistance/euc_point_ter")
plot(euc_point)
euc_poly   <- raster("S:/GIS_data_processed/EuclideanDistance/euc_terr")
plot(euc_poly)
euc_poly   <- resample(euc_poly, euc_point)
euc        <- overlay(euc_point, euc_poly, fun = min)
plot(euc)

writeRaster(euc, "S:/GIS_data_processed/EuclideanDistance/euc_point_poly.tif")
