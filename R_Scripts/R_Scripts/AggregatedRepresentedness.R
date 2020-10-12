# Producing Aggregated Representedness Raster #
rm(list = ls())
library(raster)
library(rgdal)
library(ggplot2)
library(ggspatial)
library(sf)
library(dplyr)
library(ggfortify)
library(rasterVis)
library(lattice)

# Load Data
rep_ED              <- raster("S:/Results/rep_maps/Rep_ED.tif")
rep_Landforms       <- raster("S:/Results/rep_maps/Rep_Landforms.tif")
rep_Anthromes       <- raster("S:/Results/rep_maps/Rep_Anthromes.tif")
rep_BioClim         <- raster("S:/Results/rep_maps/Rep_BioClim.tif")
rep_LandCover       <- raster("S:/Results/rep_maps/Rep_LandCover.tif")
rep_Biome           <- raster("S:/Results/rep_maps/Rep_Biome.tif")

# Compute Aggregated Representedness
Agg_rep       <- rep_ED + rep_Anthromes + rep_Biome + rep_Landforms + rep_BioClim + rep_LandCover
plot(Agg_rep)
Table         <- freq(Agg_rep)
Table         <- data.frame(Table)
Table         <- Table[1:(length(Table$value)-1),] #remove NA
Table$percent <- Table$count/sum(Table$count)*100
Table

writeRaster(Agg_rep, "S:/Results/rep_maps/Aggregated_Representedness.tif", overwrite = T)
