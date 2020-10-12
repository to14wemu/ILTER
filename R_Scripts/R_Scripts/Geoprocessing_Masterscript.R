### Geoprocessing of the GLOBAL PARAMETERS, i.e. the input data for ILTER representativity analysis ###
rm(list = ls())
write("TMPDIR = 'S:/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

library(raster)
library(rgdal)
library(foreign)
#####

#         Geoprocessing Anthromes         #
Anthromes                  <- raster("S:/GIS_data_original/anthromes/2015AD_anthromes/anthromes2015AD.asc")
crs(Anthromes)             <- '+proj=longlat +datum=WGS84'
mollproj                   <- '+proj=moll +ellps=WGS84'
Anthromes                  <- projectRaster(Anthromes, crs = mollproj, method = "ngb") #reprojecting to Mollweide
Anthromes[Anthromes == 70] <- NA  #remove NoData category 70
writeRaster(Anthromes, "S:/GIS_data_processed/anthromes/anthromes2015_mollweide.tif", overwrite = T)

##############
#         Geoprocessing Biomes        #
Biome_shp   <- readOGR(dsn=path.expand("S:/GIS_data_original/Biomes"), layer = "Ecoregions2017")
Anthromes   <- raster("S:/GIS_data_processed/anthromes/anthromes2015_mollweide.tif")
Biome_shp   <- spTransform(Biome_shp, CRSobj = crs(Anthromes))
Biome       <- rasterize(Biome_shp, Anthromes, field = "BIOME_NUM") #transform shapefile to raster
plot(Biome)
writeRaster(Biome, "S:/GIS_data_processed/Biomes/Biome.tif")
Biome       <- mask(Biome, Anthromes,
                    filename = "S:/GIS_data_processed/Biomes/Biome_noA.tif") #removing Antarctica

###############
#       Geoprocessing Bioclimate        #
setwd("S:/GIS_data_original/Bioclimate/GEnSv3")
Bioclimate                    <- raster("gens_v3.tif")       
mollproj                      <- '+proj=moll +ellps=WGS84'
Bioclimate                    <- projectRaster(Bioclimate, crs = mollproj, method = "ngb")
writeRaster(Bioclimate, "S:/GIS_data_processed/Bioclimate/Bioclimate_Mollweide.tif", overwrite = T)

Bioclimate_shp_table          <- read.csv("S:/Analysis/Bioclimate/Bioclimate_shp_table.csv") 
#Attribute table of shapefile, exported from ArcGIS
Bioclimate_shp_table          <- Bioclimate_shp_table[,c(2,6)]
#Reclassifying the 125 strata to 18 global environmental zones
Bioclimate                    <- reclassify(Bioclimate, as.matrix(Bioclimate_shp_table), include.lowest = T, 
                                            filename = "S:/GIS_data_processed/Bioclimate/Bioclimate_MW_rec.tif")
plot(Bioclimate)
Bioclimate                    <- raster("S:/GIS_data_processed/Bioclimate/Bioclimate_MW_rec.tif")
#Resampling: get same origin, extent and resolution as anthromes data
Bioclimate                    <- resample(Bioclimate, Anthromes, method = "ngb",
                                          filename = "S:/GIS_data_processed/Bioclimate/Bioclimate_MW_res.tif",
                                          overwrite = T)
plot(Bioclimate)

##############
#         Geoprocessing Economic Density        #
setwd("S:/GIS_data_original/Economic_Density/GriddedPopulation/gpw-v4-population-density-rev11_2020_15_min_tif")
GriddedPopulation             <- raster("gpw_v4_population_density_rev11_2020_15_min.tif")  
### Population Density for 2020, could also use 2015    Unit: people/km²
mollproj                      <- '+proj=moll +ellps=WGS84'
GriddedPopulation             <- projectRaster(GriddedPopulation, crs = mollproj, method = "ngb")

### Economic Power = GDP per Capita
setwd("S:/GIS_data_original/Economic_Density/GriddedGDP")
GriddedGDP                    <- raster("GDP_per_capita_PPP_1990_2015_v2.nc")  # GDP for years 1990 - 2015 in 2011 international US$
GriddedGDP                    <- projectRaster(GriddedGDP, crs = mollproj, method = "ngb")
plot(GriddedGDP)
### Economic Density
# Defined as GDP per Capita times Population Density
GriddedPopulation             <- resample(GriddedPopulation, GriddedGDP)
EconomicDensity               <- raster::overlay(GriddedPopulation, GriddedGDP, fun = function(r1,r2){return(r1 * r2)})  # Unit: GDP in US§/km²
plot(EconomicDensity)
writeRaster(EconomicDensity, "S:/GIS_data_processed/Economic_Density/EconomicDensity_Mollweide.tif", overwrite = T)

EconomicDensity               <- raster("S:/GIS_data_processed/Economic_Density/EconomicDensity_Mollweide.tif")
EconomicDensity               <- resample(EconomicDensity, Anthromes, method = "bilinear")
EconomicDensity               <- cut(EconomicDensity, breaks = c(-1, 10000, 100000, 1000000, 10000000, 30000000, 10000000000)) # 8 breaks produce 7 classes
plot(EconomicDensity)
writeRaster(EconomicDensity, "S:/GIS_data_processed/Economic_Density/EconomicDensity_MW_res_30M.tif", overwrite = T)

##############
#         Geoprocessing Landforms         #
setwd("S:/GIS_data_original/Landforms")
Landforms                <- raster("EcoLandforms2017.tif")
Landforms                <- projectRaster(Landforms, crs = mollproj, method = "ngb", 
                                        filename = "S:/GIS_data_processed/Landforms/Landforms_Moll.tif",
                                        overwrite = T)
res(Landforms)
Anthromes                <- raster("S:/GIS_data_processed/anthromes/anthromes2015_mollweide.tif")
Landforms                <- resample(Landforms, Anthromes, method = "ngb", 
                                     filename = "S:/GIS_data_processed/Landforms/Landforms_res.tif", overwrite = T)
plot(Landforms)

###############
#       Land Cover      #
## Reprojecting
setwd("S:/GIS_data_original/GlobalLandCover")
LandCover              <- raster("ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif")
plot(LandCover)
mollproj               <- '+proj=moll +ellps=WGS84'
LandCover_rep          <- projectRaster(LandCover, crs = mollproj, method = "ngb")

## Transitional Waters
Land                   <- readOGR(dsn = path.expand("S:/GIS_data_original/GlobalLandCover/ne_50m_land"), 
                                  layer = "ne_50m_land")
Land                   <- spTransform(Land, CRS = crs(LandCover_rep), method = "ngb")
Land                   <- buffer(Land, width = 0, dissolve = T)
plot(Land)
LandCover              <- mask(LandCover_rep, Land) # remove "world ocean"
bufferedLand           <- buffer(Land, width = 3*1852, dissolve = T) # 3 nautical miles buffer
plot(LandCover);plot(bufferedLand, add = T)

TransitionalWaters     <- erase(bufferedLand, Land) # produces transitional waters shapefile
plot(TransitionalWaters)
Transitional_Raster    <- rasterize(TransitionalWaters, LandCover_rep, field = 5) # produces raster with value 5
plot(Transitional_Raster)
writeRaster(Transitional_Raster, "S:/GIS_data_processed/GlobalLandCover/TransitionalWaters/Trans_Raster.tif")

Transitional_Raster    <- overlay(Transitional_Raster, LandCover, fun = function(r1,r2){return (r1+r2)})
plot(Transitional_Raster)
Transitional_Raster[Transitional_Raster != 215] <- NA # remove all categories besides 215 = transitional waters
plot(Transitional_Raster)

LandCover_trans        <- raster::merge(Transitional_Raster, LandCover,  filename = "S:/data_spatial/LandCoverGlobal/LandCover_trans.tif",
                                         overwrite = T)
plot(LandCover_trans)
LandCover_trans        <- projectRaster(LandCover_trans, crs = mollproj, method = "ngb")
writeRaster(LandCover_trans, "S:/GIS_data_processed/GlobalLandCover/LandCover_trans_MW.tif", overwrite = T)

## Reclassifying categories 
reclass_table          <- as.matrix(data.frame(from = c(80,81,82,200,201,202,71,72,150,151,152,153),
                                        to   = c(80,80,80,200,200,200,71,71,150,150,150,150)))
LandCover_rc           <- raster::reclassify(LandCover_trans,  reclass_table, 
                                      filename = "S:/GIS_data_processed/GlobalLandCover/LandCover_rec.tif", 
                                      overwrite = T)

## Removing Black and Caspian Seas
# Done in ArcGIS v. 10.7
# produced with "extract by mask" of the LandCover_trans raster and the land with 3 nm buffer shapefile
# Resulting file: LandCover.tif

## For the processing saving result after every step was important
## with the finalized LandCover-Raster most of the intermediate results have been deleted

## Resampling to Anthromes
setwd("S:/GIS_data_processed/GlobalLandCover")
LandCover              <- raster("S:/GIS_data_processed/GlobalLandCover/LandCover.tif")
LandCover              <- resample(LandCover, Anthromes, method = "ngb",
                                   filename = "S:/GIS_data_processed/GlobalLandCover/LandCover_res.tif",
                                   overwrite = T)
LandCover              <- mask(LandCover, Anthromes, 
                               filename = "S:/GIS_data_processed/GlobalLandCover/LandCover_res_noA.tif")
plot(LandCover)
