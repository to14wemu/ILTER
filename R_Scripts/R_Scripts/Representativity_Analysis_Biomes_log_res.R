## Analyzing Geographic Representedness of the ILTER network for the Biome dataset ##
rm(list = ls())
write("TMPDIR = 'S:/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

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

#####
# 2. Load Data
# 2.1 Accredited ILTER sites
ILTER_acc_poly  <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                           layer = "ilter_boundaries")
ILTER_point_buff<- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
                           layer = "Points_with_buffer")

# Point site is only one marine site, excluded from the analysis
#ILTER_points    <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
#                          layer = "Points_without_size")

# 2.2 Global Biomes
Biome                <- raster("S:/GIS_data_processed/Biomes/Biome.tif")
Anthromes            <- raster("S:/GIS_data_processed/Anthromes/anthromes2015_mollweide.tif")
Biome                <- mask(Biome, Anthromes) # remove Antarctica

ILTER_acc_poly                         <- spTransform(ILTER_acc_poly, CRSobj = crs(Anthromes))
#ILTER_points                           <- spTransform(ILTER_points, CRSobj = crs(Anthromes))
ILTER_point_buff                       <- spTransform(ILTER_point_buff, CRSobj = crs(Anthromes))
plot(Biome);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

ILTER_acc_poly$IntendedPercent_Poly    <- log10(raster::area(ILTER_acc_poly))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
ILTER_point_buff$IntendedPercent_Poly  <- log10(raster::area(ILTER_point_buff))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
# intendedPercent is the weight that every polygon shall have in the representativity analysis based on LOG distribution

sum(ILTER_acc_poly$IntendedPercent_Poly) + sum(ILTER_point_buff$IntendedPercent_Poly)
# Build-in check: value has to be 100 or very close to 100

# 3.Analysis
#   3.1 Extraction of values
#     3.1.1    Global Distribution
Biome_Global         <- freq(Biome)
Biome_Global
Biome_Global         <- data.frame(Biome_Global)
Biome_Global         <- Biome_Global[1:(length(Biome_Global$value)-1),]   #remove NA
Biome_Global$percent <- Biome_Global$count*100/sum(Biome_Global$count)
Biome_Global
write.csv(Biome_Global, "S:/Analysis/Biomes/Biome_Global_freq.csv")

#     3.1.2     Distribution on ILTER sites
## Polygon Sites ##
Biome_ILTER                        <- raster::extract(Biome, ILTER_acc_poly, small = T, df = T) # extracting all values in the shapefile extents
count                              <- count(Biome_ILTER, ID) # n is the amount of raster pixels/values in a polygon
ILTER_acc_poly$ID                  <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
Biome_ILTER                        <- full_join(Biome_ILTER, count, "ID")
Biome_ILTER                        <- full_join(Biome_ILTER, ILTER_acc_poly@data, "ID")
Biome_ILTER$IntendedPercent_Cell   <- Biome_ILTER$IntendedPercent_Poly/Biome_ILTER$n 
sum(Biome_ILTER$IntendedPercent_Cell, na.rm = T) # check: value ~ 70 - 71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be 
# applied for intendedPercent to be reached for representativity analysis

ILTER_poly <- Biome_ILTER %>% group_by(Biome) %>% summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight  ####sum(n, na.rm = T),
count_cat  <- count(Biome_ILTER, Biome)
ILTER_poly <- left_join(ILTER_poly, count_cat, "Biome")

colnames(ILTER_poly) <- c("Biome", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly           <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),]

## Point sites ##
Biome_ILTER_point                        <- raster::extract(Biome, ILTER_point_buff, small = T, df = T) 
count                                    <- count(Biome_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                      <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
Biome_ILTER_point                        <- full_join(Biome_ILTER_point, count, "ID")
Biome_ILTER_point                        <- full_join(Biome_ILTER_point, ILTER_point_buff@data, "ID")
Biome_ILTER_point$IntendedPercent_Cell   <- Biome_ILTER_point$IntendedPercent_Poly/Biome_ILTER_point$n 
sum(Biome_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: should be ~ 29-30

ILTER_point               <- Biome_ILTER_point %>% group_by(Biome) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(Biome_ILTER_point, Biome)
ILTER_point               <- left_join(ILTER_point, count_cat, "Biome")

colnames(ILTER_point)     <- c("Biome", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),]

## Combining Polygon and Point Site Data ##
ILTER                           <- full_join(ILTER_poly, ILTER_point, "Biome")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$desired_n.x[is.na(ILTER$desired_n.x)] <- 0
ILTER$vis_n                     <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(ILTER, "S:/Analysis/Biomes/Biome_ILTER_categorized.csv")       
write.csv(Biome_ILTER, "S:/Analysis/Biomes/Biome_ILTER.csv") 

#####
#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
# Tidy tables
df_Biome           <- full_join(Biome_Global[,1:2], ILTER[,c(1,8)], by = c( "value" = "Biome"))
colnames(df_Biome) <- c("value", "count_global", "count_ilter")

# Produce matrix for chisq.test
Biome_vis                     <- df_Biome                         # Necessary for subsequent visualization
Biome_chi                     <- t(df_Biome)                     # t() also produces matrix
colnames(Biome_chi)           <- Biome_chi[1,]
Biome_chi                     <- Biome_chi[-1,]
Biome_chi[is.na(Biome_chi)]   <- 0
chisq                         <- chisq.test(Biome_chi)                  
chisq

#####
# Accounting for representedness 
obs <- data.frame(chisq$observed)
exp <- data.frame(chisq$expected)

rep <- c()            # shall in the end contain value for representedness for every category
for (i in 1:length(obs)){
  mat    <- matrix(c(obs[2,i], exp[2,i], sum(obs[2,])-obs[2,i], sum(exp[2,])-exp[2,i]), nrow=2)
  chi    <- chisq.test(mat)
  rep[i] <- 
    ifelse(obs[2,i] == exp[2,i], 0,
           ifelse(obs[2,i] <  exp[2,i], -1+chi$p.value,
                  ifelse(obs[2,i] >  exp[2,i], 1-chi$p.value,
                         ifelse(obs[2,i] > 0 & exp[2,i] == 0, -999, -1000))))
}
rep

#   3.3 Visualization
## Histogram
Biome_vis <- data.frame(Biome_vis, rep)
Biome_vis$count_ilter[is.na(Biome_vis$count_ilter)] <- 0
Biome_vis$rep[Biome_vis$count_ilter == 0] <- NA
Biome_vis$category_name <- c("Tropical & Subtropical Moist Broadleaf Forests",
"Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Coniferous Forests","Temperate Broadleaf & Mixed Forests",
"Temperate Conifer Forests","Boreal Forests or Taiga","Tropical & Subtropical Grasslands, Savannas & Shrublands",
"Temperate Grasslands, Savannas & Shrublands","Flooded Grasslands & Savannas",
"Montane Grasslands & Shrublands","Tundra","Mediterranean Forests, Woodlands & Scrub","Deserts & Xeric Shrublands","Mangroves")
write.table(Biome_vis, "S:/Analysis/Biomes/Biome_vis.csv")

# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)

Representativity_Biomes <-
  ggplot(Biome_vis, aes(x = reorder(Biome_vis$category, Biome_vis$value))) + 
  geom_bar(aes(y = Biome_vis$count_global/sum(Biome_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = Biome_vis$count_ilter/sum(Biome_vis$count_ilter, na.rm = T)*100, color = Biome_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Biome Category based on Dinerstein et al. (2017)", y = "Relative Spatial Coverage [%]", color = "Representedness") +
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1),
                        labels = c("Underrepresented", "Well represented", "Overrepresented")) +
  scale_x_discrete(labels = Biome_vis$category_name) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=25, xmin = 10)

ggsave(filename = "Representativity_Biomes_1000.jpg", plot = Representativity_Biomes,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geograpic Representedness
Rep_Biome <- reclassify(Biome, Biome_vis[,c(1,4)], filename = "S:/Results/rep_maps/Rep_Biome.tif", overwrite = T) 
plot(Rep_Biome) 
