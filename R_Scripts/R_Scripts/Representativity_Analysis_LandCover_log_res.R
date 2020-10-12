## Analyzing Geographic Representedness of the ILTER network for the Land Cover dataset ##
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
#detach("package:plyr", unload = T)
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

# 2.2 Global Land Cover
LandCover              <- raster("S:/GIS_data_processed/GlobalLandCover/LandCover_res.tif")
Anthromes              <- raster("S:/GIS_data_processed/Anthromes/anthromes2015_mollweide.tif")
LandCover              <- mask(LandCover, Anthromes) # removing Antarctica

ILTER_acc_poly         <- spTransform(ILTER_acc_poly, CRSobj = crs(Anthromes))
#ILTER_points           <- spTransform(ILTER_points, CRSobj = crs(Anthromes))
ILTER_point_buff       <- spTransform(ILTER_point_buff, CRSobj = crs(Anthromes))
plot(LandCover);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

ILTER_acc_poly$IntendedPercent_Poly    <- log10(raster::area(ILTER_acc_poly))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
ILTER_point_buff$IntendedPercent_Poly  <- log10(raster::area(ILTER_point_buff))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
# intendedPercent is the weight that every polygon shall have in the representativity analysis
# for analysis of LOG10 distribution see CumulativeArea_ILTER_Sites.R-Script

sum(ILTER_acc_poly$IntendedPercent_Poly) + sum(ILTER_point_buff$IntendedPercent_Poly)
# Build-in check: value has to be 100 or very close to 100
#####

# 3. Analysis
#   3.1 Extraction of values
#     3.1.1 Global Distribution
LandCover_Global                <- freq(LandCover)
LandCover_Global
LandCover_Global                <- data.frame(LandCover_Global)
LandCover_Global                <- LandCover_Global[1:(length(LandCover_Global$value)-1),]
LandCover_Global$percent        <- LandCover_Global$count*100/sum(LandCover_Global$count)
LandCover_Global
write.csv(LandCover_Global, "S:/Analysis/LandCover/LandCover_Global_freq.csv")

#     3.1.2 Distribution on ILTER sites
## Polygon Sites  ##
LandCover_ILTER                        <- raster::extract(LandCover, ILTER_acc_poly, small = T, df = T)
count                                  <- count(LandCover_ILTER, ID)
ILTER_acc_poly$ID                      <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
LandCover_ILTER                        <- full_join(LandCover_ILTER, count, "ID")
LandCover_ILTER                        <- full_join(LandCover_ILTER, ILTER_acc_poly@data, "ID")
LandCover_ILTER$IntendedPercent_Cell   <- LandCover_ILTER$IntendedPercent_Poly/LandCover_ILTER$n 
sum(LandCover_ILTER$IntendedPercent_Cell, na.rm = T) # check: value ~ 70-71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be 
# applied for intendedPercent to be reached for representativity analysis

ILTER_poly <- LandCover_ILTER %>% group_by(LandCover_res) %>% summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight  ####sum(n, na.rm = T),
count_cat  <- count(LandCover_ILTER, LandCover_res)
ILTER_poly <- left_join(ILTER_poly, count_cat, "LandCover_res")

colnames(ILTER_poly) <- c("LandCover", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly           <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),] # remove NA

## Point sites ##
LandCover_ILTER_point                        <- raster::extract(LandCover, ILTER_point_buff, small = T, df = T) 
count                                        <- count(LandCover_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                          <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
LandCover_ILTER_point                        <- full_join(LandCover_ILTER_point, count, "ID")
LandCover_ILTER_point                        <- full_join(LandCover_ILTER_point, ILTER_point_buff@data, "ID")
LandCover_ILTER_point$IntendedPercent_Cell   <- LandCover_ILTER_point$IntendedPercent_Poly/LandCover_ILTER_point$n 
sum(LandCover_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: should be ~ 29-30

ILTER_point               <- LandCover_ILTER_point %>% group_by(LandCover_res) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(LandCover_ILTER_point, LandCover_res)
ILTER_point               <- left_join(ILTER_point, count_cat, "LandCover_res")

colnames(ILTER_point)     <- c("LandCover", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),] #remove NA

## Combining Polygon and Point Site Data ##
ILTER                           <- full_join(ILTER_poly, ILTER_point, "LandCover")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$vis_n                     <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(LandCover_ILTER, "S:/Analysis/LandCover/LandCover_ILTER.csv")
write.csv(ILTER, "S:/Analysis/LandCover/Anthromes_ILTER_categorized.csv") 

#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
# Tidy tables
df_LandCover           <- full_join(LandCover_Global[,1:2], ILTER[,c(1,8)], by = c( "value" = "LandCover"))
colnames(df_LandCover) <- c("value", "count_global", "count_ilter")


# Produce matrix for chisq.test
LandCover_vis                       <- df_LandCover                         # Necessary for subsequent visualization
LandCover_chi                       <- t(df_LandCover)                     # t() also produces matrix
colnames(LandCover_chi)             <- LandCover_chi[1,]
LandCover_chi                       <- LandCover_chi[-1,]
LandCover_chi[is.na(LandCover_chi)] <- 0
chisq                               <- chisq.test(LandCover_chi)                  
chisq

#####
# Computing Geographic Representedness
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
LandCover_vis                                               <- data.frame(LandCover_vis, rep)
LandCover_vis$count_ilter[is.na(LandCover_vis$count_ilter)] <- 0
LandCover_vis$rep[LandCover_vis$count_ilter == 0]           <- NA
LandCover_vis$category <- 	c("Cropland: Rainfed", "Herbaceous Cover", "Tree or Shrub Cover",
"Cropland: Irrigated/Post-Flooding", "Mosaic: > 50% Cropland",
"Mosaic: < 50% Cropland",  "Tree Cover: Broadleafed Evergreen",
"Tree Cover: Broadleafed Deciduous < 15%",
"Tree Cover: Broadleafed Deciduous > 40%","Tree Cover: Broadleafed Deciduous 15 - 40%",
"Tree Cover: Needleleafed Evergreen < 15%","Tree Cover: Needleleafed Evergreen > 15%","Tree Cover: Needleleafed Deciduous",
"Tree Cover: Mixed Leaf", "Mosaic: > 50% Tree and Shrub",
"Mosaic: > 50% Herbaceous", "Shrubland", "Shrubland Evergreen", "Shrubland Deciduous",
"Grassland", "Lichens and Mosses", "Sparse Vegetation", "Flooded: Tree - Fresh/Brakish Water", "Flooded: Tree - Saline Water",
"Flooded: Shrub/Herbaceous", "Urban Areas", "Bare Areas", 
"Water Bodies", "Transitional Waters", "Permanent Snow and Ice")
write.table(LandCover_vis, "S:/Analysis/LandCover/LandCover_vis_rc.csv")



# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)


Representativity_LandCover <-
ggplot(LandCover_vis, aes(x = reorder(LandCover_vis$category, LandCover_vis$value))) + 
  geom_bar(aes(y = LandCover_vis$count_global/sum(LandCover_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = LandCover_vis$count_ilter/sum(LandCover_vis$count_ilter, na.rm = T)*100, color = LandCover_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Landcover Category based on European Space Agency (2015)", y = "Relative Spatial Coverage [%]", color = "Geographic\nRepresentedness") +
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1),
                        labels = c("Underrepresented", "Well represented", "Overrepresented"))+
  scale_x_discrete(labels = LandCover_vis$category) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=13.2, xmin = 21.5)

ggsave(filename = "Representativity_LandCover_1000.jpg", plot = Representativity_LandCover,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geographic Representedness
LandCover_vis$rep[is.na(LandCover_vis$rep)] <- -1 #set values to underrepresented for categories not present on ILTER sites
Rep_LandCover <- reclassify(LandCover, LandCover_vis[,c(1,4)], 
                            filename = "S:/Results/rep_maps/Rep_LandCover.tif", overwrite = T) 
plot(Rep_LandCover) 
