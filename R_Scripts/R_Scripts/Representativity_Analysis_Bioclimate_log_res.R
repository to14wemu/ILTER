## Analyzing Geographic Representedness of the ILTER network for the Bioclimate dataset ##

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

# 2.2 Global Bioclimate
Bioclimate                             <- raster("S:/GIS_data_processed/Bioclimate/Bioclimate_MW_res.tif")
ILTER_acc_poly                         <- spTransform(ILTER_acc_poly, CRSobj = crs(Bioclimate))
#ILTER_points                           <- spTransform(ILTER_points, CRSobj = crs(Bioclimate))
ILTER_point_buff                       <- spTransform(ILTER_point_buff, CRSobj = crs(Bioclimate))
plot(Bioclimate);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

ILTER_acc_poly$IntendedPercent_Poly    <- log10(raster::area(ILTER_acc_poly))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
ILTER_point_buff$IntendedPercent_Poly  <- log10(raster::area(ILTER_point_buff))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
# intendedPercent is the weight that every polygon shall have in the representativity analysis based on LOG distribution

sum(ILTER_acc_poly$IntendedPercent_Poly) + sum(ILTER_point_buff$IntendedPercent_Poly)
# Build-in check: value needs to be 100 or very close to 100
# 3. Analysis
#   3.1 Extraction of values
#     3.1.1    Global Distribution
BioClim_Global         <- freq(Bioclimate)
BioClim_Global
BioClim_Global         <- data.frame(BioClim_Global)
BioClim_Global         <- BioClim_Global[1:(length(BioClim_Global$value)-1),]  #remove NA values
BioClim_Global$percent <- BioClim_Global$count*100/sum(BioClim_Global$count)
BioClim_Global
write.csv(BioClim_Global, "S:/Analysis/Bioclimate/BioClim_Global_freq.csv")

#     3.1.2     Distribution on ILTER sites
## Polygon Sites ##
BioClim_ILTER                        <- raster::extract(Bioclimate, ILTER_acc_poly, small = T, df = T) 
count                                <- count(BioClim_ILTER, ID)
ILTER_acc_poly$ID                    <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
BioClim_ILTER                        <- full_join(BioClim_ILTER, count, "ID")
BioClim_ILTER                        <- full_join(BioClim_ILTER, ILTER_acc_poly@data, "ID")
BioClim_ILTER$IntendedPercent_Cell   <- BioClim_ILTER$IntendedPercent_Poly/BioClim_ILTER$n 
sum(BioClim_ILTER$IntendedPercent_Cell, na.rm = T) # check: value ~ 70 - 71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be 
# applied for intendedPercent to be reached for representativity analysis

ILTER_poly           <- BioClim_ILTER %>% group_by(Bioclimate_MW_res) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat            <- count(BioClim_ILTER, Bioclimate_MW_res)
ILTER_poly           <- left_join(ILTER_poly, count_cat, "Bioclimate_MW_res")

colnames(ILTER_poly) <- c("Bioclimate", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly           <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),]

## Point Sites ##
BioClim_ILTER_point                        <- raster::extract(Bioclimate, ILTER_point_buff, small = T, df = T) 
count                                      <- count(BioClim_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                        <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
BioClim_ILTER_point                        <- full_join(BioClim_ILTER_point, count, "ID")
BioClim_ILTER_point                        <- full_join(BioClim_ILTER_point, ILTER_point_buff@data, "ID")
BioClim_ILTER_point$IntendedPercent_Cell   <- BioClim_ILTER_point$IntendedPercent_Poly/BioClim_ILTER_point$n 
sum(BioClim_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: has to be ~ 29-30

ILTER_point               <- BioClim_ILTER_point %>% group_by(Bioclimate_MW_res) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(BioClim_ILTER_point, Bioclimate_MW_res)
ILTER_point               <- left_join(ILTER_point, count_cat, "Bioclimate_MW_res")

colnames(ILTER_point)     <- c("Bioclimate", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),] #remove NA

## Combining Polygon and Point Site Data ##
ILTER                           <- full_join(ILTER_poly, ILTER_point, "Bioclimate")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$vis_n                     <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(ILTER, "S:/Analysis/Bioclimate/BioClim_ILTER_categorized.csv")       
write.csv(BioClim_ILTER, "S:/Analysis/Bioclimate/BioClim_ILTER.csv")  
#####
#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
# Tidy tables
df_BioClim           <- full_join(BioClim_Global[,1:2], ILTER[,c(1,8)], by = c( "value" = "Bioclimate"))
colnames(df_BioClim) <- c("value", "count_global", "count_ilter")
# Produce matrix for chisq.test
BioClim_vis                     <- df_BioClim                         # Necessary for subsequent visualization
BioClim_chi                     <- t(df_BioClim)                     # t() also produces matrix
colnames(BioClim_chi)           <- BioClim_chi[1,]
BioClim_chi                     <- BioClim_chi[-1,]
BioClim_chi[is.na(BioClim_chi)] <- 0
chisq                      <- chisq.test(BioClim_chi)                  
chisq

#####
# Calculating Geographic Representedness
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
BioClim_vis <- data.frame(BioClim_vis, rep)
BioClim_vis$count_ilter[is.na(BioClim_vis$count_ilter)] <- 0
BioClim_vis$rep[BioClim_vis$count_ilter == 0] <- NA
Bioclimate_shp_table <- read.csv("S:/Analysis/Bioclimate/Bioclimate_shp_table.csv")
BioClim_vis$Category <- levels(Bioclimate_shp_table$GEnZname)
write.table(BioClim_vis, "S:/Analysis/Bioclimate/Bioclimate_vis.csv")
# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)

Representativity_Bioclimate<-
  ggplot(BioClim_vis, aes(x = as.factor(BioClim_vis$value))) + 
  geom_bar(aes(y = BioClim_vis$count_global/sum(BioClim_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = BioClim_vis$count_ilter/sum(BioClim_vis$count_ilter, na.rm = T)*100, color = BioClim_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Bioclimate Category based on Metzger et al. (2013)", 
       y = "Relative Spatial Coverage [%]", color = "Geographic\nRepresentedness") +
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1),
                        labels = c("Underrepresented", "well represented", "Overrepresented")) +
    scale_x_discrete(labels = BioClim_vis$Category) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=14, xmin = 12)

  ggsave(filename = "Representativity_Bioclimate_1000.jpg", plot = Representativity_Bioclimate,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geographic Representedness
BioClim_vis$rep[is.na(BioClim_vis$rep)] <- -1 # Asssign underrepresented value to categories not present on ILTER sites
Rep_Bioclimate <- reclassify(Bioclimate, BioClim_vis[,c(1,4)], filename = "S:/Results/rep_maps/Rep_BioClim.tif",
                              overwrite = T) 
plot(Rep_Bioclimate) 
