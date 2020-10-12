## Analyzing Geographic Representedness of the ILTER network for the Landforms dataset ##
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
# 2.1 Accreted ILTER sites
# 2.1 Accredited ILTER sites
ILTER_acc_poly  <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                           layer = "ilter_boundaries")
ILTER_point_buff<- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
                           layer = "Points_with_buffer")

# Point site is only one marine site, excluded from the analysis
#ILTER_points    <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
#                          layer = "Points_without_size")

# 2.2 Landform Units
Landforms                              <- raster("S:/GIS_data_processed/Landforms/Landforms_res.tif")
ILTER_acc_poly                         <- spTransform(ILTER_acc_poly, CRSobj = crs(Landforms))
#ILTER_points                           <- spTransform(ILTER_points, CRSobj = crs(Landforms))
ILTER_point_buff                       <- spTransform(ILTER_point_buff, CRSobj = crs(Landforms))
plot(Landforms);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

ILTER_acc_poly$IntendedPercent_Poly    <- log10(raster::area(ILTER_acc_poly))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
ILTER_point_buff$IntendedPercent_Poly  <- log10(raster::area(ILTER_point_buff))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
# intendedPercent is the weight that every polygon shall have in the representativity analysis
# for analysis of LOG10 distribution see CumulativeArea_ILTER_Sites.R-Script

sum(ILTER_acc_poly$IntendedPercent_Poly) + sum(ILTER_point_buff$IntendedPercent_Poly)
# Build-in check: value needs to be 100 or very close to 100

# 3. Analysis
#   3.1 Extraction of values
#     3.1.1    Global Distribution
Landforms                       <- reclassify(Landforms, cbind(0,0, NA), include.lowest = T) # remove ocean
Landforms_Global                <- freq(Landforms)
Landforms_Global
Landforms_Global                <- data.frame(Landforms_Global)
Landforms_Global                <- Landforms_Global[1:(length(Landforms_Global$value)-1),]
Landforms_Global$percent        <- Landforms_Global$count*100/sum(Landforms_Global$count)
Landforms_Global
write.csv(Landforms_Global, "S:/Analysis/Landforms/Landforms_Global_freq.csv")

#     3.1.2    Weighted Distribution on ILTER sites
## Polygon sites ##
Landforms_ILTER                         <- raster::extract(Landforms, ILTER_acc_poly, small = T, df = T)
count                                   <- count(Landforms_ILTER, ID)
ILTER_acc_poly$ID                       <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
Landforms_ILTER                         <- full_join(Landforms_ILTER, count, "ID")
Landforms_ILTER                         <- full_join(Landforms_ILTER, ILTER_acc_poly@data, "ID")
Landforms_ILTER$IntendedPercent_Cell    <- Landforms_ILTER$IntendedPercent_Poly/Landforms_ILTER$n 
sum(Landforms_ILTER$IntendedPercent_Cell, na.rm = T) # build in check: should be ~ 70 - 71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be 
# applied for intendedPercent to be reached for representativity analysis

ILTER_poly               <- Landforms_ILTER %>% group_by(Landforms_res) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                <- count(Landforms_ILTER, Landforms_res)
ILTER_poly               <- left_join(ILTER_poly, count_cat, "Landforms_res")

colnames(ILTER_poly)     <- c("Landform", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n     <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly               <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),] # remove NA

## Point sites ##
Landforms_ILTER_point                        <- raster::extract(Landforms, ILTER_point_buff, small = T, df = T) 
count                                        <- count(Landforms_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                          <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
Landforms_ILTER_point                        <- full_join(Landforms_ILTER_point, count, "ID")
Landforms_ILTER_point                        <- full_join(Landforms_ILTER_point, ILTER_point_buff@data, "ID")
Landforms_ILTER_point$IntendedPercent_Cell   <- Landforms_ILTER_point$IntendedPercent_Poly/Landforms_ILTER_point$n 
sum(Landforms_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: should be ~ 29-30

ILTER_point               <- Landforms_ILTER_point %>% group_by(Landforms_res) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(Landforms_ILTER_point, Landforms_res)
ILTER_point               <- left_join(ILTER_point, count_cat, "Landforms_res")

colnames(ILTER_point)     <- c("Landform", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),] # remove NA

## Combining Polygon and Point Site Data ##
ILTER                           <- full_join(ILTER_poly, ILTER_point, "Landform")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$desired_n.x[is.na(ILTER$desired_n.x)] <- 0
ILTER$vis_n                     <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(ILTER, "S:/Analysis/Landforms/Landforms_ILTER_categorized.csv")             
write.csv(Landforms_ILTER, "S:/Analysis/Landforms/Landforms_ILTER.csv") 

#####
#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
df_Landforms              <- full_join(Landforms_Global[,1:2], ILTER[,c(1,8)], by = c( "value" = "Landform"))
colnames(df_Landforms)    <- c("value", "count_global", "count_ilter")
Landforms_vis             <- df_Landforms                        # used for visualization
Landforms_chi             <- t(df_Landforms)                     # t() also produces matrix
colnames(Landforms_chi)   <- Landforms_chi[1,]
Landforms_chi             <- Landforms_chi[-1,]
Landforms_chi[is.na(Landforms_chi)] <- 0
chisq <- chisq.test(Landforms_chi)                  
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
Landforms_vis               <- data.frame(Landforms_vis, rep)
Landforms_vis$count_ilter[is.na(Landforms_vis$count_ilter)] <- 0
Landforms_vis$rep[Landforms_vis$count_ilter == 0] <- NA
Landforms_vis$category_name <- c("Nearly Flat Plains", "Smooth Plains With Some Local Relief", "Irregular Plains With Moderate Relief", 
                                 "Irregular Plains With Low Hills", "Scattered Moderate Hills", "Scattered High Hills", 
                                 "Scattered Low Mountains", "Scattered High Mountains", "Moderate Hills", "High Hills",
                                 "Tablelands With Moderate Relief", "Tableland with Considerable Relief", "Tablelands With High Relief", 
                                 "Tablelands With Very High Relief", "Low Mountains", "High Mountains")
write.table(Landforms_vis, "S:/Analysis/Landforms/Landforms_vis.csv")

# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)

Representativity_Landforms<-
  ggplot(Landforms_vis, aes(x = reorder(Landforms_vis$category_name, Landforms_vis$value))) + 
  geom_bar(aes(y = Landforms_vis$count_global/sum(Landforms_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = Landforms_vis$count_ilter/sum(Landforms_vis$count_ilter, na.rm = T)*100, color = Landforms_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Landforms based on Karagulle et al. (2017)", y = "Relative Spatial Coverage [%]", 
       color = "Geographic\nRepresentedness") + 
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1), 
                        labels = c("Underrepresented", "Well represented", "Overrepresented"))+
  scale_x_discrete(labels = Landforms_vis$category_name) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  
 # annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=35, xmin = 12) 
ggsave(filename = "Representativity_Landforms_1000.jpg", plot = Representativity_Landforms,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geographic Representedness
Landforms_vis$rep[is.na(Landforms_vis$rep)] <- -1 #set values to underrepresented for categories not present on ILTER sites
Rep_Landforms <- reclassify(Landforms, Landforms_vis[,c(1,4)], 
                            filename = "S:/Results/rep_maps/Rep_Landforms.tif", overwrite = T) 
plot(Rep_Landforms) 
