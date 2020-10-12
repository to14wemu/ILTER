## Analyzing Geographic Representedness of the ILTER network for the anthrome dataset ##
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

# 2.2 Global Anthromes
Anthromes                              <- raster("S:/GIS_data_processed/anthromes/anthromes2015_mollweide.tif")
ILTER_acc_poly                         <- spTransform(ILTER_acc_poly, CRSobj = crs(Anthromes))
#ILTER_points                           <- spTransform(ILTER_points, CRSobj = crs(Anthromes))
ILTER_point_buff                       <- spTransform(ILTER_point_buff, CRSobj = crs(Anthromes))
plot(Anthromes);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

ILTER_acc_poly$IntendedPercent_Poly      <- log10(raster::area(ILTER_acc_poly))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
ILTER_point_buff$IntendedPercent_Poly    <- log10(raster::area(ILTER_point_buff))/
  (sum(log10(raster::area(ILTER_acc_poly))) + sum(log10(raster::area(ILTER_point_buff))))*100
# intendedPercent is the weight that every polygon shall have in the representativity analysis based on LOG distribution

sum(ILTER_acc_poly$IntendedPercent_Poly) + sum(ILTER_point_buff$IntendedPercent_Poly)
# Build-in check: value has to be 100 or very close to 100

# 3. Analysis
#   3.1 Extraction of values
#     3.1.1    Global Distribution
Anthromes_Global         <- freq(Anthromes)
Anthromes_Global
Anthromes_Global         <- data.frame(Anthromes_Global)
Anthromes_Global         <- Anthromes_Global[1:(length(Anthromes_Global$value)-1),] #remove NA values
Anthromes_Global$percent <- Anthromes_Global$count*100/sum(Anthromes_Global$count)
Anthromes_Global
write.csv(Anthromes_Global, "S:/Analysis/Anthromes/Anthromes_Global_freq.csv")

#     3.1.2    Weighted Distribution on ILTER sites
## Polygon sites ##
Anthromes_ILTER                        <- raster::extract(Anthromes, ILTER_acc_poly, small = T, df = T) 
count                                  <- count(Anthromes_ILTER, ID) # n is the amount of raster pixels/values in a polygon
ILTER_acc_poly$ID                      <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
Anthromes_ILTER                        <- full_join(Anthromes_ILTER, count, "ID")
Anthromes_ILTER                        <- full_join(Anthromes_ILTER, ILTER_acc_poly@data, "ID")
Anthromes_ILTER$IntendedPercent_Cell   <- Anthromes_ILTER$IntendedPercent_Poly/Anthromes_ILTER$n 
sum(Anthromes_ILTER$IntendedPercent_Cell, na.rm = T) # build in check: has to be ~ 70-71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs 
# to be applied for intendedPercent to be reached for representativity analysis

ILTER_poly               <- Anthromes_ILTER %>% group_by(anthromes2015_mollweide) %>% 
  summarize( sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                <- count(Anthromes_ILTER, anthromes2015_mollweide)
ILTER_poly               <- left_join(ILTER_poly, count_cat, "anthromes2015_mollweide")

colnames(ILTER_poly)     <- c("anthromes", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n     <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly               <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),]

## Point sites ##
Anthromes_ILTER_point                        <- raster::extract(Anthromes, ILTER_point_buff, small = T, df = T) 
count                                        <- count(Anthromes_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                          <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
Anthromes_ILTER_point                        <- full_join(Anthromes_ILTER_point, count, "ID")
Anthromes_ILTER_point                        <- full_join(Anthromes_ILTER_point, ILTER_point_buff@data, "ID")
Anthromes_ILTER_point$IntendedPercent_Cell   <- Anthromes_ILTER_point$IntendedPercent_Poly/Anthromes_ILTER_point$n 
sum(Anthromes_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: has to be ~ 29-30

ILTER_point               <- Anthromes_ILTER_point %>% group_by(anthromes2015_mollweide) %>% 
                              summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(Anthromes_ILTER_point, anthromes2015_mollweide)
ILTER_point               <- left_join(ILTER_point, count_cat, "anthromes2015_mollweide")

colnames(ILTER_point)     <- c("anthromes", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),] #remove NA

## Combining Polygon and Point Site Data ##
ILTER                     <- full_join(ILTER_poly, ILTER_point, "anthromes")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$vis_n               <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(ILTER, "S:/Analysis/Anthromes/Anthromes_ILTER_categorized.csv")             
write.csv(Anthromes_ILTER, "S:/Analysis/Anthromes/Anthromes_ILTER.csv")       
#####
#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
df_Anthromes              <- full_join(Anthromes_Global[,1:2], ILTER[,c(1,8)], by = c("value" = "anthromes"))
colnames(df_Anthromes)    <- c("value", "count_global", "count_ilter")

Anthromes_vis             <- df_Anthromes                        # used for visualization
Anthromes_chi             <- t(df_Anthromes)                     # t() also produces matrix
colnames(Anthromes_chi)   <- Anthromes_chi[1,]
Anthromes_chi             <- Anthromes_chi[-1,]
chisq                     <- chisq.test(Anthromes_chi)                  
chisq

#####
# Calculate Geographic Representedness 
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
Anthromes_vis               <- data.frame(Anthromes_vis, rep)
Anthromes_vis$rep[Anthromes_vis$count.ilter == 0] <- NA
Anthromes_vis$category_name <- c("Urban","Mixed Settlements","Rice Villages","Irrigated Villages",
                                 "Rainfed Villages","Pastoral Villages","Residential Irrigated Croplands",
                                 "Residential Rainfed Croplands","Populated Croplands","Remote Croplands",
"Residential Rangelands","Populated Rangelands","Remote Rangelands",
"Residential Woodlands","Popoulated Woodlands","Remote Woodlands","Inhabited Treeless & Barren Lands",
"Wild Woodlands", "Wild Treeless & Barren Lands","Ice, Uninhabited")
write.table(Anthromes_vis, "S:/Analysis/Anthromes/Anthromes_vis.csv")

# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)

Representativity_Anthromes <-
  ggplot(Anthromes_vis, aes(x = as.factor(Anthromes_vis$value))) + 
  geom_bar(aes(y = Anthromes_vis$count_global/sum(Anthromes_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = Anthromes_vis$count_ilter/sum(Anthromes_vis$count_ilter)*100, color = Anthromes_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Anthrome Category based on Ellis et al. (2020)", y = "Relative Spatial Coverage [%]", 
       color = "Geographic\nRepresentedness") + 
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1), 
                        labels = c("Underrepresented", "Well represented", "Overrepresented"))+
  scale_x_discrete(labels = Anthromes_vis$category_name) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=13, xmin = 15) + ylim(0,15)
ggsave(filename = "Representativity_Anthromes_1000.jpg", plot = Representativity_Anthromes,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geographic Representedness
Rep_Anthromes <- reclassify(Anthromes, Anthromes_vis[,c(1,4)], filename = "S:/Results/rep_maps/Rep_Anthromes.tif", overwrite = T) 
plot(Rep_Anthromes) 