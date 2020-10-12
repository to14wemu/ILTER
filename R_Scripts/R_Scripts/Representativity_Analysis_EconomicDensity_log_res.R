## Analyzing Geographic Representedness of the ILTER network for the Economic Density dataset ##
rm(list = ls())
write("TMPDIR = 'S:/temp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

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

# 2. Load Data
# 2.1 Accredited ILTER sites
ILTER_acc_poly  <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                           layer = "ilter_boundaries")
ILTER_point_buff<- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
                           layer = "Points_with_buffer")

# Point site is only one marine site, excluded from the analysis
#ILTER_points    <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
#                          layer = "Points_without_size")

# 2.2 Economic Density
EconomicDensity                        <- raster("S:/GIS_data_processed/Economic_Density/EconomicDensity_MW_res_30M.tif")
ILTER_acc_poly                         <- spTransform(ILTER_acc_poly, CRSobj = crs(EconomicDensity))
#ILTER_points                           <- spTransform(ILTER_points, CRSobj = crs(EconomicDensity))
ILTER_point_buff                       <- spTransform(ILTER_point_buff, CRSobj = crs(EconomicDensity))
plot(EconomicDensity);plot(ILTER_acc_poly, add = T);plot(ILTER_point_buff, add = T)#;plot(ILTER_points, add = T)

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
ED_Global         <- freq(EconomicDensity)
ED_Global
ED_Global         <- data.frame(ED_Global)
ED_Global         <- ED_Global[1:(length(ED_Global$value)-1),]  #remove NA
ED_Global$percent <- ED_Global$count*100/sum(ED_Global$count)
ED_Global
write.csv(ED_Global, "S:/Analysis/EconomicDensity/ED_Global_freq.csv")

#     3.1.2    Weighted Distribution on ILTER sites
## Polygon sites ##
ED_ILTER                        <- raster::extract(EconomicDensity, ILTER_acc_poly, small = T, df = T) 
count                           <- count(ED_ILTER, ID) # n is the amount of raster pixels/values in a polygon
ILTER_acc_poly$ID               <- as.numeric(rownames(ILTER_acc_poly@data))+1 # +1 to match IDs in count and anthromes ILTER
ED_ILTER                        <- full_join(ED_ILTER, count, "ID")
ED_ILTER                        <- full_join(ED_ILTER, ILTER_acc_poly@data, "ID")
ED_ILTER$IntendedPercent_Cell   <- ED_ILTER$IntendedPercent_Poly/ED_ILTER$n 
sum(ED_ILTER$IntendedPercent_Cell, na.rm = T) # check: value ~ 70 - 71
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be 
# applied for intendedPercent to be reached for representativity analysis

ILTER_poly <- ED_ILTER %>% group_by(EconomicDensity_MW_res_30M) %>% summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight  ####sum(n, na.rm = T),
count_cat  <- count(ED_ILTER, EconomicDensity_MW_res_30M)
ILTER_poly <- left_join(ILTER_poly, count_cat, "EconomicDensity_MW_res_30M")

colnames(ILTER_poly) <- c("ED", "IntendedPercent_Cell", "actual_n")
ILTER_poly
ILTER_poly$desired_n <- sum(ILTER_poly$actual_n) * ILTER_poly$IntendedPercent_Cell/100
ILTER_poly           <- ILTER_poly[1:(length(ILTER_poly$desired_n)-1),] #remove NA

## Point sites ##
ED_ILTER_point                        <- raster::extract(EconomicDensity, ILTER_point_buff, small = T, df = T) 
count                                 <- count(ED_ILTER_point, ID) # n is the amount of raster pixels/values in a polygon
ILTER_point_buff$ID                   <- as.numeric(rownames(ILTER_point_buff@data))+1 # +1 to match IDs in count and anthromes ILTER
ED_ILTER_point                        <- full_join(ED_ILTER_point, count, "ID")
ED_ILTER_point                        <- full_join(ED_ILTER_point, ILTER_point_buff@data, "ID")
ED_ILTER_point$IntendedPercent_Cell   <- ED_ILTER_point$IntendedPercent_Poly/ED_ILTER_point$n 
sum(ED_ILTER_point$IntendedPercent_Cell, na.rm = T) # build in check: should be ~ 29-30
# knowing both the count (n) and the intendedPercent allows for calculation of a weight that needs to be applied for intendedPercent
# to be reached for representativity analysis

ILTER_point               <- ED_ILTER_point %>% group_by(EconomicDensity_MW_res_30M) %>% 
  summarize(sum(IntendedPercent_Cell, na.rm = T))
# counting all raster cells on the sites and assigning the mean weight 
count_cat                 <- count(ED_ILTER_point, EconomicDensity_MW_res_30M)
ILTER_point               <- left_join(ILTER_point, count_cat, "EconomicDensity_MW_res_30M")

colnames(ILTER_point)     <- c("ED", "IntendedPercent_Cell", "actual_n")
ILTER_point
ILTER_point$desired_n     <- sum(ILTER_point$actual_n) * ILTER_point$IntendedPercent_Cell/100
ILTER_point               <- ILTER_point[1:(length(ILTER_point$desired_n)-1),] #remove NA

## Combining Polygon and Point Site Data ##
ILTER                           <- full_join(ILTER_poly, ILTER_point, "ED")
ILTER$desired_n.y[is.na(ILTER$desired_n.y)] <- 0
ILTER$vis_n                     <- ILTER$desired_n.x + ILTER$desired_n.y

write.csv(ILTER, "S:/Analysis/EconomicDensity/ED_ILTER_categorized.csv")       
write.csv(ED_ILTER, "S:/Analysis/EconomicDensity/ED_ILTER.csv")               

#####
#   3.2 Representativity Analysis: chisq.test
#         A matrix with both the global and ILTER distribution is required
# Tidy tables
df_ED           <- full_join(ED_Global[,1:2], ILTER[,c(1,8)], by = c( "value" = "ED"))
colnames(df_ED) <- c("value", "count_global", "count_ilter")

# Produce matrix for chisq.test
ED_vis                     <- df_ED                         # Necessary for subsequent visualization
ED_chi                     <- t(df_ED)                     # t() also produces matrix
colnames(ED_chi)           <- ED_chi[1,]
ED_chi                     <- ED_chi[-1,]
ED_chi[is.na(ED_chi)]      <- 0
chisq                      <- chisq.test(ED_chi)                  
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
ED_vis <- data.frame(ED_vis, rep)
ED_vis$rep[ED_vis$count.ILTER == 0] <- NA
ED_vis$category_name <- c("<0.01","0.01-0.1","0.1-1","1-10","10-30",">30")
write.table(ED_vis, "S:/Analysis/EconomicDensity/ED_vis.csv")

# produce table to add to ggplot layout
chisq_vis           <- data.frame(c(chisq$p.value, chisq$statistic))
colnames(chisq_vis) <- "Representativity"
rownames(chisq_vis) <- c("p value","x squared")
chisq_vis           <- round(chisq_vis, digits = 0)
mytheme             <- ttheme_default(base_size = 10)


Representativity_ED<-
  ggplot(ED_vis, aes(x = as.factor(ED_vis$value))) + 
  geom_bar(aes(y = ED_vis$count_global/sum(ED_vis$count_global)*100), stat = "identity") + 
  geom_point(aes(y = ED_vis$count_ilter/sum(ED_vis$count_ilter)*100, color = ED_vis$rep), 
             alpha = 0.8, size = 3.5) + theme_bw() +
  labs(x = "Economic Density [Mio. US$/km²]", y = "Relative Spatial Coverage [%]", color = "Geographic\nRepresentedness") +
  scale_color_gradient2(midpoint = 0, high = rgb(30,158,75, max = 255), mid = rgb(255,255,191, max = 255),
                        low = rgb(158,30,133, max = 255), space = "Lab", breaks=c(-1,0,1),
                        labels = c("Underrepresented", "Well represented", "Overrepresented")) +
  scale_x_discrete(labels = ED_vis$category_name) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  annotation_custom(tableGrob(chisq_vis, theme = mytheme), ymin=29, xmin = 5)

ggsave(filename = "Representativity_ED_1000.jpg", plot = Representativity_ED,
       dpi = 1000, path = "S:/Results/rep_histograms", width = 28, height = 20, units = "cm", device = "jpg")

## Mapping Geographic Representedness
Rep_ED <- reclassify(EconomicDensity, ED_vis[,c(1,4)], filename = "S:/Results/rep_maps/Rep_ED.tif", overwrite = T) 
plot(Rep_ED) 
