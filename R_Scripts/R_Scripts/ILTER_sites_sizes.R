## Analysis of the effect of different weighting methods on the influence of single ILTER sites ##
rm(list = ls())

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
# Accreted ILTER sites
ILTER_acc_poly  <- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"), 
                           layer = "ilter_boundaries")
ILTER_point_buff<- readOGR(dsn=path.expand("S:/GIS_data_original/ILTER/2020_09_28"),
                           layer = "Points_with_buffer")
mollproj               <- '+proj=moll +ellps=WGS84'
ILTER_acc_poly         <- spTransform(ILTER_acc_poly, CRSobj = mollproj)
ILTER_acc_point        <- spTransform(ILTER_point_buff, CRSobj = mollproj)
plot(ILTER_acc_poly);plot(ILTER_acc_point, add = T)

# 3. Analysis
# 3.1 Original areal extents
ILTER_acc_poly$area  <- raster::area(ILTER_acc_poly) / 1000000
ILTER_acc_point$area <- raster::area(ILTER_acc_point)/ 1000000

area <- c(ILTER_acc_poly$area, ILTER_acc_point$area)
area <- sort(area, decreasing = F)
min(raster::area(ILTER_acc_poly))
max(area)
ID   <- 1:(length(area))
area <- data.frame(ID,area)

area$percent    <- area$area/sum(area$area)*100
area$cumPercent <- rep(area$percent[1], length(area$area))
for (i in 2:length(area$area)){
  area$cumPercent[i] <- area$percent[i] + area$cumPercent[i - 1]
}
ggplot(area, aes(x = ID, y = cumPercent)) + geom_line(color = "steelblue", size = 1.2) + theme_bw() + labs(y = "Cumulative Percent")

# 3.2 LOG
area_samp6 <- log(area$area)
area_samp6 <- area_samp6 + abs(min(area_samp6)) + 1
area_samp6 <- sort(area_samp6, decreasing = F)
area_samp6 <- data.frame(ID,area_samp6)

area_samp6$percent     <- area_samp6$area_samp6/sum(area_samp6$area_samp6)*100
area_samp6$cumPercent  <- rep(area_samp6$percent[1], length(area$area))
for (i in 2:length(area$area)){
  area_samp6$cumPercent[i] <- area_samp6$percent[i] + area_samp6$cumPercent[i - 1]
}
ggplot(area_samp6, aes(x = ID, y = cumPercent)) + geom_line(color = "steelblue", size = 1.2) + theme_bw() + 
  labs(y = "Cumulative Percent", x = "ID - Sorted by Area Ascending")

# 3.3 Same Weight
Same_Weight    <- c()
Same_Weight$ID <- seq(1, 742, 1)
Same_Weight$percent <- rep(1/742*100, 742)
Same_Weight <- data.frame(Same_Weight)
Same_Weight$cumPercent  <- rep(Same_Weight$percent[1], length(area$area))
for (i in 2:length(area$area)){
  Same_Weight$cumPercent[i] <- Same_Weight$percent[i] + Same_Weight$cumPercent[i - 1]
}

# 3.4 combine to one plot
## Comparison of the different weighting methods 
areal_plots <-
ggplot(data = area, aes(x = ID, y = cumPercent)) + 
  geom_line(aes(color = "No Weighting"), size = 1.2) +
  # geom_line(data = area_samp, aes(color = "10,000 km²"), size = 1.2) +
  # geom_line(data = area_samp2, aes(color = "1,000 km²"), size = 1.2) + 
  #geom_line(data = area_samp3, aes(color = "100 km²"), size = 1.2) +
  # geom_line(data = area_samp4, aes(color = "10 km²"), size = 1.2) +
  #geom_line(data = area_samp5, aes(color = "1 km²"), size = 1.2) +
  geom_line(data = area_samp6, aes(color = "Log-Weighted"), size = 1.2) +
  geom_line(data = Same_Weight, aes(color = "Same Weight"), size = 1.2)+
  #geom_abline(aes(color = "Same Weight", intercept = 0, slope = 100/length(ILTER_acc_poly)), size = 1.2) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 360, vjust = 0.5),
                     legend.position = c(0.2,0.8)) +
  scale_color_manual("Weighting",
                    # breaks = c("1:1 Line", "Original",   "Log" ),
                     values = c("Same Weight" = "grey", "No Weighting" = "coral", "Log-Weighted" = "steelblue")) +
  labs(y = "Cumulative \nPercentage \nof total \nArea within \nILTER",
       x = "ILTER sites \nSorted by Area Ascending")
ggsave(areal_plots, path = "S:/Results", dpi = 600, height = 10, width = 15, units = "cm", device = "png",
       filename = "areal_plots.png")
ggsave(areal_plots, path = "S:/Results", dpi = 1000, height = 10, width = 15, units = "cm", device = "jpg",
       filename = "areal_plots_1000.jpg")

## Comparison of single polygon weights in log distribution
log_plot <-
  ggplot(area_samp6, aes(x = ID, y = percent)) + geom_line(color = "steelblue", size = 1.2) + theme_bw() + 
  labs(y = "Percent", x = "ID - Sorted by Area Ascending")
ggsave(log_plot, path = "S:/Results", dpi = 600, height = 10, width = 15, units = "cm", device = "png",
       filename = "log_plot.png")
