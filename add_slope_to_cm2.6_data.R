rm(list = ls())

library(biomod2)
library(raster)
library(readr)
library(colorRamps)
library(automap)
library(beepr)
library(gstat)
library(cowplot)
library(maps)
library(Rmisc)
library(svMisc)
library(pals)
library(doParallel)
library(data.table)
library(grid)
library(ENMeval)
library(ggpubr)
library(mgcv)

dir = paste0("/Users/", Sys.info()[7], "/")

load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
load(paste0(dir, "biomod_migclim/data/slope.RData"))
myExpl = stack(myExpl, slope)
names(myExpl)[6] = "var1.pred.6"

d = as.data.frame(rasterToPoints(myExpl))
d = d[,c("var1.pred.3","var1.pred.5", "var1.pred.4", "var1.pred.6")]
colnames(d) = c("Depth", "Lon", "Lat", "Slope")
qplot(d$Lon, d$Lat, colour = d$Slope) + 
  
summary(d)

d = d[which(complete.cases(d$Slope)),]

g = gam(Slope ~ s(Lon, Lat) + s(Depth), data = d)

setwd(paste0(dir, "Google Drive/R/Biomod/data/"))
static = read_csv("CM_0.05.csv", col_names  = T)
names(static)

static$Depth = static$Depth*-1

static$Slope = predict(g, static)

qplot(static$Lon, static$Lat, color = static$Slope)

ggplot(subset(static, Slope > -1))+
  geom_point(aes(Lon, Lat, color = Slope))

write.csv(static, file = "CM_0.05_with_slope.csv")
