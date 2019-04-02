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

dir = paste0("/Users/", Sys.info()[7], "/")
load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
depth = myExpl[[3]]
plot(depth)
proj4string(depth) = CRS("+init=epsg:4326")
slope_asp = terrain(depth, opt=c('slope', 'aspect'), unit='degrees', neighbors=8)
slope = slope_asp[[1]]
plot(slope_asp[[1]])

sp = c("lobster", "scallop")[1]
setwd(paste0(dir, "biomod_migclim/", sp))
load("data_and_variables_for_response_curves.RData")
depth = Data[, c(5,4,3)]
qplot(depth$var1.pred.5, depth$var1.pred.4, color = depth$var1.pred.3)
colnames(depth) = c("Lon","Lat","Depth")

coordinates(depth) = ~ Lon + Lat
auto = autofitVariogram(depth$Depth ~ 1, depth)
depth <- depth[-zerodist(depth)[,1],] 
g = gstat(formula = Depth ~ 1, model = auto$var_model, data = depth, maxdist = 0.08) #kisei is using 0.08
xrange = range(depth$Lon); yrange = range(depth$Lat)
grid= expand.grid(Lon = seq(from = xrange[1], to = xrange[2], by = .05), #0.01 
                  Lat = seq(from = yrange[1], to = yrange[2], by = .05)) #0.01
gridded(grid) = ~Lon + Lat
p = predict(g, newdata = grid)
depth = raster(p)

plot(depth)
proj4string(depth) = CRS("+init=epsg:4326")
slope_asp = terrain(depth, opt=c('slope', 'aspect'), unit='degrees', neighbors=8)
slope = slope_asp[[1]]
plot(slope_asp[[1]])

depth_slope = stack(depth, slope)

df = rasterToPoints(depth_slope)
df = as.data.frame(df)
df = df[,c("x", "y", "slope")]
qplot(df$x, df$y, colour = df$slope )

ggplot(subset(df, slope < 1.1))+
  geom_point(aes(x, y, color = slope), alpha = 0.5, size = 1.5) + 
  scale_x_continuous("") + 
  scale_y_continuous("") + 
  scale_colour_viridis_c("")
  
qplot(Data$var1.pred.5, Data$var1.pred.4, colour = Data$var1.pred.3)

colnames(df) = c("var1.pred.5", "var1.pred.4", "var1.pred.6")

d = merge(Data, df, all = T)
# d = d[!duplicated(d$var1.pred.1), ]

Data = d

save(Data, file = "data_and_variables_for_response_curves_with_slope.RData")
