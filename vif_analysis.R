dir = paste0("/Users/", Sys.info()[7], "/")

load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
load(paste0(dir, "biomod_migclim/data/slope.RData"))
myExpl = stack(myExpl, slope)
names(myExpl)[6] = "var1.pred.6"

v1 = rasterToPoints(myExpl$var1.pred.1)
v2 = rasterToPoints(myExpl$var1.pred.2)
v3 = rasterToPoints(myExpl$var1.pred.3)
v4 = rasterToPoints(myExpl$var1.pred.4)
v5 = rasterToPoints(myExpl$var1.pred.5)
v6 = rasterToPoints(myExpl$var1.pred.6)

v = merge(v1, v2)
v = merge(v, v3)
v = merge(v, v4)
v = merge(v, v5)
v = merge(v, v6)

colnames(v) = c("lat", "lon", "temp", "salt", "depth", "lat", "lon", "slope")

source("/Users/ktanaka/Google Drive/R/misc/vif.R")

vif <- v[, c("temp", "salt","depth", "slope")]
vif_func(in_frame = vif, thresh=3, trace=T) 
