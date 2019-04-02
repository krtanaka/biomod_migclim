# load cm2.6 data and project habitat change year 1-80

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

sp = c("lobster", "scallop")[1]
# op = c("tuned_lobster", "tuned_scallop", "default")[3]

setwd(paste0(dir, "/Google Drive/R/biomod/", sp))
setwd(paste0(dir, "/Desktop/", sp))
load("Biomod_Results.RData")
load("final_sdms.RData")

dir = paste0("/Users/", Sys.info()[7], "/")

# static = read_csv(paste0(dir, "Google Drive/R/Biomod/data/CM_0.05.csv"), col_names  = T)
static = read_csv(paste0(dir, "Google Drive/R/Biomod/data/CM_0.05_with_slope.csv"), col_names  = T)
load(paste0(dir, "Google Drive/R/Biomod/data/distance_offshore_for_static_data_low_res.rda"))

static$Distant_Offshore = df$distance

month = c("apr", "may", "jun", "sep", "oct", "nov")
season = c("spring", "fall")[1]

period = c(9:18) #first 10 years
period = c(79:88) #last 10 years

for (i in 1:length(month)) {
  
  bt = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_", month[i], ".csv")), col_names  = T)
  bs = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_", month[i], ".csv")), col_names  = T)
  
  bt = merge(static, bt, by = c('Lat', 'Lon'))
  bs = merge(static, bs, by = c('Lat', 'Lon'))
  
  bt$Depth = bt$Depth*-1
  bs$Depth = bs$Depth*-1
  bt = subset(bt, Depth>0)
  bs = subset(bs, Depth>0)
  bt = subset(bt, Depth<400)
  bs = subset(bs, Depth<400)
  
  btmed = cbind(bt[,c(1,2,7)], apply(bt[,period],1, median)); names(btmed)[4] = 'temp'
  bsmed = cbind(bs[,c(1,2,7)], apply(bs[,period],1, median)); names(bsmed)[4] = 'sal'
  
  depth = btmed[,c(1,2,3)]
  btmed = btmed[,c(1,2,4)]
  bsmed = bsmed[,c(1,2,4)]
  lat = btmed[,c(1,2,1)]
  lon = bsmed[,c(1,2,2)]
  
  colnames(depth)[3] = "var1.pred"
  colnames(btmed)[3] = "var1.pred"
  colnames(bsmed)[3] = "var1.pred"
  colnames(lat)[3] = "var1.pred"
  colnames(lon)[3] = "var1.pred"
  
  bt <- rasterFromXYZ(btmed[, c("Lon", "Lat", "var1.pred")])
  bs <- rasterFromXYZ(bsmed[, c("Lon", "Lat", "var1.pred")])
  depth <- rasterFromXYZ(depth[, c("Lon", "Lat", "var1.pred")])
  lat <- rasterFromXYZ(lat[, c("Lon", "Lat", "var1.pred")])
  lon <- rasterFromXYZ(lon[, c("Lon", "Lat", "var1.pred")])
  myExplFuture = stack(bt,bs,depth,lat,lon)
  # rm(bs, bsmed, bt, btmed, depth, df, lat, lon, static)
  
  myBiomodProjFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExplFuture,
                                          proj.name = 'future',
                                          selected.models = final_sdms,
                                          binary.meth = 'TSS',
                                          compress = 'xz',
                                          clamping.mask = T,
                                          output.format = '.grd')
  
  # projection with final SDMs, weighted mean
  d = myBiomodProjFuture@proj@val@layers
  dd = stack(d)
  ddd = subset(dd,final_sdms)
  
  avg = weighted.mean(ddd, mw$weight)
  
  if(period %in% c(9:18)){
    
    jpeg(filename = paste0(dir, "Google Drive/R/Biomod/lobster/future_", month[i], "_1_10.jpg"), 
         res = 500, height = 7, width = 7, units = "in")
    plot(avg, zlim = c(0,1000), col = parula(100))
    map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
    dev.off()
    
  } else { 
    
    jpeg(filename = paste0(dir, "Google Drive/R/Biomod/lobster/future_", month[i], "_70_80.jpg"), 
         res = 500, height = 7, width = 7, units = "in")
    plot(avg, zlim = c(0,1000), col = parula(100))
    map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
    dev.off()
    
  }
  
} #running biomod by monthly time step
for (i in 1:length(season)) {
  
  if (season[i] == "spring") {
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_apr.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_may.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_jun.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]; rm(bt_1, bt_2, bt_3)
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_apr.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_may.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_jun.csv")), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  if (season[i] == "fall") {
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_sep.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_oct.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_nov.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_sep.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_oct.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_nov.csv")), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  bt = merge(static, bt, by = c('Lat', 'Lon'))
  bs = merge(static, bs, by = c('Lat', 'Lon'))
  
  bt$Depth = bt$Depth*-1
  bs$Depth = bs$Depth*-1
  bt = subset(bt, Depth>0)
  bs = subset(bs, Depth>0)
  bt = subset(bt, Depth<400)
  bs = subset(bs, Depth<400)
  
  btmed = cbind(bt[,c(1,2,7)], apply(bt[,period],1, median)); names(btmed)[4] = 'temp'
  bsmed = cbind(bs[,c(1,2,7)], apply(bs[,period],1, median)); names(bsmed)[4] = 'sal'
  
  depth = btmed[,c(1,2,3)]
  btmed = btmed[,c(1,2,4)]
  bsmed = bsmed[,c(1,2,4)]
  lat = btmed[,c(1,2,1)]
  lon = bsmed[,c(1,2,2)]
  
  colnames(depth)[3] = "var1.pred"
  colnames(btmed)[3] = "var1.pred"
  colnames(bsmed)[3] = "var1.pred"
  colnames(lat)[3] = "var1.pred"
  colnames(lon)[3] = "var1.pred"
  
  bt <- rasterFromXYZ(btmed[, c("Lon", "Lat", "var1.pred")])
  bs <- rasterFromXYZ(bsmed[, c("Lon", "Lat", "var1.pred")])
  depth <- rasterFromXYZ(depth[, c("Lon", "Lat", "var1.pred")])
  lat <- rasterFromXYZ(lat[, c("Lon", "Lat", "var1.pred")])
  lon <- rasterFromXYZ(lon[, c("Lon", "Lat", "var1.pred")])
  myExplFuture = stack(bt,bs,depth,lat,lon)
  # rm(bs, bsmed, bt, btmed, depth, df, lat, lon, static)
  
  myBiomodProjFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExplFuture,
                                          proj.name = 'future',
                                          selected.models = final_sdms,
                                          binary.meth = 'TSS',
                                          compress = 'xz',
                                          clamping.mask = T,
                                          output.format = '.grd')
  
  # projection with final SDMs, weighted mean
  d = myBiomodProjFuture@proj@val@layers
  dd = stack(d)
  ddd = subset(dd,final_sdms)
  
  avg = weighted.mean(ddd, mw$weight)
  
  if(period %in% c(9:18)){
    
    jpeg(filename = paste0(dir, "Google Drive/R/Biomod/lobster/future_", season[i], "_1_10.jpg"), 
         res = 500, height = 7, width = 7, units = "in")
    plot(avg, zlim = c(0,1000), col = parula(100))
    map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
    dev.off()
    
  } else { 
    
    jpeg(filename = paste0(dir, "Google Drive/R/Biomod/lobster/future_", season[i], "_70_80.jpg"), 
         res = 500, height = 7, width = 7, units = "in")
    plot(avg, zlim = c(0,1000), col = parula(100))
    map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
    dev.off()
    
  }
  
} #running biomod by seasonal time step

# Loop 80 years, create dynamic environment dataset for migclim analysis with CM.2.6-------------------------------------

setwd(paste0(dir, "/Google Drive/R/Biomod/", sp))load("Biomod_Results.RData"); load("final_sdms.RData")
setwd(paste0(dir, "/Desktop/")); load(paste0(dir, "/Desktop/", sp,"/Biomod_Results.RData")); load(paste0(dir, "/Desktop/", sp,"/final_sdms.RData"))

season = c("spring", "fall", "annual")[1]
for (k in 1:length(season)){
  
  season = season[k]
  
  if (season == "spring") {
    
    bt_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_apr.csv"), col_names  = T)
    bt_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_may.csv"), col_names  = T)
    bt_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_jun.csv"), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]; rm(bt_1, bt_2, bt_3)
    
    bs_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_apr.csv"), col_names  = T)
    bs_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_may.csv"), col_names  = T)
    bs_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_jun.csv"), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]; rm(bs_1, bs_2, bs_3)
    
  }
  
  if (season == "fall") {
    
    bt_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_sep.csv"), col_names  = T)
    bt_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_oct.csv"), col_names  = T)
    bt_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_nov.csv"), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]; rm(bt_1, bt_2, bt_3)
    
    bs_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_sep.csv"), col_names  = T)
    bs_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_oct.csv"), col_names  = T)
    bs_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_nov.csv"), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]; rm(bs_1, bs_2, bs_3)
    
  }
  
  if (season == "annual") {
    
    bt = read_csv("/Users/kisei/Google Drive/R/Biomod/data/temp/Delta_Annual_temp.csv", col_names  = T)
    bs = read_csv("/Users/kisei/Google Drive/R/Biomod/data/salt/Delta_Annual_sal.csv", col_names  = T)
    
  }
  
  static = read_csv("/Users/ktanaka/Google Drive/R/Biomod/data/CM_0.05.csv", col_names  = T)
  static = read_csv("/Users/ktanaka/Google Drive/R/Biomod/data/CM_0.05_with_slope.csv", col_names  = T)
  
  load("/Users/kisei/Google Drive/R/Biomod/data/distance_offshore_for_static_data_low_res.rda")
  
  static$Distant_Offshore = df$distance
  
  bt = merge(static, bt, by = c('Lat', 'Lon'))
  bs = merge(static, bs, by = c('Lat', 'Lon'))
  
  bt$Depth = bt$Depth*-1
  bs$Depth = bs$Depth*-1
  bt = subset(bt, Depth>0)
  bs = subset(bs, Depth>0)
  bt = subset(bt, Depth<400)
  bs = subset(bs, Depth<400)
  
  for (i in 1:80){
    
    depth = bt[,c(1,2,7)];colnames(depth)[3] = "var1.pred"
    btmed = bt[,c(1,2,(i+8))];colnames(btmed)[3] = "var1.pred"
    bsmed = bs[,c(1,2,(i+8))];colnames(bsmed)[3] = "var1.pred"
    lat = bt[,c(1,2,1)];colnames(lat)[3] = "var1.pred"
    lon = bt[,c(1,2,2)];colnames(lon)[3] = "var1.pred"
    
    #coordinates(bt) = ~Lon + Lat; 
    #rasterize(bt$Lon, bt$Lat, bt$LMZ)
    btemp <- rasterFromXYZ(btmed[, c("Lon", "Lat", "var1.pred")])
    bsalt <- rasterFromXYZ(bsmed[, c("Lon", "Lat", "var1.pred")])
    depth <- rasterFromXYZ(depth[, c("Lon", "Lat", "var1.pred")])
    lat <- rasterFromXYZ(lat[, c("Lon", "Lat", "var1.pred")])
    lon <- rasterFromXYZ(lon[, c("Lon", "Lat", "var1.pred")])
    
    myExplFuture = stack(btemp,bsalt,depth,lat,lon); rm(btemp,bsalt,depth,lat,lon)
    
    myBiomodProjFuture <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = myExplFuture,
      proj.name = 'future',
      selected.models = final_sdms,
      binary.meth = 'TSS',
      compress = 'xz',
      clamping.mask = T,
      output.format = '.grd')
    
    d = myBiomodProjFuture@proj@val@layers
    dd = stack(d)
    ddd = subset(dd,final_sdms)
    
    avg = weighted.mean(ddd, mw$weight)
    
    spts <- rasterToPoints(avg, spatial = TRUE)
    x <- as.data.frame(spts)
    x$layer = as.integer(x$layer)
    if (i==1){
      POS = x[,c(2,3,1)]
      colnames(POS)[3] = "hsmap1"
    }
    
    if (i>1){
      POS[,(i+2)] = x[,1]
      colnames(POS)[i+2] = paste("hsmap", i, sep="")
    }
    
  }
  
  # depth <- rasterToPoints(depth, spatial = TRUE)
  # depth <- as.data.frame(depth); colnames(depth)[1] = "depth"
  # 
  # btemp <- rasterToPoints(btemp, spatial = TRUE)
  # btemp <- as.data.frame(btemp);colnames(btemp)[1] = "btemp"
  # 
  # bsalt <- rasterToPoints(bsalt, spatial = TRUE)
  # bsalt <- as.data.frame(bsalt);colnames(bsalt)[1] = "bsalt"
  # 
  # POS1 = merge(POS, depth)
  # POS1 = merge(POS1, btemp)
  # POS1 = merge(POS1, bsalt)
  
  # POS = POS1[,c(1:2,83:85, 3:82)]
  
  if (op == "default") {
    
    write_csv(POS, paste0("Biomod_1_80_", season, "_with_default_options.csv"), col_names = T)
    
  }else{
    
    write_csv(POS, paste0("Biomod_1_80_", season, ".csv"), col_names = T)
  }
  
} #run by seasonal timesteps

