#running biomod with best sdms (tss > 0.5, ROC > 0.8) as well as best and worst sdm

#best losbter sdm = GBM Run 3
#worlst lobster sdm = ANN run 1

#best scallop sdm = GBM run 3
#worst scallop sdm = Maxent 1 run 1

library(readr)
library(data.table)
library(raster)
library(biomod2)

sp = c("lobster", "scallop")[2]

dir = "/Users/Kisei/"

op = c("tuned_lobster", "tuned_scallop", "default")[2]

if (sp == "lobster") {
  
  setwd(paste0(dir, "/Desktop"))
  load(paste0(dir, "/Desktop/lobster/Biomod_Results.RData"))
  load(paste0(dir, "/Desktop/lobster/final_sdms.RData"))
  
} 

if (sp == "scallop") {
  
  setwd(paste0(dir, "/Desktop"))
  load(paste0(dir, "/Desktop/scallop/Biomod_Results.RData"))
  load(paste0(dir, "/Desktop/scallop/final_sdms.RData"))
  
}

# load(paste0(dir, "/Google Drive/R/Biomod/lobster/Biomod_Results_with_Default_Options.RData"))

season = c("spring", "fall", "annual")[3]

for (s in 1:length(season)){
  
  season = season[s]
  
  if (season == "spring") {
    
    bt_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_apr.csv"), col_names  = T)
    bt_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_may.csv"), col_names  = T)
    bt_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_jun.csv"), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_apr.csv"), col_names  = T)
    bs_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_may.csv"), col_names  = T)
    bs_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_jun.csv"), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  if (season == "fall") {
    
    bt_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_sep.csv"), col_names  = T)
    bt_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_oct.csv"), col_names  = T)
    bt_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/temp/Delta_nov.csv"), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_sep.csv"), col_names  = T)
    bs_2 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_oct.csv"), col_names  = T)
    bs_3 = read_csv(paste0(dir, "Google Drive/R/Biomod/data/salt/Delta_nov.csv"), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  if (season == "annual") {
    
    bt = read_csv("/Users/kisei/Google Drive/R/Biomod/data/temp/Delta_Annual_temp.csv", col_names  = T)
    bs = read_csv("/Users/kisei/Google Drive/R/Biomod/data/salt/Delta_Annual_sal.csv", col_names  = T)
    
  }
  
  static = read_csv("/Users/Kisei/Google Drive/R/Biomod/data/CM_0.05.csv", col_names  = T)
  load("/Users/Kisei/Google Drive/R/Biomod/data/distance_offshore_for_static_data_low_res.rda")
  
  static$Distant_Offshore = df$distance
  
  bt = merge(static, bt, by = c('Lat', 'Lon'))
  bs = merge(static, bs, by = c('Lat', 'Lon'))
  
  bt$Depth = bt$Depth*-1
  bs$Depth = bs$Depth*-1
  bt = subset(bt, Depth>0)
  bs = subset(bs, Depth>0)
  bt = subset(bt, Depth<400)
  bs = subset(bs, Depth<400)
  
  for (m in 1:length(final_sdms)) {
    
    # m = 1
    
    for (i in 1:80){
      
      # i = 1
      
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
      
      myExplFuture = stack(btemp,bsalt,depth,lat,lon)
      
      myBiomodProjFuture <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = myExplFuture,
        proj.name = 'future',
        selected.models = final_sdms[m],
        binary.meth = 'TSS',
        compress = 'xz',
        clamping.mask = T,
        output.format = '.grd')
      
      d = myBiomodProjFuture@proj@val@layers
      dd = stack(d)
      ddd = subset(dd,final_sdms)
      spts = rasterToPoints(ddd, spatial = T)
      x <- as.data.frame(spts)
      
      if (i==1){
        
        POS = x[,c(2,3,1)]
        colnames(POS)[3] = "hsmap1"
        
      } else {
        
        POS[,(i+2)] = x[,1]
        colnames(POS)[i+2] = paste("hsmap", i, sep="")
        
      }
    }
    
    POS$Model = final_sdms[m]
    
    POS = POS[,c(1:2, 83, 3:82)]
    
    if (op == "default") {
      
      write_csv(POS, paste0("Biomod_1_80_", final_sdms[m], "_", season, "_with_default_options.csv"), col_names = T)
      
    }else{
      
      write_csv(POS, paste0("Biomod_1_80_", final_sdms[m], "_", season, ".csv"), col_names = T)
    }
    
  }
  
}
