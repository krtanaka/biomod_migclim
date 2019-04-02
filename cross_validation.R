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
op = c("tuned", "default")[1]

if (sp == "lobster") {
  
  load(paste0(dir, "biomod_migclim/lobster/lobster_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data
  
  DataSpecies = lobster[,c(7,8,24)] #Lat, Lon, Catch Number
  DataSpecies$expcatnum = ifelse(DataSpecies$expcatnum > 1, 1, 0) #try chaniging presence/absence threashold to make it less sensitive
  names(DataSpecies)[1] = 'Y_WGS84'
  names(DataSpecies)[2] = 'X_WGS84'
  names(DataSpecies)[3] = 'lobster'
  
  myRespName = 'lobster' # the name of studied species
  myResp = as.numeric(DataSpecies[,myRespName]) # the presence/absences data for our species
  myRespXY = DataSpecies[,c("X_WGS84","Y_WGS84")] # the XY coordinates of species data
  
  # Load raster stacked explanatory variables
  # Prepared by kriging interpolation, maxdist = 0.08, resolution = 0.05
  load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
  load(paste0(dir, "biomod_migclim/data/slope.RData"))
  myExpl = stack(myExpl, slope)
  names(myExpl)[6] = "var1.pred.6"
  
  # par(mfrow = c(1,3))
  # plot(myExpl$var1.pred.1, main = "Bottom Temp (deg C)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  # plot(myExpl$var1.pred.2, main = "Bottom Salt (ppt)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  # plot(myExpl$var1.pred.3, main = "Depth (m)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName)
  
  rm(DataSpecies, lobster); myBiomodData; plot(myBiomodData)
  
}

if (sp == "scallop") {
  
  load(paste0(dir, "biomod_migclim/scallop/scallop_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data
  scallop$year = substr(scallop$cruise6, 1, 4); scallop = subset(scallop, year %in% c(1984:2016))
  
  DataSpecies = scallop[,c(19, 20, 16)] #Lat, Lon, Catch Number
  DataSpecies$expcatnum = ifelse(DataSpecies$expcatnum > 1, 1, 0) #try chaniging presence/absence threashold to make it less sensitive
  names(DataSpecies)[1] = 'Y_WGS84'
  names(DataSpecies)[2] = 'X_WGS84'
  names(DataSpecies)[3] = 'scallop'
  
  myRespName = 'scallop' # the name of studied species
  myResp = as.numeric(DataSpecies$scallop) # the presence/absences data for our species
  myRespXY = DataSpecies[,c("X_WGS84","Y_WGS84")] # the XY coordinates of species data
  
  # Load raster stacked explanatory variables
  # Prepared by kriging interpolation, maxdist = 0.08, resolution = 0.05
  load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
  # load(paste0(dir, "biomod_migclim/data/slope.RData"))
  # myExpl = stack(myExpl, slope)
  # names(myExpl)[6] = "var1.pred.6"
  
  # par(mfrow = c(1,3))
  # plot(myExpl$var1.pred.1, main = "Bottom Temp (deg C)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  # plot(myExpl$var1.pred.2, main = "Bottom Salt (ppt)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  # plot(myExpl$var1.pred.3, main = "Depth (m)")
  # map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName)
  
  rm(DataSpecies, lobster); myBiomodData; plot(myBiomodData)
  
}

myBiomodOption <- BIOMOD_ModelingOptions()

DataSplitTable <- BIOMOD_cv(myBiomodData, k = 5, rep = 2, do.full.models = F)
DataSplitTable.y <- BIOMOD_cv(myBiomodData, stratified.cv = F, stratify = "y", k = 2)
colnames(DataSplitTable.y)[1:2] <- c("RUN11", "RUN12")
DataSplitTable <- cbind(DataSplitTable, DataSplitTable.y)
head(DataSplitTable)

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
load("final_sdms.RData")
model = unique(mw$models)
model <- gsub(x = model, pattern = "MAXENT1", replacement = "MAXENT.Phillips")
model <- gsub(x = model, pattern = "MAXENT1", replacement = "MAXENT")

i = 2

for (i in 1:length(model)) {
  
  myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                      models = model[i],
                                      models.options = myBiomodOption,
                                      DataSplitTable = DataSplitTable,
                                      VarImport=0,
                                      models.eval.meth = c('ROC'),
                                      do.full.models = FALSE,
                                      modeling.id="test")
  
  eval <- get_evaluations(myBiomodModelOut,as.data.frame=T)
  
  eval$strat <- NA
  eval$strat[grepl("13",eval$Model.name)] <- "Full"
  eval$strat[!(grepl("11",eval$Model.name)|
                 grepl("12",eval$Model.name)|
                 grepl("13",eval$Model.name))] <- "Random"
  eval$strat[grepl("11",eval$Model.name)|grepl("12",eval$Model.name)] <- "Strat"
  
  # png(paste0("/Users/Kisei/Desktop/", model[i], ".png"), height = 5, width = 5, units = "in", res = 500)
  pdf(paste0("/Users/Kisei/Desktop/", sp, "_", model[i], ".pdf"), height = 5, width = 5)
  boxplot(eval$Testing.data ~ eval$strat, ylab="AUC", main = paste0(sp, "_", model[i]), pch = 20)
  dev.off()
  
}
