# save response surves

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

sp = c("lobster", "scallop")[2]

#if you need to plot results with slope
load(paste0(dir, "Desktop/", sp, "/data_and_variables_for_response_curves.RData"))
load(paste0(dir, "Desktop/", sp, "/Biomod_Results.RData"))
setwd(paste0(dir, "Desktop"))#you must point to the parent directory, e.g. Desktop, not Desktop/lobster

#if you need to plot results without slope
load(paste0(dir, "Google Drive/R/Biomod/", sp, "/data_and_variables_for_response_curves.RData"))
load(paste0(dir, "Google Drive/R/Biomod/", sp, "/Biomod_Results.RData"))
setwd("~/Google Drive/R/Biomod") #you must point to the parent directory, e.g. Desktop, not Desktop/lobster

myModel <- BIOMOD_LoadModels(myBiomodModelOut)

interval = length(myBiomodModelOut@models.computed)

load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
load(paste0(dir, "biomod_migclim/data/slope.RData"))
myExpl = stack(myExpl, slope)
names(myExpl)[6] = "var1.pred.6"
Data = rasterToPoints(myExpl)
Data = Data[,c("var1.pred.1", "var1.pred.2", "var1.pred.3", "var1.pred.4", "var1.pred.5", "var1.pred.6")]
Data = as.data.frame(Data)
Data = Data[complete.cases(Data), ]

show.variables = c("var1.pred.1", "var1.pred.2", "var1.pred.3", "var1.pred.4", "var1.pred.5", "var1.pred.6")

for( i in 1:length(myBiomodModelOut@models.computed)){
  
  i = 10
  
  d = response.plot2(models  = myModel[seq(i,interval,10)], #model, total # of runs, total # of models
                     Data = Data, 
                     show.variables= show.variables,
                     do.bivariate = FALSE,
                     fixed.var.metric = 'median',
                     col = pals::parula(interval/6),
                     legend = TRUE)
  
  if(i == 1) model = "GLM"
  if(i == 2) model = "GAM"
  if(i == 3) model = "GBM"
  if(i == 4) model = "CTA"
  if(i == 5) model = "ANN"
  if(i == 6) model = "SRE"
  if(i == 7) model = "FDA"
  if(i == 8) model = "MARS"
  if(i == 9) model = "RF"
  if(i == 10) model = "MAXENT"
  # if(i == 10) model = "MAXENT.Tsuruoka"
  # if(i == 11) model = "MAXENT.Phillips"
  
  write.csv(d, paste0("Biomod_", model, ".csv"))
  
}

myRespPlot2D <- response.plot2(models = myModel,
                               # models  = myModel[seq(2,33,11)], #model, total # of runs, interval
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'mean',
                               col = parula(length(myModel)),
                               legend = T,
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'))

myRespPlot3D <- response.plot2(models  = myModel[1],
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = TRUE,
                               fixed.var.metric = 'median',
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                               display_title=FALSE)


