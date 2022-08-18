rm(list = ls())

# load the library
library(biomod2)
library(raster)
library(readr)
library(colorRamps)
library(automap)
library(beepr)
library(gstat)
library(cowplot)
library(maps)
# library(Rmisc)
library(svMisc)
library(pals)
library(doParallel)
library(data.table)
library(grid)
library(ENMeval)
library(ggpubr)
library(devtools)

#install.packages("devtools")

#devtools::install_github('biomodhub/biomod2')
##

op = c("tuned", "default")[2]

select = dplyr::select

sp = read_csv("Loligo_1992_present_NEFSC_SPRING.csv") #load survey data
fl = read_csv("Loligo_1992_present_NEFSC_FALL.csv") #load survey data

sp = sp %>% select(lon, lat, depth, salinity, btemp, num) %>% as.data.frame()
fl = fl %>% select(lon, lat, depth, salinity, btemp, num) %>% as.data.frame()

# df = rbind(sp, fl); rm(sp, fl)
df = sp; rm(sp, fl)


btemp = df %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>% 
  group_by(lon, lat) %>% 
  summarise(btemp = median(btemp, na.rm = T)) %>% 
  rasterFromXYZ()

salinity = df %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>% 
  group_by(lon, lat) %>% 
  summarise(salinity = median(salinity, na.rm = T)) %>% 
  rasterFromXYZ()

depth = df %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>% 
  group_by(lon, lat) %>% 
  summarise(depth = median(depth, na.rm = T)) %>% 
  rasterFromXYZ()

DataSpecies = df[,c("lat", "lon", "num")] #Lat, Lon, Catch
DataSpecies$num = ifelse(DataSpecies$num > 1, 1, 0) 
names(DataSpecies)[1] = 'Y_WGS84'
names(DataSpecies)[2] = 'X_WGS84'
names(DataSpecies)[3] = 'squid'

myRespName = 'squid' # the name of studied species
myResp = DataSpecies[,myRespName] # the presence/absences vector
myRespXY = DataSpecies[,c("X_WGS84","Y_WGS84")] # the XY coordinates of species data

load("myExpl_Kriged_maxDist0.08_res0.05.RData")
load("slope.RData")
myExpl = stack(myExpl, slope)
names(myExpl)[6] = "var1.pred.6"
myExpl = stack(btemp, salinity, depth)
names(myExpl) = c("var1.pred.1", "var1.pred.2", "var1.pred.3")

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,##error: no psuedo absences selection! ! no data has been set aside for modeling evaluation! Some NAs have been automatically removed from your data###
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

rm(DataSpecies, squid); myBiomodData; plot(myBiomodData)




# Tune and save biomod modeling options --------------------------------------------

file.copy(from = "maxent.jar",
          to = paste0("/Users/", Sys.info()[7], "/Desktop/maxent.jar"))

myBiomodOption = BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = paste0("/Users/", Sys.info()[7], "/Desktop/maxent.jar")))

source("Biomod_Tuning_KRT.R")

time.seq <- system.time(
  Biomod.tuning <- BIOMOD_tuning_KRT(myBiomodData,
                                     models = c('GLM', 'GAM', 'GBM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips'), ###kisei had this hashtagged out
                                     #models = "CTA",
                                     models.options = myBiomodOption,
                                     env.ME = myExpl,
                                     # metric = c("ROC", "TSS"),
                                     metric = c("ROC"),
                                     # metric = c("TSS"),
                                     metric.ME = "Mean.AUC",
                                     n.bg.ME = ncell(myExpl))
)


beep(2)

save(Biomod.tuning, file = paste0("/Users/", Sys.info()[7], "/Desktop/Squid_Biomod_Tuning_Results.RData"))

# load tuned results ------------------------------------------------------
load(paste0("/Users/", Sys.info()[7], "/Desktop/Squid_Biomod_Tuning_Results.RData"))

#see difference
plot(Biomod.tuning$tune.CTA.rpart)
plot(Biomod.tuning$tune.CTA.rpart2)
plot(Biomod.tuning$tune.ANN)
plot(Biomod.tuning$tune.GAM)
plot(Biomod.tuning$tune.RF)
plot(Biomod.tuning$tune.MARS)
plot(Biomod.tuning$tune.GBM)
plot(Biomod.tuning$tune.GLM)
plot(Biomod.tuning$tune.FDA)
#eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results)

Biomod.tuning$tune.GLM
Biomod.tuning$tune.GAM
Biomod.tuning$tune.GBM
Biomod.tuning$tune.ANN
Biomod.tuning$tune.MARS
Biomod.tuning$tune.FDA
Biomod.tuning$tune.RF



Biomod.tuning$models.options

#####set biomod_modeling options for squid ---------------------
if (op == "default") {
  
  myBiomodOption <- BIOMOD_ModelingOptions(
    MAXENT.Phillips = list(path_to_maxent.jar = paste0("/Users/", Sys.info()[7], "/Desktop/maxent.jar"))
  )
}

if (op == "tuned") {
  
  myBiomodOption <- BIOMOD_ModelingOptions(
    
    GLM = list( type = 'simple',
                interaction.level = 2,#auto says 0
                myFormula = lobster ~ var1.pred.1 + var1.pred.2 + var1.pred.3 + var1.pred.4 + var1.pred.5,
                test = 'none',
                family = binomial(link = 'logit'),
                mustart = 0.5,
                control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE) ),
    
    GBM = list( distribution = 'bernoulli',
                n.trees = 5000,
                interaction.depth = 9,
                n.minobsinnode = 5,
                shrinkage = 0.005,
                bag.fraction = 0.5,
                train.fraction = 1,
                cv.folds = 3,
                keep.data = FALSE,
                verbose = FALSE,
                perf.method = 'cv'),
    
    GAM = list( algo = 'GAM_mgcv',
                type = 's_smoother',
                k = -1,
                interaction.level = 2, #auto says 0
                myFormula = NULL,
                family = binomial(link = 'logit'),
                method = 'GCV.Cp',
                optimizer = c('outer','newton'),
                select = FALSE,
                knots = NULL,
                paraPen = NULL,
                control = list(nthreads = 1, 
                               irls.reg = 0, 
                               epsilon = 1e-07, 
                               maxit = 200, 
                               trace = FALSE,
                               mgcv.tol = 1e-07,
                               mgcv.half = 15, 
                               rank.tol = 1.49011611938477e-08, 
                               nlm = list(ndigit=7, 
                                          gradtol=1e-06, 
                                          stepmax=2, 
                                          steptol=1e-04, 
                                          iterlim=200, 
                                          check.analyticals=0),
                               optim = list(factr=1e+07), 
                               newton = list(conv.tol=1e-06, 
                                             maxNstep=5, 
                                             maxSstep=2, 
                                             maxHalf=30, 
                                             use.svd=0), 
                               outerPIsteps = 0, 
                               idLinksBases = TRUE, 
                               scalePenalty = TRUE, 
                               efs.lspmax = 15, 
                               efs.tol = 0.1, 
                               keepData = FALSE, 
                               # scale.est = fletcher, 
                               edge.correct = FALSE) ),
    
    CTA = list( method = 'class',
                parms = 'default',
                cost = NULL,
                control = list(xval = 5, 
                               minbucket = 5, 
                               minsplit = 5, 
                               cp = list(cp=0.000211202162710146), 
                               maxdepth = list(maxdepth=7)) ),
    
    ANN = list( NbCV = 5,
                size = 8,
                decay = 0.001,
                rang = 0.1,
                maxit = 500),
    
    SRE = list( quant = 0.025),
    
    FDA = list( method = 'mars',
                add_args = list(degree = 2, nk = 20)),
    
    MARS = list( type = 'simple',
                 interaction.level = 0,
                 myFormula = NULL,
                 nk = 20,
                 penalty = 2,
                 thresh = 0.001,
                 nprune = NULL,
                 pmethod = 'backward'),
    
    RF = list( do.classif = TRUE,
               ntree = 500,
               mtry = 4,
               nodesize = 5,
               maxnodes = NULL),
    
    #MAXENT.Tsuruoka = list(l1_regularizer = 0,
    # l2_regularizer = 0,
    #use_sgd = FALSE,
    #set_heldout = 0,
    #verbose = FALSE),
    
    MAXENT.Phillips = list( path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar"),
                            memory_allocated = 2048, #auto says 512
                            background_data_dir = 'default',
                            maximumbackground = 'default',
                            maximumiterations = 200, #originally 10000
                            visible = FALSE,
                            linear = TRUE,
                            quadratic = TRUE,
                            product = TRUE,
                            threshold = TRUE,
                            hinge = TRUE,
                            lq2lqptthreshold = 80,
                            l2lqthreshold = 10,
                            hingethreshold = 15,
                            beta_threshold = -1,
                            beta_categorical = -1,
                            beta_lqp = -1,
                            beta_hinge = -1,
                            betamultiplier = 1,
                            defaultprevalence = 0.5)
  )
  
} 

# Run Biomod, make sure to save results ---------------------------

#save on your desktop first
setwd(paste0("/Users/", Sys.info()[7], "/Desktop"))

# Computing the models
myBiomodModelOut = BIOMOD_Modeling(
  myBiomodData,
  models = c(
    'GLM',
    'GAM',
    'GBM',
    'CTA',
    'ANN',
    'SRE',
    'FDA',
    'MARS',
    'RF'),
  # 'MAXENT.Tsuruoka',
  #'MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval = 3,
  DataSplit = 80,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling", sep = ""))

# projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myCurrentProj <- get_predictions(myBiomodProj)
plot(myCurrentProj, zlim = c(0,1000))

myBiomodModelEval <- get_evaluations(myBiomodModelOut)

# save everything
save.image("Biomod_Results.RData")


# Load the saved Biomod results, compare results --------------------------------------------------------
load("Biomod_Results.RData")

# Get TSS and ROC summaries
# summarize TSS and ROC scores of all models
tss = myBiomodModelEval["TSS","Testing.data",,,]; tss = as.data.frame(tss)
roc = myBiomodModelEval["ROC","Testing.data",,,]; roc = as.data.frame(roc)

# roc$mean = apply(roc, 1, mean)
# tss$mean = apply(tss, 1, mean)

get_variables_importance(myBiomodModelOut)

models = c(
  'GLM',
  'GAM',
  'GBM',
  'CTA',
  'ANN',
  'SRE',
  'FDA',
  'MARS',
  'RF')
# 'MAXENT.Tsuruoka',
# 'MAXENT.Phillips')

tss = cbind(models, tss)
colnames(tss) = c("SDM", "TSS_RUN1","TSS_RUN2","TSS_RUN3")
roc = cbind(models, roc)
colnames(roc) = c("SDM", "ROC_RUN1","ROC_RUN2","ROC_RUN3")
tss_roc = merge(tss, roc)
#tss_roc = subset(tss_roc, SDM!="MAXENT.Tsuruoka")
#tss_roc$SDM = gsub("MAXENT.Phillips", "MAXENT", tss_roc$SDM)
tss_roc






####
##
# ######Compare model performance and select final models #######
#organize all SDMs by TSS and ROC
model_eval = t(as.data.frame(myBiomodModelEval))
seq = seq(1, nrow(model_eval), by = 4)
model_eval = model_eval[seq,]
model_eval = as.data.frame(model_eval)
model_eval$combo = model_eval$TSS + model_eval$ROC
model_eval=model_eval[order(-model_eval$combo),]
#row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Phillips", replacement = "MAXENT1")  
#row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Tsuruoka", replacement = "MAXENT2")  

sdm = data.frame(model = rep(NA,nrow(model_eval)),run = rep(NA,nrow(model_eval)))

for (j in 1:nrow(model_eval)){
  ff = strsplit(row.names(model_eval[j,]), '[._]')[[1]]
  sdm[j,1] = ff[3]
  sdm[j,2] = gsub("RUN", "", ff[length(ff)-1]) 
} 

evals = cbind(sdm,model_eval);rownames(evals) <- c()

d1 = summarySE(evals, measurevar = "TSS", groupvars = c("model"))
d2 = summarySE(evals, measurevar = "ROC", groupvars = c("model"))

eval = merge(d1[c("model","TSS","se")],d2[c("model","ROC","se")], by = "model")
names(eval) = c("Model", "TSS", "SE2", "ROC", "SE1")





png(filename = "ensemble_model_evaluation.png", res = 600, height = 6, width = 6, units = "in")
ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
  geom_point(size=4) +
  scale_shape_manual(values = rep(18, 11))+
  geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
  geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
  ggtitle("Ensemble Model Evaluation")
#dev.off()

#remove SDMs with TSS < 0.5
evals = subset(evals, TSS > 0.5)

# weighting all models by TSS, ROC, and combined weight
weight = evals[c(1:nrow(evals)),5] #4 is ROC, 5 is combo
scale <- function(x){(x-min(x))/(max(x)-min(x))} #scale weight to 0-1
weight = scale(weight)
models = evals$model
run = evals$run
mw = data.frame(models, run, weight)####error msg. arguments imply differing number of rows: 0,1
mw$SDM = paste0(mw$models, "_", mw$run)


#jpeg("model_weights_tunedS3.jpg", res = 500, height = 10, width = 25, units = "in")
ggplot(mw, aes(x=SDM, y=weight, color=models)) + geom_point(size = 5) 
#dev.off()

#mww = mw[-c(4, 5, 6, 8, 9, 13), ]       #for plotting purposes only 
#png("model_weightsSCAL.png", res = 500, height = 5, width = 7, units = "in")


#dev.off()



# good_sdms = models with TSS > 0.5
a = "squid_AllData_RUN"
good_sdms = paste0(a, mw$run, "_", mw$models)
good_sdms = gsub("MAXENT1", "MAXENT.Phillips", good_sdms)
good_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", good_sdms)

d = myBiomodProj@proj@val@layers
dd = stack(d)
ddd = subset(dd,good_sdms)
#plot(ddd[[21:24]])

#final_sdms = models with realistic outputs
#mw = mw[mw$models!="CTA" & mw$models!="GBM",] 
#mw = mw[mw$models!="CTA",] 

final_sdms = paste0(a, mw$run, "_", mw$models)
final_sdms = gsub("MAXENT1", "MAXENT.Phillips", final_sdms)
final_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", final_sdms)
d = myBiomodProj@proj@val@layers
dd = stack(d)
ddd = subset(dd,final_sdms)

#do weighted mean
avg = weighted.mean(ddd, mw$weight)

#jpeg("Biomod_Parsimonious_Model_Prediction_tunedS3.jpg", res = 500, height = 10, width = 10, units = "in")
plot(avg, zlim = c(0,1000), col = parula(100), xaxt='n', yaxt='n')
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
box()
degAxis(1, cex.axis = 1)
degAxis(2, cex.axis = 1)
#dev.off()

#save list of final SDMs and weights. 
#save(final_sdms, mw, file = "scallop_final_sdms_tuned.RData")

# save response curves ----------------------------------------------------

for( i in 1:11){
  
  if(i == 1) model = "GLM"
  if(i == 2) model = "GAM"
  if(i == 3) model = "GBM"
  if(i == 4) model = "CTA"
  if(i == 5) model = "ANN"
  if(i == 6) model = "SRE"
  if(i == 7) model = "FDA"
  if(i == 8) model = "MARS"
  if(i == 9) model = "RF"
  # if(i == 10) model = "MAXENT.Tsuruoka"
  #if(i == 11) model = "MAXENT.Phillips"
  
  myModel <- BIOMOD_LoadModels(myBiomodModelOut, models=model)
  
  myRespPlot2D <- response.plot2(models  = myModel,
                                 Data = get_formal_data(myBiomodModelOut,'expl.var'),
                                 show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                                 do.bivariate = FALSE,
                                 fixed.var.metric = 'mean',
                                 col = c("blue", "red"),
                                 legend = TRUE,
                                 data_species = get_formal_data(myBiomodModelOut,'resp.var'))
  
  
  
  write.csv(myRespPlot2D$var1.pred.1, paste0(model, "L1.csv"))
  write.csv(myRespPlot2D$var1.pred.2, paste0(model, "L2.csv"))
  write.csv(myRespPlot2D$var1.pred.3, paste0(model, "L3.csv"))
  write.csv(myRespPlot2D$var1.pred.4, paste0(model, "L4.csv"))
  write.csv(myRespPlot2D$var1.pred.5, paste0(model, "L5.csv"))
  
  
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


png("C:/Users/claire.ober/Desktop/SDM_Loligo_NY_BIGHT/shapes.png",
    width = 18, height = 12, units = "in", res = 300)

myRespPlot3D <- response.plot2(models  = myModel[1],
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = TRUE,
                               fixed.var.metric = 'median',
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                               display_title=FALSE)

dev.off()













#####haven't started this section yet####
###----environmental data----

######read your own in. dont need bottom stress or sediment. 
bStress = readOGR("C:/Users/claire.ober/Desktop/SDM_Loligo_NY_BIGHT/Environmental Data/Bottom Stress/MAB_median.shp")
sed = readOGR("C:/Users/claire.ober/Desktop/SDM_Loligo_NY_BIGHT/Environmental Data/Bottom Sediment/conmapsg.shp")




#####read in the environmental data files
ts1992 = read.csv("envirodata_csvfiles_SPRING/ts1992.csv")
ts1993 = read.csv("envirodata_csvfiles_SPRING/ts1993.csv")
ts1994 = read.csv("envirodata_csvfiles_SPRING/ts1994.csv")
ts1995 = read.csv("envirodata_csvfiles_SPRING/ts1995.csv")
ts1996 = read.csv("envirodata_csvfiles_SPRING/ts1996.csv")
ts1997 = read.csv("envirodata_csvfiles_SPRING/ts1997.csv")
ts1998 = read.csv("envirodata_csvfiles_SPRING/ts1998.csv")
ts1999 = read.csv("envirodata_csvfiles_SPRING/ts1999.csv")
ts2000 = read.csv("envirodata_csvfiles_SPRING/ts2000.csv")
ts2001 = read.csv("envirodata_csvfiles_SPRING/ts2001.csv")
ts2002 = read.csv("envirodata_csvfiles_SPRING/ts2002.csv")
ts2003 = read.csv("envirodata_csvfiles_SPRING/ts2003.csv")
ts2004 = read.csv("envirodata_csvfiles_SPRING/ts2004.csv")
ts2005 = read.csv("envirodata_csvfiles_SPRING/ts2005.csv")
ts2006 = read.csv("envirodata_csvfiles_SPRING/ts2006.csv")
ts2007 = read.csv("envirodata_csvfiles_SPRING/ts2007.csv")
ts2008 = read.csv("envirodata_csvfiles_SPRING/ts2008.csv")
ts2009 = read.csv("envirodata_csvfiles_SPRING/ts2009.csv")
ts2010 = read.csv("envirodata_csvfiles_SPRING/ts2010.csv")
ts2011 = read.csv("envirodata_csvfiles_SPRING/ts2011.csv")
ts2012 = read.csv("envirodata_csvfiles_SPRING/ts2012.csv")
ts2013 = read.csv("envirodata_csvfiles_SPRING/ts2013.csv")
ts2014 = read.csv("envirodata_csvfiles_SPRING/ts2014.csv")
ts2015 = read.csv("envirodata_csvfiles_SPRING/ts2015.csv")
ts2016 = read.csv("envirodata_csvfiles_SPRING/ts2016.csv")

##############################################################


myExpl = stack( system.file( "external/bioclim/current/bio3.grd",
                             
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             
                             package="biomod2"),
                
                system.file( "external/bioclim/current/bio7.grd",
                             
                             package="biomod2"),
                
                system.file( "external/bioclim/current/bio11.grd",
                             
                             package="biomod2"),
                
                system.file( "external/bioclim/current/bio12.grd",
                             
                             package="biomod2"))


#####################

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

rm(DataSpecies, squid); myBiomodData; plot(myBiomodData)

##############################


# Computing the models
myBiomodModelOut = BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM',
             'GAM',
             'GBM',
             'CTA',
             'ANN',
             'SRE',
             'FDA',
             'MARS',
             'RF',
             #'MAXENT.Tsuruoka',
             'MAXENT.Phillips'),
  # models = "CTA",
  #models.options = myBiomodOption,
  NbRunEval = 3,
  DataSplit = 80,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep = ""))





# projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myCurrentProj <- get_predictions(myBiomodProj)
plot(myCurrentProj, zlim = c(0,1000))

myBiomodModelEval <- get_evaluations(myBiomodModelOut)

#save everything
save.image("Biomod_Results_Scallop.RData")
# save.image("Biomod_Results_with_Default_Options.RData")



