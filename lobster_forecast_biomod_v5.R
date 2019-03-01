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
library(Rmisc)
library(svMisc)
library(pals)
library(doParallel)
library(data.table)

dir = "/Users/Kisei/"
# dir = "/Users/Chenlab/"

# Prepare Biomod input data -----------------------------------------------
load(paste0(dir, "Google Drive/R/Biomod/lobster/lobster_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data

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
load(paste0(dir, "Google Drive/R/Biomod/lobster/myExpl_Kriged_maxDist0.08_res0.05.RData"))

par(mfrow = c(1,3))
plot(myExpl$var1.pred.1, main = "Bottom Temp (deg C)")
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
plot(myExpl$var1.pred.2, main = "Bottom Salt (ppt)")
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()
plot(myExpl$var1.pred.3, main = "Depth (m)")
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE); box()

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

rm(DataSpecies, lobster); myBiomodData; plot(myBiomodData)

# Tune biomod modeling options --------------------------------------------
dir = "/Users/Kisei/"
dir = "/Users/Chenlab/"

cl = makeCluster(16); registerDoParallel(cl)
cl = makeCluster(12); registerDoParallel(cl)

myBiomodOption = BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar")))

source(paste0(dir, "Google Drive/Research/Charlie/scripts/BIOMOD_tuning_KRT.R"))

time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning_KRT(myBiomodData,
                                                         # models = "MAXENT.Phillips",
                                                         models.options = myBiomodOption,
                                                         env.ME = myExpl,
                                                         metric = "ROC",
                                                         metric.ME = "Mean.AUC",
                                                         n.bg.ME = ncell(myExpl)))
beep(2)

setwd(paste0(dir, "Google Drive/R/Biomod/lobster"))
save(Biomod.tuning, file = "Lobster_Biomod_Tuning_Results.RData")
stopCluster(cl)

load(paste0(dir, "Google Drive/R/Biomod/lobster/Lobster_Biomod_Tuning_Results.RData"))

#see difference
plot(Biomod.tuning$tune.CTA.rpart)
plot(Biomod.tuning$tune.CTA.rpart2)
plot(Biomod.tuning$tune.ANN)
plot(Biomod.tuning$tune.GAM)
plot(Biomod.tuning$tune.RF)
plot(Biomod.tuning$tune.MARS)
plot(Biomod.tuning$tune.GBM)
plot(Biomod.tuning$tune.FDA)

Biomod.tuning$models.options

# Run Biomod, make sure to save results ---------------------------

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
  
  MAXENT.Tsuruoka = list(l1_regularizer = 0,
                         l2_regularizer = 0,
                         use_sgd = FALSE,
                         set_heldout = 0,
                         verbose = FALSE),
  
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
#maxent.jar file cannot be in google drive, choose desktop or download folder

# set directory for saved data
setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))

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
             'MAXENT.Tsuruoka',
             'MAXENT.Phillips'),
  # models = "MAXENT.Phillips",
  models.options = myBiomodOption,
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

# myCurrentProj <- get_predictions(myBiomodProj)
# plot(myCurrentProj, zlim = c(0,1000))

myBiomodModelEval <- get_evaluations(myBiomodModelOut)

#save everything
save.image("Biomod_Results.RData")
# save.image("Biomod_Results_with_Default_Options.RData")

# Load the saved Biomod results, both lobster and scallop --------------------------------------------------------
sp = c("lobster", "scallop")[2]

setwd(paste0(dir, "/Google Drive/R/Biomod/", sp))
load("Biomod_Results.RData")

# Get TSS and ROC summaries -------------------------------------------------------
# summarize TSS and ROC scores of all models
tss = myBiomodModelEval["TSS","Testing.data",,,]; tss = as.data.frame(tss)
roc = myBiomodModelEval["ROC","Testing.data",,,]; roc = as.data.frame(roc)

# roc$mean = apply(roc, 1, mean)
# tss$mean = apply(tss, 1, mean)

get_variables_importance(myBiomodModelOut)

models = c('GLM',
           'GAM',
           'GBM',
           'CTA',
           'ANN',
           'SRE',
           'FDA',
           'MARS',
           'RF',
           'MAXENT.Tsuruoka',
           'MAXENT.Phillips')

tss = cbind(models, tss)
tss
# write_csv(tss, "TSS.csv")

# Compare model performance and select final models, both lobster and scallop -------------------------------
sp = c("lobster", "scallop")

for (i in 1:length(sp)) {
  
  setwd(paste0(dir, "/Google Drive/R/Biomod/", sp[[i]]))
  load("Biomod_Results.RData")
  
  #organize all SDMs by TSS and ROC
  model_eval = t(as.data.frame(myBiomodModelEval))
  seq = seq(1, nrow(model_eval), by = 4)
  model_eval = model_eval[seq,]
  model_eval = as.data.frame(model_eval)
  model_eval$combo = model_eval$TSS + model_eval$ROC
  model_eval=model_eval[order(-model_eval$combo),]
  row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Phillips", replacement = "MAXENT1")  
  row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Tsuruoka", replacement = "MAXENT2")  
  
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
  
  jpeg(paste0("/Users/Kisei/Desktop/model_evaluation_legend.jpg"), res = 500, height = 7, width = 2, units = "in")
  
  legend = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
    geom_point(size=4) +
    scale_shape_manual(values = rep(18, 11))+
    geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
    geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
    ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))
  
  legend = cowplot::get_legend(legend)
  grid.newpage()
  grid.draw(legend)
  dev.off()

  jpeg(paste0("/Users/Kisei/Desktop/model_evaluation_", sp[[i]], ".jpg"), res = 500, height = 7, width = 5, units = "in")
  
  if (sp[[i]] == "lobster") {
    p = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
      geom_point(size=4) +
      scale_shape_manual(values = rep(18, 11))+
      geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
      geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
      ggtitle("H. americanus")+
      ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+ 
      theme(legend.position="none")
  }
  
  if (sp[[i]] == "scallop") {
    p =  ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
      geom_point(size=4) +
      scale_shape_manual(values = rep(18, 11))+
      geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
      geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
      ggtitle("P. magellanicus")+
      ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+ 
      theme(legend.position="none")
    
  }
  
  print(p)
  
  # ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, label = Model)) +
  #   geom_point(size=2) + geom_text(aes(label = Model), hjust = -0.1, vjust = -0.8, size = 5) +
  #   scale_shape_manual(values = rep(18, 11))+
  #   geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
  #   geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) +
  #   ggtitle("") +
  #   xlim(0.5, 1) +
  #   ylim(0.3, 0.8) +
  #   theme_classic(base_size = I(20)) +
  #   theme(legend.position="none")
  
  dev.off()
  
  #remove SDMs with TSS < 0.5
  evals = subset(evals, TSS > 0.5)
  
  # weighting all models by TSS, ROC, and combined weight
  weight = evals[c(1:nrow(evals)),5] #4 is ROC, 5 is combo
  scale <- function(x){(x-min(x))/(max(x)-min(x))} #scale weight to 0-1
  weight = scale(weight)
  models = evals$model
  run = evals$run
  mw = data.frame(models, run, weight)
  mw$SDM = paste0(mw$models, "_", mw$run)
  
  jpeg(paste0("/Users/Kisei/Desktop/model_weights_legend.jpg"), res = 500, height = 7, width = 2, units = "in")
  
  legend = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
    geom_point(aes(size = weight)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  legend = cowplot::get_legend(legend)
  grid.newpage()
  grid.draw(legend)
  dev.off()
  
  jpeg(paste0("/Users/Kisei/Desktop/model_weights_", sp[[i]], ".jpg"), res = 500, height = 7, width = 5, units = "in")
  
  if (sp[[i]] == "lobster") {
    p = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
      # geom_point(size = 5) + 
      geom_point(aes(size = weight)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle("H. americanus")+ 
      theme(legend.position="none")
  }
  
  if (sp[[i]] == "scallop") {
    p = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
      # geom_point(size = 5) + 
      geom_point(aes(size = weight)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle("P. magellanicus")+ 
      theme(legend.position="none")
  }
  
  print(p)
  dev.off()
  
  # good_sdms = models with TSS > 0.5
  a = paste0(sp[[i]], "_AllData_RUN")
  good_sdms = paste0(a, mw$run, "_", mw$models)
  good_sdms = gsub("MAXENT1", "MAXENT.Phillips", good_sdms)
  good_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", good_sdms)
  
  d = myBiomodProj@proj@val@layers
  dd = stack(d)
  ddd = subset(dd,good_sdms)
  
  jpeg(paste0("/Users/Kisei/Desktop/good_sdms_", sp[[i]], ".jpg"), res = 500, height = 10, width = 15, units = "in")
  plot(ddd)
  dev.off()
  
  # create final_sdms = visually check each model output and remove any unrealistic models
  mw = mw[mw$models!="CTA",]
  mw = mw[mw$models!="MAXENT1",]
  mw = mw[mw$models!="GBM",]
  
  final_sdms = paste0(a, mw$run, "_", mw$models)
  final_sdms = gsub("MAXENT1", "MAXENT.Phillips", final_sdms)
  final_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", final_sdms)
  d = myBiomodProj@proj@val@layers
  dd = stack(d)
  ddd = subset(dd,final_sdms)
  
  #do weighted mean
  avg = weighted.mean(ddd, mw$weight); dev.off(); 
  plot(avg, zlim = c(0,1000), col = parula(100), xaxt='n', yaxt='n')
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1)
  
  avg@data@values = (avg@data@values - avg@data@min)/(avg@data@max - avg@data@min)
  
  jpeg(paste0("/Users/Kisei/Desktop/Biomod_Parsimonious_Model_Prediction", sp[[i]], ".jpg"), res = 500, height = 10, width = 10, units = "in")
  plot(avg, zlim = c(0,1), col = parula(100), xaxt='n', yaxt='n')
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1, las = 1)
  legend("topleft", ifelse(sp[[i]] == "lobster", "H. americanus", "P. magellanicus"), bty = "n", cex = 2)
  dev.off()
  
  #save list of final SDMs and weights. 
  save(final_sdms, mw, file = paste0(sp[[i]], "_final_sdms.RData"))
  
}

# save response surves ----------------------------------------------------
dir = "/Users/Kisei/"
dir = "/Users/Chenlab/"

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
load("Biomod_Results.RData")

myModel <- BIOMOD_LoadModels(myBiomodModelOut)

for( i in 1:length(myBiomodModelOut@models.computed)){
  
  interval = length(myBiomodModelOut@models.computed)
  
  d <- response.plot2(models  = myModel[seq(i,interval,11)], #model, total # of runs, total # of models
                      Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                      show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                      do.bivariate = FALSE,
                      fixed.var.metric = 'median',
                      col = pals::parula(interval/6),
                      legend = TRUE,
                      data_species = get_formal_data(myBiomodModelOut,'resp.var'))
  
  if(i == 1) model = "GLM"
  if(i == 2) model = "GAM"
  if(i == 3) model = "GBM"
  if(i == 4) model = "CTA"
  if(i == 5) model = "ANN"
  if(i == 6) model = "SRE"
  if(i == 7) model = "FDA"
  if(i == 8) model = "MARS"
  if(i == 9) model = "RF"
  if(i == 10) model = "MAXENT.Tsuruoka"
  if(i == 11) model = "MAXENT.Phillips"
  
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


# load cm2.6 data and project habitat change year 1-80 -------------------------------------
dir = "/Users/Kisei/"
dir = "/Users/Chenlab/"

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
load("Biomod_Results.RData")
# load(paste0(dir, "/Google Drive/R/Biomod/lobster/Biomod_Results_with_Default_Options.RData"))
load(paste0(dir, "/Google Drive/R/Biomod/lobster/lobster_final_SDMs.RData"))

static = read_csv("CM_0.05.csv", col_names  = T)
load("distance_offshore_for_static_data_low_res.rda")

static$Distant_Offshore = df$distance

month = c("apr", "may", "jun", "sep", "oct", "nov")
season = c("spring", "fall")

period = c(9:18) #first 10 years
period = c(79:88) #last 10 years

for (i in 1:length(month)) {
  
  bt = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_", month[i], ".csv")), col_names  = T)
  bs = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_", month[i], ".csv")), col_names  = T)
  
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
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_apr.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_may.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_jun.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_apr.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_may.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_jun.csv")), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  if (season[i] == "fall") {
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_sep.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_oct.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_nov.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_sep.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_oct.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_nov.csv")), col_names  = T)
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
setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
load(paste0(dir, "/Google Drive/R/Biomod/lobster/Biomod_Results.RData"))
# load(paste0(dir, "/Google Drive/R/Biomod/lobster/Biomod_Results_with_Default_Options.RData"))
load(paste0(dir, "/Google Drive/R/Biomod/lobster/lobster_final_SDMs.RData"))

timestep = c("apr", "may", "jun", "sep", "oct", "nov")
season = c("spring", "fall")

for (k in 1:length(timestep)){
  
  month = timestep[k]
  
  # month = "apr"
  
  static = read_csv("CM_0.05.csv", col_names  = T)
  load("distance_offshore_for_static_data_low_res.rda")
  
  static$Distant_Offshore = df$distance
  
  bt = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_", month, ".csv")), col_names  = T)
  bs = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_", month, ".csv")), col_names  = T)
  
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
    
    myExplFuture = stack(btemp,bsalt,depth,lat,lon)
    
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
  
  depth <- rasterToPoints(depth, spatial = TRUE)
  depth <- as.data.frame(depth); colnames(depth)[1] = "depth"
  
  btemp <- rasterToPoints(btemp, spatial = TRUE)
  btemp <- as.data.frame(btemp);colnames(btemp)[1] = "btemp"
  
  bsalt <- rasterToPoints(bsalt, spatial = TRUE)
  bsalt <- as.data.frame(bsalt);colnames(bsalt)[1] = "bsalt"
  
  POS1 = merge(POS, depth)
  POS1 = merge(POS1, btemp)
  POS1 = merge(POS1, bsalt)
  
  POS = POS1[,c(1:2,83:85, 3:82)]
  write_csv(POS1, paste0(getwd(), "/", "Biomod_1_80_", month, ".csv"), col_names = T)
  
}#run by monthly timestep
for (k in 1:length(season)){
  
  season = season[k]
  
  if (season == "spring") {
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_apr.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_may.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_jun.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_apr.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_may.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_jun.csv")), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  if (season == "fall") {
    
    bt_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_sep.csv")), col_names  = T)
    bt_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_oct.csv")), col_names  = T)
    bt_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/temp/Delta_nov.csv")), col_names  = T)
    bt = rbindlist(list(bt_1,bt_2,bt_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
    bs_1 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_sep.csv")), col_names  = T)
    bs_2 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_oct.csv")), col_names  = T)
    bs_3 = read_csv(paste0(paste0(dir, "Google Drive/R/Biomod/lobster/salt/Delta_nov.csv")), col_names  = T)
    bs = rbindlist(list(bs_1,bs_2,bs_3))[,lapply(.SD,mean), list(Lon, Lat)]
    
  }
  
  static = read_csv("CM_0.05.csv", col_names  = T)
  load("distance_offshore_for_static_data_low_res.rda")
  
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
    
    myExplFuture = stack(btemp,bsalt,depth,lat,lon)
    
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
  
  depth <- rasterToPoints(depth, spatial = TRUE)
  depth <- as.data.frame(depth); colnames(depth)[1] = "depth"
  
  btemp <- rasterToPoints(btemp, spatial = TRUE)
  btemp <- as.data.frame(btemp);colnames(btemp)[1] = "btemp"
  
  bsalt <- rasterToPoints(bsalt, spatial = TRUE)
  bsalt <- as.data.frame(bsalt);colnames(bsalt)[1] = "bsalt"
  
  POS1 = merge(POS, depth)
  POS1 = merge(POS1, btemp)
  POS1 = merge(POS1, bsalt)
  
  POS = POS1[,c(1:2,83:85, 3:82)]
  write_csv(POS1, paste0(getwd(), "/", "Biomod_1_80_", season, ".csv"), col_names = T)
  
} #run by seasonal timesteps

# Loop 30 years, create dynamic environment dataset for migclim analysis with FVCOM-------------------------------------
setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
load("Biomod_Results.RData")
load("lobster_final_SDMs.RData")
load("fvcom_ts_1978-2013.RData")

fvcom_raster = function(df){
  
  coordinates(df) = ~Lon + Lat
  
  v = variogram(var1.pred ~ 1, df)
  auto = autofitVariogram(var1.pred ~ 1, df)
  g = gstat(formula = var1.pred ~ 1, model = auto$var_model, data = df, maxdist = 0.1)  
  
  xrange = range(df$Lon)
  yrange = range(df$Lat)
  grid= expand.grid(Lon = seq(from = xrange[1], to = xrange[2], by = 0.05), 
                    Lat = seq(from = yrange[1], to = yrange[2], by = 0.05))
  gridded(grid) = ~Lon + Lat
  
  p = predict(g, newdata = grid)
  p = raster(p)
  
  return(p)
}

timestep = c("apr", "may", "jun", "sep", "oct", "nov")

for (k in 1:length(timestep)){
  
  month = timestep[k]
  
  # month = "apr"
  # month = "may"
  # month = "jun"
  # month = "sep"
  # month = "oct"
  # month = "nov"
  
  if (month == "apr") {
    tt = t[, c(1:9, (rep(c(3), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(3), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  if (month == "may") {
    tt = t[, c(1:9, (rep(c(4), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(4), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  if (month == "jun") {
    tt = t[, c(1:9, (rep(c(5), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(5), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  if (month == "sep") {
    tt = t[, c(1:9, (rep(c(8), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(8), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  if (month == "oct") {
    tt = t[, c(1:9, (rep(c(9), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(9), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  if (month == "nov") {
    tt = t[, c(1:9, (rep(c(10), 36) + rep((1:36)*12 - 2, each = 1)))]
    ss = s[, c(1:9, (rep(c(10), 36) + rep((1:36)*12 - 2, each = 1)))]
    tt = t[,c(1:9, 16:45)]; ss = s[,c(1:9, 16:45)]
  }
  
  bt = tt; bs = ss; rm(tt,ss)
  
  num.cols <- c(6,10:39)
  bt[num.cols] <- sapply(bt[num.cols], as.numeric); bs[num.cols] <- sapply(bs[num.cols], as.numeric); rm(num.cols)
  
  bt$depth = bt$depth*-1
  bs$depth = bs$depth*-1
  bt = subset(bt, depth>0)
  bs = subset(bs, depth>0)
  bt = subset(bt, depth<400)
  bs = subset(bs, depth<400)
  
  for (i in 1:30){
    
    colnames(bt)[1:2] = c("Lat", "Lon"); colnames(bs)[1:2] = c("Lat", "Lon")
    
    bt = subset(bt, Lat <= 45.22589 & Lat >= 38.97589 & Lon <= -66.0226, Lon >= -75.6226)
    bs = subset(bs, Lat <= 45.22589 & Lat >= 38.97589 & Lon <= -66.0226, Lon >= -75.6226)
    
    depth = bt[,c(1,2,6)];colnames(depth)[3] = "var1.pred"
    btmed = bt[,c(1,2,(i+9))];colnames(btmed)[3] = "var1.pred"
    bsmed = bs[,c(1,2,(i+9))];colnames(bsmed)[3] = "var1.pred"
    lat = bt[,c(1,2,1)];colnames(lat)[3] = "var1.pred"
    lon = bt[,c(1,2,2)];colnames(lon)[3] = "var1.pred"
    
    btemp <- fvcom_raster(btmed)
    bsalt <- fvcom_raster(bsmed)
    
    load("lon-lat-depth_raster_layers_fvcom.RData")
    
    myExplFuture = stack(btemp,bsalt,depth,lat,lon)
    
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
    
    progress(i, progress.bar = T)
    Sys.sleep(0.01)
    if (i == 32) cat("Done!\n")
    
  }
  
  depth <- rasterToPoints(depth, spatial = TRUE)
  depth <- as.data.frame(depth); colnames(depth)[1] = "depth"
  
  btemp <- rasterToPoints(btemp, spatial = TRUE)
  btemp <- as.data.frame(btemp);colnames(btemp)[1] = "btemp"
  
  bsalt <- rasterToPoints(bsalt, spatial = TRUE)
  bsalt <- as.data.frame(bsalt);colnames(bsalt)[1] = "bsalt"
  
  POS1 = merge(POS, depth)
  POS1 = merge(POS1, btemp)
  POS1 = merge(POS1, bsalt)
  
  POS = POS1[,c(1:2, 33:35, 3:32)]
  
  write_csv(POS1, paste0("Biomod_1_30_", month, ".csv"), col_names = T)
  
}