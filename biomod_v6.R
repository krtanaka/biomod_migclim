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

# cl = makeCluster(detectCores()-1); registerDoParallel(cl)
op = c("tuned_lobster", "tuned_scallop", "default")[1]

# cl = makeCluster(64); registerDoParallel(cl)
# cl = makeCluster(8); registerDoParallel(cl)
# stopCluster(cl)

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
<<<<<<< HEAD
  load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
  load(paste0(dir, "biomod_migclim/data/slope.RData"))
  myExpl = stack(myExpl, slope)
  names(myExpl)[6] = "var1.pred.6"
=======
  load(paste0(dir, "Google Drive/R/Biomod/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
  
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
<<<<<<< HEAD
  load(paste0(dir, "biomod_migclim/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
  load(paste0(dir, "biomod_migclim/data/slope.RData"))
  myExpl = stack(myExpl, slope)
  names(myExpl)[6] = "var1.pred.6"
=======
  load(paste0(dir, "Google Drive/R/Biomod/data/myExpl_Kriged_maxDist0.08_res0.05.RData"))
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
  
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

# Tune and save biomod modeling options --------------------------------------------

myBiomodOption = BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar")))

source(paste0(dir, "biomod_migclim/BIOMOD_tuning_KRT.R"))

time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning_KRT(myBiomodData,
                                                         # models = "MAXENT.Phillips",
                                                         # models = "CTA",
                                                         models.options = myBiomodOption,
                                                         env.ME = myExpl,
                                                         # metric = c("ROC", "TSS"),
                                                         metric = c("ROC"),
                                                         # metric = c("TSS"),
                                                         metric.ME = "Mean.AUC",
                                                         n.bg.ME = ncell(myExpl)))
beep(2)

setwd(paste0(dir, "Desktop"))
save(Biomod.tuning, file = paste0(sp, "_Biomod_Tuning_Results.RData"))


# load tuned results ------------------------------------------------------
<<<<<<< HEAD
load("~/Google Drive/R/Biomod/lobster/lobster_Biomod_Tuning_Results_1.RData") #tuning without slope
load("~/Google Drive/R/Biomod/lobster/lobster_Biomod_Tuning_Results_2.RData") #tuning with slope
load("~/Google Drive/R/Biomod/scallop/scallop_Biomod_Tuning_Results_1.RData") #tuning without slope
load("~/Google Drive/R/Biomod/scallop/scallop_Biomod_Tuning_Results_2.RData") #tuning with slope
=======
setwd(paste0(dir, "Desktop"))
load(paste0(sp, "_Biomod_Tuning_Results.RData"))
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e

#see difference
plot(Biomod.tuning$tune.CTA.rpart)
plot(Biomod.tuning$tune.CTA.rpart2)
plot(Biomod.tuning$tune.ANN)
plot(Biomod.tuning$tune.GAM)
plot(Biomod.tuning$tune.RF)
plot(Biomod.tuning$tune.MARS)
plot(Biomod.tuning$tune.GBM)
<<<<<<< HEAD
plot(Biomod.tuning$tune.GLM)
plot(Biomod.tuning$tune.FDA)
eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results)

Biomod.tuning$tune.GLM
Biomod.tuning$tune.GAM
Biomod.tuning$tune.GBM
Biomod.tuning$tune.ANN
Biomod.tuning$tune.MARS
Biomod.tuning$tune.FDA
Biomod.tuning$tune.RF
=======
plot(Biomod.tuning$tune.FDA)
eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results)

>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
Biomod.tuning$models.options

# set biomod_modeling options for lobster and scallop ---------------------

<<<<<<< HEAD
if (op == "tuned") {
  
  if (so == "lobster") {
    
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
    
  } 
  
  if (sp == "scallop") {
    
    myBiomodOption <- BIOMOD_ModelingOptions(
      
      GLM = list( type = 'simple',
                  interaction.level = 0,#auto says 0
                  myFormula = scallop ~ var1.pred.1 + var1.pred.2 + var1.pred.3 + var1.pred.4 + var1.pred.5,
                  test = 'none',
                  family = binomial(link = 'logit'),
                  mustart = 0.5,
                  control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE) ),
      
      GBM = list( distribution = 'bernoulli',
                  n.trees = 2500,
                  interaction.depth = 9,
                  n.minobsinnode = 5,
                  shrinkage = 0.5,
                  bag.fraction = 0.5,
                  train.fraction = 1,
                  cv.folds = 3,
                  keep.data = FALSE,
                  verbose = FALSE,
                  perf.method = 'cv'),
      
      GAM = list( algo = 'GAM_mgcv',
                  type = 's_smoother',
                  k = -1,
                  interaction.level = 0, #auto says 0
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
                                 cp = list(cp=0.000525578135949544), 
                                 maxdepth = list(maxdepth=7)) ),
      
      ANN = list( NbCV = 5,
                  size = 8,
                  decay = 0.001,
                  rang = 0.1,
                  maxit = 500),
      
      SRE = list( quant = 0.025),
      
      FDA = list( method = 'mars',
                  add_args = list(degree = 2, nk = 18)),
      
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
                 mtry = 5,
                 nodesize = 5,
                 maxnodes = NULL),
      
      MAXENT.Tsuruoka = list(l1_regularizer = 0,
                             l2_regularizer = 0,
                             use_sgd = FALSE,
                             set_heldout = 0,
                             verbose = FALSE),
      
      MAXENT.Phillips = list( path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar"),
                              memory_allocated = 512, #auto says 512
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
  }
  
=======
if (op == "tuned_lobster") {
  
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
}

if (op == "tuned_scallop") {
  
  myBiomodOption <- BIOMOD_ModelingOptions(
    
    GLM = list( type = 'simple',
                interaction.level = 0,#auto says 0
                myFormula = scallop ~ var1.pred.1 + var1.pred.2 + var1.pred.3 + var1.pred.4 + var1.pred.5,
                test = 'none',
                family = binomial(link = 'logit'),
                mustart = 0.5,
                control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE) ),
    
    GBM = list( distribution = 'bernoulli',
                n.trees = 2500,
                interaction.depth = 9,
                n.minobsinnode = 5,
                shrinkage = 0.5,
                bag.fraction = 0.5,
                train.fraction = 1,
                cv.folds = 3,
                keep.data = FALSE,
                verbose = FALSE,
                perf.method = 'cv'),
    
    GAM = list( algo = 'GAM_mgcv',
                type = 's_smoother',
                k = -1,
                interaction.level = 0, #auto says 0
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
                               cp = list(cp=0.000525578135949544), 
                               maxdepth = list(maxdepth=7)) ),
    
    ANN = list( NbCV = 5,
                size = 8,
                decay = 0.001,
                rang = 0.1,
                maxit = 500),
    
    SRE = list( quant = 0.025),
    
    FDA = list( method = 'mars',
                add_args = list(degree = 2, nk = 18)),
    
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
               mtry = 5,
               nodesize = 5,
               maxnodes = NULL),
    
    MAXENT.Tsuruoka = list(l1_regularizer = 0,
                           l2_regularizer = 0,
                           use_sgd = FALSE,
                           set_heldout = 0,
                           verbose = FALSE),
    
    MAXENT.Phillips = list( path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar"),
                            memory_allocated = 512, #auto says 512
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
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
}

if (op == "default") {
  
<<<<<<< HEAD
  myBiomodOption <- BIOMOD_ModelingOptions(
=======
  myBiomodOption_default <- BIOMOD_ModelingOptions(
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
    
    MAXENT.Phillips = list( path_to_maxent.jar = paste0(dir, "Desktop/maxent.jar"))
    
  )
}

# Run Biomod, make sure to save results ---------------------------

#save on your desktop first
<<<<<<< HEAD
setwd(paste0(dir, "Desktop"))
=======
setwd("/Users/Kisei/Desktop")
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e

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
    'RF',
<<<<<<< HEAD
    # 'MAXENT.Tsuruoka',
    'MAXENT.Phillips'),
  # models = "GAM",
=======
    'MAXENT.Tsuruoka',
    'MAXENT.Phillips'),
  # models = "MAXENT.Tsuruoka",
  # models.options = myBiomodOption_default,
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
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
<<<<<<< HEAD
if(op == "default"){
  save.image("Biomod_Results_with_Default_Options.RData")
}else{
  save.image("Biomod_Results.RData")
}

=======
save.image("Biomod_Results.RData")
save.image("Biomod_Results_with_Default_Options.RData")

# Load the saved Biomod results, compare results --------------------------------------------------------

for (i in 1:length(sp)) {
  
  i = 1
  
  setwd(paste0(dir, "/biomod_migclim/", sp[[i]]))
  
  if (op == "default") {
    
    load("Biomod_Results_with_Default_Options.RData")
    
  }  else {
    
    load("Biomod_Results.RData")
    
  }
  
  # Get TSS and ROC summaries
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
  colnames(tss) = c("SDM", "TSS_RUN1","TSS_RUN2","TSS_RUN3")
  roc = cbind(models, roc)
  colnames(roc) = c("SDM", "ROC_RUN1","ROC_RUN2","ROC_RUN3")
  tss_roc = merge(tss, roc)
  tss_roc = subset(tss_roc, SDM!="MAXENT.Tsuruoka")
  tss_roc$SDM = gsub("MAXENT.Phillips", "MAXENT", tss_roc$SDM)
  tss_roc
  # write_csv(tss_roc, paste0("/Users/Kisei/Desktop/", sp[[i]], "_ROC_TSS.csv"))
  
  # Compare model performance and select final models, both lobster and scallop
  # organize all SDMs by TSS and ROC
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
  eval = subset(eval, Model!="MAXENT2")
  eval$Model = gsub("MAXENT1", "MAXENT", eval$Model)
  
  # jpeg(paste0("/Users/Kisei/Desktop/model_evaluation_legend.jpg"), res = 500, height = 2.5, width = 1, units = "in")
  # legend = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
  #   geom_point(size=4) +
  #   scale_shape_manual(values = rep(18, 11))+
  #   geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
  #   geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
  #   ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95)) + 
  #   theme(legend.position = "right", legend.title = element_blank())
  # 
  # legend = cowplot::get_legend(legend)
  # grid.newpage()
  # grid.draw(legend)
  # dev.off()
  
  jpeg(paste0("/Users/Kisei/Desktop/model_evaluation_", sp[[i]], ".jpg"), res = 500, height = 5, width = 3, units = "in")
  
  if (sp[[i]] == "lobster") {
    
    p = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
      geom_point(size=4) +
      scale_shape_manual(values = rep(18, 11))+
      geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
      geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
      ggtitle("American lobster")+
      ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+ 
      theme_classic() + 
      theme(legend.position="none")    
  }
  
  if (sp[[i]] == "scallop") {
    
    p =  ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
      geom_point(size=4) +
      scale_shape_manual(values = rep(18, 11))+
      geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
      geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
      ggtitle("Sea scallop")+
      ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+
      theme_classic() + 
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
  evals = subset(evals, ROC > 0.8)
  
  # weighting all models by TSS, ROC, and combined weight
  weight = evals[c(1:nrow(evals)),5] #4 is ROC, 5 is combo
  scale <- function(x){(x-min(x))/(max(x)-min(x))} #scale weight to 0-1
  weight = scale(weight)
  models = evals$model
  run = evals$run
  mw = data.frame(models, run, weight)
  mw$SDM = paste0(mw$models, "_", mw$run)
  
  # jpeg(paste0("/Users/Kisei/Desktop/model_weights_legend.jpg"), res = 500, height = 5, width = 1, units = "in")
  # legend = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
  #   geom_point(aes(size = weight)) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # legend = cowplot::get_legend(legend)
  # grid::grid.newpage()
  # grid::grid.draw(legend)
  # dev.off()
  
  jpeg(paste0("/Users/Kisei/Desktop/model_weights_", sp[[i]], ".jpg"), res = 500, height = 7, width = 5, units = "in")
  
  if (sp[[i]] == "lobster") {
    p = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
      # geom_point(size = 5) + 
      geom_point(aes(size = weight)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle("American lobster")+ 
      theme_classic() + 
      theme(legend.position="none")  
  }
  
  if (sp[[i]] == "scallop") {
    p = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
      # geom_point(size = 5) + 
      geom_point(aes(size = weight)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ggtitle("Sea scallop")+
      theme_classic() + 
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
  
  #all sdm > 0.5 TSS
  # jpeg(paste0("/Users/Kisei/Desktop/good_sdms_", sp[[i]], ".jpg"), res = 500, height = 10, width = 20, units = "in")
  plot(ddd, col = matlab.like(100), nc = 6)
  # dev.off()
  
  # create final_sdms = visually check each model output and remove any unrealistic models
  
  if (sp == "lobster") {
    mw = mw[mw$models!="CTA",]
    mw = mw[mw$models!="MAXENT1",]
    # mw = mw[mw$models!="GBM",]
  }
  
  if (sp == "scallop") {
    mw = mw[mw$models!="CTA",]
    # mw = mw[mw$models!="MAXENT1",]
    # mw = mw[mw$models!="GBM",]
  }
  
  final_sdms = paste0(a, mw$run, "_", mw$models)
  final_sdms = gsub("MAXENT1", "MAXENT.Phillips", final_sdms)
  final_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", final_sdms)
  d = myBiomodProj@proj@val@layers
  dd = stack(d)
  ddd = subset(dd,final_sdms)
  
  #do weighted mean
  avg = weighted.mean(ddd, mw$weight); dev.off(); 
  avg@data@values = (avg@data@values - avg@data@min)/(avg@data@max - avg@data@min)
  
  jpeg(paste0("/Users/Kisei/Desktop/Biomod_Parsimonious_Model_Prediction", sp[[i]], ".jpg"), res = 500, height = 5, width = 10, units = "in")
  
  par(mfrow = c(1,2), mar = c(3,3,3,5))
  plot(avg, zlim = c(0,1), col = matlab.like(100), 
       main = ifelse(sp[[i]] == "lobster", "American lobster", "Sea scallop"),
       xaxt='n', yaxt='n', 
       legend = F)
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1, las = 1)
  legend("bottomright", "Mean", bty = "n", cex = 1.4)
  
  #plot standard error from final SDMs
  source("/Users/Kisei/Google Drive/R/misc/color palette function.R")
  steps = c("blue", "red")
  pal = color.palette(steps, space="rgb")
  error_col = pal
  
  plot(calc(ddd, fun=sd)*0.001, 
       col = matlab.like(100), 
       zlim = c(0,1),
       axes = F)
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1, las = 2)
  legend("bottomright", "Standard deviation", bty = "n", cex = 1.4)
  
  dev.off()
  
  # #save list of final SDMs and weights. 
  # if (op == "default") {
  #   
  #   save(final_sdms, mw, file = "final_sdms_default.RData")
  #   
  # }else{
  #   
  #   save(final_sdms, mw, file = "final_sdms.RData")
  #   
  # }
  
}

# save response surves ----------------------------------------------------
sp = c("lobster", "scallop")[2]

load(paste0(dir, "biomod_migclim/", sp[[i]], "/data_and_variables_for_response_curves.RData"))
load(paste0(dir, "biomod_migclim/", sp[[i]], "/Biomod_Results.RData"))

setwd(paste0(dir, "Google Drive/R/Biomod/", sp[[i]]))
myModel <- BIOMOD_LoadModels(myBiomodModelOut)

interval = length(myBiomodModelOut@models.computed)

for( i in 1:length(myBiomodModelOut@models.computed)){
  
  ## this was used with old biomod projections
  # d <- response.plot2(models  = myModel[seq(i,interval,11)], #model, total # of runs, total # of models
  #                     Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  #                     show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  #                     do.bivariate = FALSE,
  #                     fixed.var.metric = 'median',
  #                     col = pals::parula(interval/6),
  #                     legend = TRUE,
  #                     data_species = get_formal_data(myBiomodModelOut,'resp.var'))
  
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


# load cm2.6 data and project habitat change year 1-80 -------------------------------------
setwd(paste0(dir, "/Google Drive/R/biomod/", sp))
load("Biomod_Results.RData")
load("final_sdms.RData")

static = read_csv("/Users//kisei/Google Drive/R/Biomod/data/CM_0.05.csv", col_names  = T)
load("/Users//kisei/Google Drive/R/Biomod/data/distance_offshore_for_static_data_low_res.rda")

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

setwd(paste0(dir, "/Google Drive/R/Biomod/", sp)); load("Biomod_Results.RData"); load("final_sdms.RData")
setwd(paste0(dir, "/Desktop/")); load(paste0(dir, "/Desktop/", sp,"/Biomod_Results.RData")); load(paste0(dir, "/Desktop/", sp,"/final_sdms.RData"))

season = c("spring", "fall", "annual")[3]
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
  
  static = read_csv("/Users/kisei/Google Drive/R/Biomod/data/CM_0.05.csv", col_names  = T)
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

>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
