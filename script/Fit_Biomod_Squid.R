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

op = c("tuned", "default")[2]

select = dplyr::select

sp = read_csv("data/Loligo_1992_present_NEFSC_SPRING.csv") #load survey data
fl = read_csv("data/Loligo_1992_present_NEFSC_FALL.csv") #load survey data

sp = sp %>% select(lon, lat, depth, salinity, btemp, num) %>% as.data.frame()
fl = fl %>% select(lon, lat, depth, salinity, btemp, num) %>% as.data.frame()

df = rbind(sp, fl); rm(sp, fl)

DataSpecies = df[,c("lat", "lon", "num")] #Lat, Lon, Catch
DataSpecies$num = ifelse(DataSpecies$num > 1, 1, 0) #try chaniging presence/absence threashold to make it less sensitive
names(DataSpecies)[1] = 'Y_WGS84'
names(DataSpecies)[2] = 'X_WGS84'
names(DataSpecies)[3] = 'squid'

myRespName = 'squid' # the name of studied species
myResp = DataSpecies[,myRespName] # the presence/absences vector
myRespXY = DataSpecies[,c("X_WGS84","Y_WGS84")] # the XY coordinates of species data

load("data/myExpl_Kriged_maxDist0.08_res0.05.RData")
load("data/slope.RData")
myExpl = stack(myExpl, slope)
names(myExpl)[6] = "var1.pred.6"

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

rm(DataSpecies, lobster); myBiomodData; plot(myBiomodData)

# Tune and save biomod modeling options --------------------------------------------

file.copy(from = "maxent.jar",
          to = paste0("/Users/", Sys.info()[7], "/Desktop/maxent.jar"))

myBiomodOption = BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = paste0("/Users/", Sys.info()[7], "/Desktop/maxent.jar")))

source("script/Biomod_Tuning_KRT.R")

time.seq <- system.time(
  Biomod.tuning <- BIOMOD_tuning_KRT(myBiomodData,
                                     # models = "MAXENT.Phillips",
                                     # models = "CTA",
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
eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results)

Biomod.tuning$tune.GLM
Biomod.tuning$tune.GAM
Biomod.tuning$tune.GBM
Biomod.tuning$tune.ANN
Biomod.tuning$tune.MARS
Biomod.tuning$tune.FDA
Biomod.tuning$tune.RF

plot(Biomod.tuning$tune.FDA)
eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results)

Biomod.tuning$models.options

# set biomod_modeling options for lobster and scallop ---------------------
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

# Run Biomod, make sure to save results ---------------------------

#save on your desktop first
setwd(paste0("/Users/", Sys.info()[7], "/Desktop"))

# Computing the models
myBiomodModelOut = BIOMOD_Modeling(
  myBiomodData,
  models = c(
    # 'GLM',
    'GAM',
    # 'GBM',
    # 'CTA',
    # 'ANN',
    # 'SRE',
    # 'FDA',
    # 'MARS',
    'RF'),
    # 'MAXENT.Tsuruoka',
    # 'MAXENT.Phillips'),
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
  # 'GLM',
  'GAM',
  # 'GBM',
  # 'CTA',
  # 'ANN',
  # 'SRE',
  # 'FDA',
  # 'MARS',
  'RF')
  # 'MAXENT.Tsuruoka',
  # 'MAXENT.Phillips')

tss = cbind(models, tss)
colnames(tss) = c("SDM", "TSS_RUN1","TSS_RUN2","TSS_RUN3")
roc = cbind(models, roc)
colnames(roc) = c("SDM", "ROC_RUN1","ROC_RUN2","ROC_RUN3")
tss_roc = merge(tss, roc)
tss_roc = subset(tss_roc, SDM!="MAXENT.Tsuruoka")
tss_roc$SDM = gsub("MAXENT.Phillips", "MAXENT", tss_roc$SDM)
tss_roc

