library(MigClim)
library(raster)
library(SDMTools)
library(maps)
library(gstat)
library(automap)
library(ggplot2)
library(pals)
library(rworldxtra); data("countriesHigh")

#MigClim.userGuide()

dir = "C:/Users/Kisei" #pc
# dir = "~" #mac

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))

analysis = "proj" #dispersion based on CM2.6 data, year 1-80
# analysis = "hind" #dispersion based on fvcom data, 1984-2013

# month specific ----------------------------------------------------------
month = "apr"

if (analysis == "proj"){
  POS = read.csv(paste0("Biomod_1_80_", month, ".csv"))[1:82] #don't need last two columns
}

if (analysis == "hind"){
  POS = read.csv(paste0("Biomod_1_30_", month, ".csv"))[1:82] #don't need last two columns
}

POS = POS[,c(2,1,3:length(POS))]

names(POS)[1] = "Y"
names(POS)[2] = "X"

names(POS)

# combine months ----------------------------------------------------------

latlon =  read.csv(paste0("Biomod_1_80_", "apr", ".csv"))[1:2]
 
POS_1 = read.csv(paste0("Biomod_1_80_", "apr", ".csv"))[3:82]; POS_1 = as.data.frame(t(POS_1)); POS_1$Index = seq(4,960,12)
POS_2 = read.csv(paste0("Biomod_1_80_", "may", ".csv"))[3:82]; POS_2 = as.data.frame(t(POS_2)); POS_2$Index = seq(5,960,12)
POS_3 = read.csv(paste0("Biomod_1_80_", "jun", ".csv"))[3:82]; POS_3 = as.data.frame(t(POS_3)); POS_3$Index = seq(6,960,12)
POS_4 = read.csv(paste0("Biomod_1_80_", "sep", ".csv"))[3:82]; POS_4 = as.data.frame(t(POS_4)); POS_4$Index = seq(9,960,12)
POS_5 = read.csv(paste0("Biomod_1_80_", "oct", ".csv"))[3:82]; POS_5 = as.data.frame(t(POS_5)); POS_5$Index = seq(10,960,12)
POS_6 = read.csv(paste0("Biomod_1_80_", "nov", ".csv"))[3:82]; POS_6 = as.data.frame(t(POS_6)); POS_6$Index = seq(11,960,12)

POS = rbind(POS_1, POS_2, POS_3, POS_4, POS_5, POS_6)
POS <- POS[order(POS$Index),] 
POS = data.frame(t(POS))

colnames(POS) = POS[10498, ] # the last row will be the header
POS = POS[-10498, ]     

POS = cbind(latlon, POS)
colnames(POS)[-c(1:2)] = paste0("hsmap_", colnames(POS[-c(1:2)]))

POS = POS[,c(2,1,3:length(POS))]

names(POS)[1] = "Y"
names(POS)[2] = "X"

names(POS); rm(POS_1, POS_2, POS_3, POS_4, POS_5, POS_6)

#----specifiy initial distirbution---- 
if (analysis == "proj") distBound = 500
if (analysis == "hind") distBound = 500

POS$InitialDist = rep(NA, nrow(POS))

for (i in 1:nrow(POS)){
  
  x = median(as.numeric(POS[i,3:8])) #chose 1st time step (3) or average of first 3 time steps (3:5); 3:8 means first year apr-nov
  
  if(x > distBound){
    
    y = rbinom(n=1, size=1, prob=(mean(x)/1000))
    
  } else {
    
    y = 0
  }
  
  POS[i,"InitialDist"] = y
  
}

summary(POS$InitialDist); qplot(POS$X, POS$Y, color = POS$InitialDist)

#set barrier
if (analysis == "proj") barrier = 300
if (analysis == "hind") barrier = 300

POS$Barrier = rep(NA, nrow(POS))

for (i in 1:nrow(POS)){
  
  x = median(as.numeric(POS[i,3:8]))
  
  if (analysis == "proj") {
    
    if(x > barrier){
      
      y = 0
      
    } else {
      
      y = 1
      # y = rbinom(n=1, size=1, prob=(mean(x)/1000))
    }
    
    POS[i,"Barrier"] = y
    
  }
  
  if (analysis == "hind") {
    
    if(x > barrier){
      
      y = 0
      
    } else {
      
      y = 1
      # y = rbinom(n=1, size=1, prob=(mean(x)/1000))
    }
    
    POS[i,"Barrier"] = y
    
  }
  
} 

summary(POS$Barrier); qplot(POS$X, POS$Y, color = POS$Barrier)

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim/", analysis))

dataframe2asc(POS)

# hsMap01 <- raster("hsmap1.asc")
# hsMap80 <- raster("hsmap_959.asc")
# 
# par(mfrow = c(1,2))
# plot(hsMap01);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
# plot(hsMap80);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)

source(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim.plot.R"))

if (analysis == "proj"){
  
  timestep = 80

  MigClim.migrate(
    iniDist = POS[,c("X","Y","InitialDist")],
    hsMap = POS[, 3:c(timestep+2)], 
    rcThreshold = 0, 
    barrier = POS[, "Barrier"], 
    barrierType = "weak", 
    envChgSteps = timestep, 
    dispSteps = 2, 
    dispKernel = c(0.01), #big infleunce by changing to 0.1
    iniMatAge = 5, 
    propaguleProd = c(0.02, 0.1, 0.5, 0.9), 
    lddFreq = 0,
    lddMinDist = 3, 
    lddMaxDist = 4,
    simulName = "MIGPOSTest", 
    replicateNb = 3, 
    overWrite = TRUE, 
    testMode = FALSE, 
    fullOutput = FALSE, 
    keepTempFiles = FALSE)
  
  migclim.plot(asciiFile="MIGPOSTest/MIGPOSTest3_raster.asc",outDir="", fileFormat="jpeg", fullOutput=FALSE)
  
}

if (analysis == "hind"){
  
  MigClim.migrate(
    iniDist = POS[,c("X","Y","InitialDist")],
    hsMap = POS[, 3:32], 
    rcThreshold = 500, 
    barrier = POS[, "Barrier"], 
    barrierType = "weak", 
    envChgSteps = length(POS[, 3:32]), 
    dispSteps = 2, 
    dispKernel = c(1.0,0.4,0.16,0.06,0.03), 
    iniMatAge = 4, 
    propaguleProd = c(0.1), 
    lddFreq = 0,
    lddMinDist = 3, 
    lddMaxDist = 4,
    simulName = "MIGPOSTest", 
    replicateNb = 3, 
    overWrite = TRUE, 
    testMode = FALSE, 
    fullOutput = FALSE, 
    keepTempFiles = FALSE)
  
  migclim.plot(asciiFile="MIGPOSTest/MIGPOSTest3_raster.asc",outDir="", fileFormat="jpeg", fullOutput=FALSE)
  
}