library(MigClim)
library(raster)
library(SDMTools)
library(maps)
library(gstat)
library(automap)
library(ggplot2)
library(pals)
library(readr)
library(rworldxtra); data("countriesHigh")

#MigClim.userGuide()

dir = "/Users/Kisei"

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
timestep = length(POS)-2

# combine months ----------------------------------------------------------

latlon =  read.csv(paste0("Biomod_1_80_", "apr", ".csv"))[1:2]

d = as.data.frame(matrix(nrow = 10497, ncol = 80))
d[is.na(d)] <- 0

for (i in 1:80) {
  
  POS_1 = read_csv(paste0("Biomod_1_80_", "apr", ".csv"))[3:82]
  POS_2 = read_csv(paste0("Biomod_1_80_", "may", ".csv"))[3:82]
  POS_3 = read_csv(paste0("Biomod_1_80_", "jun", ".csv"))[3:82]
  
  POS_1 = POS_1[i]
  POS_2 = POS_2[i]
  POS_3 = POS_3[i]
  
  POS_Spring = cbind(POS_1, POS_2, POS_3)
  POS_Spring = as.data.frame(round(rowMeans(POS_Spring),0))
  
  d[i]=POS_Spring
  
}

POS_Spring = d

d = as.data.frame(matrix(nrow = 10497, ncol = 80))
d[is.na(d)] <- 0

for (i in 1:80) {
  
  POS_1 = read_csv(paste0("Biomod_1_80_", "sep", ".csv"))[3:82]
  POS_2 = read_csv(paste0("Biomod_1_80_", "oct", ".csv"))[3:82]
  POS_3 = read_csv(paste0("Biomod_1_80_", "nov", ".csv"))[3:82]
  
  POS_1 = POS_1[i]
  POS_2 = POS_2[i]
  POS_3 = POS_3[i]
  
  POS_Fall = cbind(POS_1, POS_2, POS_3)
  POS_Fall = as.data.frame(round(rowMeans(POS_Fall),0))
  
  d[i]=POS_Fall
  
}

POS_Fall = d

rm(POS_1, POS_2, POS_3)

POS_Spring = as.data.frame(t(POS_Spring)); POS_Spring$Index = seq(1,160,2)
POS_Fall = as.data.frame(t(POS_Fall)); POS_Fall$Index = seq(2,160,2)

POS = rbind(POS_Spring, POS_Fall)
POS <- POS[order(POS$Index),] 
POS = data.frame(t(POS))

rm(POS_Spring, POS_Fall)

colnames(POS) = POS[10498, ] # the last row will be the header
POS = POS[-10498, ]     

POS = cbind(latlon, POS); rm(latlon)
colnames(POS)[-c(1:2)] = paste0("hsmap", colnames(POS[-c(1:2)]))

POS = POS[,c(2,1,3:length(POS))]

names(POS)[1] = "Y"
names(POS)[2] = "X"

names(POS)
timestep = length(POS)-2

#----specifiy initial distirbution---- 
if (analysis == "proj") distBound = 500
if (analysis == "hind") distBound = 500

POS$InitialDist = rep(NA, nrow(POS))

for (i in 1:nrow(POS)){
  
  x = median(as.numeric(POS[i,3:4])) #chose 1st time step (3) or average of first 3 time steps (3:5); 3:8 means first year apr-nov
  
  if(x > distBound){
    
    y = rbinom(n=1, size=1, prob=(mean(x)/1000))
    
  } else {
    
    y = 0
  }
  
  POS[i,"InitialDist"] = y
  
}

summary(POS$InitialDist); qplot(POS$X, POS$Y, color = POS$InitialDist)

#set barrier
if (analysis == "proj") barrier = 200
if (analysis == "hind") barrier = 200

POS$Barrier = rep(NA, nrow(POS))

for (i in 1:nrow(POS)){
  
  x = median(as.numeric(POS[i,3:4]))
  
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
# hsMap80 <- raster("hsmap160.asc")
# 
# par(mfrow = c(1,2))
# plot(hsMap01);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
# plot(hsMap80);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)

source(paste0(dir, "/Google Drive/R/Biomod/migclim.plot.R"))

#----specifiy initial distirbution using survey data---- 
ct = read.csv(paste0(dir, "/Google Drive/Research/Geostatistical_Model/VAST/data/LIS.csv"), header = T)
ri = read.csv(paste0(dir, "/Google Drive/Research/Geostatistical_Model/VAST/data/RI.csv"), header = T)
ma = read.csv(paste0(dir, "/Google Drive/Research/Geostatistical_Model/VAST/data/MA.csv"), header = T)
me = read.csv(paste0(dir, "/Google Drive/Research/Geostatistical_Model/VAST/data/ME.csv"), header = T)
fe = read.csv(paste0(dir, "/Google Drive/Research/Geostatistical_Model/VAST/data/NEFSC.csv"), header = T)

df = rbind(ct, ri, ma, me, fe); rm(ct, ri, ma, me, fe)
df = df[!is.na(df$Lat) & !is.na(df$Lon),]    

df  = df[which(df$Year %in% c(2014:2016)),]
# df$Catch_n = ifelse(df$Catch_n > 1, 1, 0) #try chaniging this to 5 to make it less sensitive

df  = df[which(df$Stat_Area %in% c(464:465, 511:515, 521:526, 537:539, 541:543, 561:562, 551:552, 611:616, 621:626)),]

qplot(df$Lon, df$Lat, color=log10(df$Catch_n))

#get survey lat lon
survey_latlon = as.data.frame(df[,c("Lon", "Lat")]) 

#get biomod lat lon
POS$biomod_grid_ID <- seq.int(nrow(POS))

biomod_latlon = POS[c("X", "Y", "biomod_grid_ID")]
coordinates(biomod_latlon) <- ~ X + Y
gridded(biomod_latlon) = T
biomod_latlon <- raster(biomod_latlon, layer = 1)

#attach depth and biomod rasters to survey resolution
biomod_grid = SDMTools::extract.data(survey_latlon, biomod_latlon)

df = cbind(df, biomod_grid)
qplot(df$Lon, df$Lat, color=df$biomod_grid) + scale_color_gradientn(colours = parula(length(df$biomod_grid)))

df = aggregate(Catch_n~biomod_grid, data=df, FUN=sum)
colnames(df) = c("biomod_grid_ID", "Obs_Catch")
POS = merge(POS, df, all = T)
qplot(POS$X, POS$Y, color=log(POS$Obs_Catch)) + scale_color_gradientn(colours = parula(length(POS$Obs_Catch)))

POS$InitialDist = ifelse(POS$Obs_Catch > 3, 1, 0)
POS$InitialDist = ifelse(is.na(POS$InitialDist) == T, 0, POS$InitialDist)

summary(POS$InitialDist); qplot(POS$X, POS$Y, color = POS$InitialDist)

#set barrier
barrier = 100

POS$Barrier = rep(NA, nrow(POS))

for (i in 1:nrow(POS)){
  
  x = median(as.numeric(POS[i,3:4]))
  
  if(x > barrier){
    
    y = 0
    
  } else {
    
    y = 1
    # y = rbinom(n=1, size=1, prob=(mean(x)/1000))
  }
  
  POS[i,"Barrier"] = y
  
} 

summary(POS$Barrier); qplot(POS$X, POS$Y, color = POS$Barrier)
rm(biomod_grid, biomod_latlon, survey_latlon, df, d, x, y)

POS = POS[c(2:163, 165:166)]
names(POS)

setwd(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim/", analysis))

dataframe2asc(POS)

# hsMap01 <- raster("hsmap1.asc")
# hsMap80 <- raster("hsmap_959.asc")
# 
# par(mfrow = c(1,2))
# plot(hsMap01);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
# plot(hsMap80);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)


# do migclim --------------------------------------------------------------

if (analysis == "proj"){
  
  source(paste0(dir, "/Google Drive/R/Biomod/migclim.plot_V2.R"))
  source(paste0(dir, "/Google Drive/R/Biomod/migclim_migrate.R"))
  
  barrierType = c("weak", "strong")[1]
  dispKernel = c(0.01, 0.05, 0.1)[1]
  iniMatAge = c(4,5,6,7,8)[1:2]
  # lddFreq = 0
  # lddMinDist = 3
  # lddMaxDist = 4
  
  setting = expand.grid(barrierType, dispKernel, iniMatAge)
  colnames(setting) = c("barrierType", "dispKernel", "iniMatAge")
  
  # setting = expand.grid(barrierType, dispKernel, iniMatAge, lddFreq, lddMinDist, lddMaxDist)
  # colnames(setting) = c("barrierType", "dispKernel", "iniMatAge", "lddFreq", "lddMinDist", "lddMaxDist")
  
  for (i in 1:nrow(setting)) {
    
    setwd(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim/proj"))
    
    temp_setting = setting[i,]
    
    par1 = as.character(temp_setting$barrierType)
    par2 = as.numeric(temp_setting$dispKernel)
    par3 = as.numeric(temp_setting$iniMatAge)
    # par4 = as.numeric(temp_setting$lddFreq)
    # par5 = as.numeric(temp_setting$lddMinDist)
    # par6 = as.numeric(temp_setting$lddMaxDist)
    
    MigClim.migrate_KRT(
      iniDist = POS[,c("X","Y","InitialDist")],
      hsMap = POS[, 3:c(timestep+2)], 
      rcThreshold = 0, 
      barrier = POS[, "Barrier"], 
      barrierType = par1, 
      envChgSteps = timestep, 
      dispSteps = 2, 
      dispKernel = par2, #big infleunce by changing to 0.1
      iniMatAge = par3, 
      propaguleProd = c(0.02, 0.1, 0.5, 0.9), 
      # lddFreq = par4,
      # lddMinDist = par5, 
      # lddMaxDist = par6,
      # simulName = "MIGPOSTest", 
      replicateNb = 3, 
      overWrite = TRUE, 
      testMode = FALSE, 
      fullOutput = FALSE, 
      keepTempFiles = FALSE)
    
    # setwd(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim/proj/", par1, "_", par2, "_", par3, "_", par4, "_", par5, "_", par6))
    
    setwd(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim/proj/", par1, "_", par2, "_", par3))
    
    migclim.plot()
    
    
  }
  
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
}