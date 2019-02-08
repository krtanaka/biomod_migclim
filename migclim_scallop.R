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

dir = "/Users/Kisei" #pc

setwd(paste0(dir, "/Google Drive/R/Biomod/scallop/"))

analysis = "proj" #dispersion based on CM2.6 data, year 1-80

# Load POS file -----------------------------------------------------------
POS <- read_csv("POS.csv")
POS = POS[,c(84, 85, 1:80)]
colnames(POS)[1:2] = c("Y","X")
names(POS)
timestep = length(POS)-2

#----specifiy initial distirbution using survey data---- 
df = read_csv("Scallop_Catch_2014_2016.csv")

df = df[!is.na(df$lat) & !is.na(df$lon),]    

qplot(df$lon, df$lat, color=df$expcatnum)

#get survey lat lon
survey_latlon = as.data.frame(df[,c("lon", "lat")]) 

#get biomod lat lon
POS$biomod_grid_ID <- seq.int(nrow(POS))

biomod_latlon = POS[c("X", "Y", "biomod_grid_ID")]
coordinates(biomod_latlon) <- ~ X + Y
gridded(biomod_latlon) = T
biomod_latlon <- raster(biomod_latlon, layer = 1)

#attach depth and biomod rasters to survey resolution
biomod_grid = SDMTools::extract.data(survey_latlon, biomod_latlon)

df = cbind(df, biomod_grid)
qplot(df$lon, df$lat, color=df$biomod_grid) + scale_color_gradientn(colours = parula(length(df$biomod_grid)))

df = aggregate(expcatnum~biomod_grid, data=df, FUN=sum)
colnames(df) = c("biomod_grid_ID", "Obs_Catch")
POS = merge(POS, df, all = T)
qplot(POS$X, POS$Y, color=POS$Obs_Catch) + scale_color_gradientn(colours = parula(length(POS$Obs_Catch)))

POS$InitialDist = ifelse(POS$Obs_Catch > 0, 1, 0)
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

POS = POS[c(2:83, 85:86)]
names(POS)

setwd(paste0(dir, "/Google Drive/R/Biomod/scallop/migclim/", analysis))

dataframe2asc(POS)

# hsMap01 <- raster("hsmap1.asc")
# hsMap80 <- raster("hsmap_959.asc")
# 
# par(mfrow = c(1,2))
# plot(hsMap01);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
# plot(hsMap80);map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)


# do migclim
source(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim.plot_V2.R"))
source(paste0(dir, "/Google Drive/R/Biomod/lobster/migclim_migrate.R"))

# Small juveniles typically remain inshore and within a home range of about 5–15 km, and do not exhibit large- scale seasonal movements (Cooper et al., 1975). Mature individuals exhibit an average annual range of 32 km (Campbell, 1986), and have a higher tolerance to deeper and cooler waters.

barrierType = c("weak", "strong")
dispKernel = c(0.01, 0.05, 0.1)
iniMatAge = c(2,3,4) # reach sexually maturity at 2 years, but research suggests that egg production is low until they reach age 4 (NMFS 2002; Hart 2001)
lddFreq = 0
lddMinDist = 3
lddMaxDist = 4

setting = expand.grid(barrierType, dispKernel, iniMatAge, lddFreq, lddMinDist, lddMaxDist)
colnames(setting) = c("barrierType", "dispKernel", "iniMatAge", "lddFreq", "lddMinDist", "lddMaxDist")

for (i in 1:nrow(setting)) {
  
  setwd(paste0(dir, "/Google Drive/R/Biomod/scallop/migclim/proj"))
  
  temp_setting = setting[i,]
  
  par1 = as.character(temp_setting$barrierType)
  par2 = as.numeric(temp_setting$dispKernel)
  par3 = as.numeric(temp_setting$iniMatAge)
  par4 = as.numeric(temp_setting$lddFreq)
  par5 = as.numeric(temp_setting$lddMinDist)
  par6 = as.numeric(temp_setting$lddMaxDist)
  
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
    lddFreq = par4,
    lddMinDist = par5, 
    lddMaxDist = par6,
    # simulName = "MIGPOSTest", 
    replicateNb = 3, 
    overWrite = TRUE, 
    testMode = FALSE, 
    fullOutput = FALSE, 
    keepTempFiles = FALSE)
  
  setwd(paste0(dir, "/Google Drive/R/Biomod/scallop/migclim/proj/", par1, "_", par2, "_", par3, "_", par4, "_", par5, "_", par6))
  
  
  migclim.plot()
  
  
}

