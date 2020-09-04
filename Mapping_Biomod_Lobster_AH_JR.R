rm(list = ls())

library(sp)
library(automap)
library(gstat)
library(rworldxtra); data(countriesHigh)
library(readr)
library(easyGgplot2)
library(colorRamps)
library(grid)
library(gridExtra)
library(ggpubr)
library(Rmisc)

setwd("/Users/Kisei/Google Drive/R/Biomod/lobster")

source("/Users/Kisei/Google Drive/R/misc/color palette function.R")
steps = c("blue", "white", "red")
col = color.palette(steps, space="rgb")

# Mapping average occurance probability -------------------------------------------
get_avg = function(time_step){
  
  POS = read.csv(paste0("Biomod_1_80_", time_step, ".csv"))
  
  POS = POS[,c(2,1,3:82)]
  
  names(POS)[1] = "Y"
  names(POS)[2] = "X"
  
  names(POS)
  
  avg = cbind(POS[,1:2], POS$hsmap1*0.001)
  colnames(avg) = c("latitude", "longitude", "beta")
  
  coordinates(avg) = ~longitude + latitude
  
  auto = autofitVariogram(beta ~ 1, avg)
  g = gstat(formula = beta ~ 1, model = auto$var_model, data = avg, maxdist = 0.05)
  
  xrange = range(avg$longitude)
  yrange = range(avg$latitude)
  grid= expand.grid(longitude = seq(from = xrange[1], to = xrange[2], by = .03), 
                    latitude = seq(from = yrange[1], to = yrange[2], by = .03))
  gridded(grid) = ~longitude + latitude
  
  p = predict(g, newdata = grid)
  
  return(p)
  
}
spring = get_avg("spring")
fall = get_avg("fall")

spring$Spring = spring$var1.pred
spring$Fall = fall$var1.pred

spplot(spring, 
       main=list(label="American lobster - fall",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:100)/100, #for slope
       col.regions = matlab.like(100),
       zcol = c("Fall"),
       scales=list(draw=T),
       colorkey = F) 

spplot(spring, 
       main=list(label="American lobster - spring",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:100)/100, #for slope
       col.regions = matlab.like(100),
       zcol = c("Spring"),
       scales=list(draw=T),
       colorkey = F) 

# Mapping trend ---------------------------------------------------------------
# time_step = "sep" #choose individual month or season (e.g. apr, sep, spring, fall)
get_trend = function(time_step){
  
  POS = read.csv(paste0("Biomod_1_80_", time_step, ".csv"))
  
  POS[,c(3:82)] = POS[,c(3:82)]*0.001
  
  POS = POS[,c(2,1,3:82)]
  
  names(POS)[1] = "Y"
  names(POS)[2] = "X"
  
  names(POS)
  
  betaf = function(vec){
    
    beta = lm(vec ~ seq(1:80))$coef[2]
    p = summary(lm(vec ~ seq(1:80)))$ coefficients [2,4]
    return(beta) # beta gives you a slope, if you want p-value, change it to p
    #   return(p) # beta gives you a slope, if you want p-value, change it to p
    
  }
  
  res = as.data.frame(apply(POS[, 3:82], 1, betaf))
  trend = cbind(POS[,1:2], res)
  colnames(trend) = c("latitude", "longitude", "beta")
  
  coordinates(trend) = ~longitude + latitude
  
  auto = autofitVariogram(beta ~ 1, trend)
  g = gstat(formula = beta ~ 1, model = auto$var_model, data = trend, maxdist = 0.05)
  
  xrange = range(trend$longitude)
  yrange = range(trend$latitude)
  grid= expand.grid(longitude = seq(from = xrange[1], to = xrange[2], by = .03), 
                    latitude = seq(from = yrange[1], to = yrange[2], by = .03))
  gridded(grid) = ~longitude + latitude
  
  p = predict(g, newdata = grid)
  
  return(p)
  
}
spring = get_trend("spring")
fall = get_trend("fall")

spring$Spring = spring$var1.pred
spring$Fall = fall$var1.pred

max = max(spring@data, na.rm = T)*12000
min = max*-1

spplot(spring, 
       main=list(label="American lobster",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/10000, #for slope
       col.regions = col,
       zcol = c("Spring", "Fall"),
       scales=list(draw=T),
       colorkey = T) 
