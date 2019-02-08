library(sp)
library(automap)
library(gstat)
library(rworldxtra); data(countriesHigh)

# see changes in HSI ------------------------------------------------------
setwd("C:/Users/Kisei/Google Drive/R/Biomod/lobster")
setwd("~/Kisei/Google Drive/R/Biomod/lobster")

# time_step = "sep" #choose individual month or season (e.g. apr, sep, spring, fall)
get_trend = function(time_step){
  
  POS = read.csv(paste0("C:/Users/Kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_", time_step, ".csv"))
  # POS = read.csv(paste0("~/Google Drive/R/Biomod/Biomod_1_80_", month, ".csv"))
  
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

# monthly results ---------------------------------------------------------
apr = get_trend("apr")
may = get_trend("may")
jun = get_trend("jun")
sep = get_trend("sep")
oct = get_trend("oct")
nov = get_trend("nov")

apr$April = apr$var1.pred
apr$May = may$var1.pred
apr$June = jun$var1.pred
apr$September = sep$var1.pred
apr$October = oct$var1.pred
apr$November = nov$var1.pred

source("C:/Users/Kisei/Google Drive/R/misc/color palette function.R")
source("~/Google Drive/R/misc/color palette function.R")

steps = c("blue", "white", "red")
pal = color.palette(steps, space="rgb")
col = pal

max = max(apr@data, na.rm = T)*100
min = max*-1

# png("C:/Users/Kisei/Google Drive/R/Biomod/Habitat_Change.png", width = 4000, height = 5000, res = 500)
png("Habitat_Change_Month.png", width = 4000, height = 5000, res = 500)
spplot(apr, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/100, #for slope
       col.regions = col,
       zcol = c("April", "May", "June", "September", "October", "November"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()

# seasonal results --------------------------------------------------------
spring = get_trend("spring")
fall = get_trend("fall")

spring$Spring = spring$var1.pred
spring$Fall = fall$var1.pred

source("C:/Users/Kisei/Google Drive/R/misc/color palette function.R")
source("~/Google Drive/R/misc/color palette function.R")

steps = c("blue", "white", "red")
pal = color.palette(steps, space="rgb")
col = pal

max = max(spring@data, na.rm = T)*100
min = max*-1

# png("C:/Users/Kisei/Google Drive/R/Biomod/Habitat_Change.png", width = 4000, height = 5000, res = 500)
png("Habitat_Change_Season.png", width = 4000, height = 5000, res = 500)
spplot(spring, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/100, #for slope
       col.regions = col,
       zcol = c("Spring", "Fall"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()


