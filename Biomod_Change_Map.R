library(sp)
library(automap)
library(gstat)
library(rworldxtra); data(countriesHigh)

# see changes in HSI
setwd("/Kisei/Google Drive/R/Biomod/lobster")

# time_step = "sep" #choose individual month or season (e.g. apr, sep, spring, fall)
get_trend = function(time_step){
  
  POS = read.csv(paste0("/Users/Kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_", time_step, ".csv"))

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

source("/Users/Kisei/Google Drive/R/misc/color palette function.R")

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

source("/Users/Kisei/Google Drive/R/misc/color palette function.R")

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

# ks test -----------------------------------------------------------------

#load biomod output (cm2.6 1-80 years presence-absence)
setwd("/Users/Kisei/Google Drive/R/Biomod/lobster")
df <- read_csv("Biomod_1_80_fall.csv") #lobster fall biomod output
df <- read_csv("Biomod_1_80_spring.csv") #lobster spring biomod output

#choose analysis
test = "first_last_10_years" #compare first and last 10 years
test = "first_last_40_years" #compare first and last 40 years

latlon = df[1:2]
df = df[3:82]

res = matrix(0, ncol = 4, nrow = nrow(df))
colnames(res) = c("D", "ks.ts", "ks.lt", "ks.gt")

for (i in 1:10497) {
  
  d = data.frame(df[i,1:80])
  
  if (test == "first_last_10_years") {
    x = as.numeric(d[,1:10]*0.001) #first 10 years
    y = as.numeric(d[,71:length(d)]*0.001) #last 10 years
  }
  
  if (test == "first_last_40_years") {
    x = as.numeric(d[,1:40]*0.001) #1st 315 months sep 1981- dec 2007
    y = as.numeric(d[,41:length(d)]*0.001) #last 120 months
  }
  
  # P_1 = ecdf(x)
  # P_2 = ecdf(y)# P is a function giving the empirical CDF of X
  # 
  # plot(P_1, col = 4)
  # lines(P_2, col = 2)
  
  ks.ts <- ks.test(x, y, alternative = "two.sided") #two-sided (equal)
  ks.lt <- ks.test(x, y, alternative = "less") #the CDF of x lies below that of y
  ks.gt <- ks.test(x, y, alternative = "greater") # the CDF of x lies above and hence to the left of that for y
  
  D = ks.ts$statistic
  ts = ks.ts$p.value
  lt = ks.lt$p.value
  gt = ks.gt$p.value
  
  # overlap <- function(x, y) {
  #   F.x <- ecdf(x); F.y <- ecdf(y)
  #   z <- uniroot(function(z) F.x(z) + F.y(z) - 1, interval<-c(min(c(x,y)), max(c(x,y))))
  #   return(list(Root=z, F.x=F.x, F.y=F.y))
  # }
  # 
  # d = overlap(x,y)
  # 
  # P_1(0.5) 
  # P_2(0.5) 
  
  # #make cdf comparison plots
  # x = data.frame(x)
  # y = data.frame(y)
  # 
  # if (test == "first_last_10_years") {
  # 
  #   x$Period = "First_10_years"
  #   y$Period = "Last_10_years"
  # 
  # }
  # 
  # if (test == "1981_2007_Baseline") {
  # 
  #   x$Period = "1981-2007"
  #   y$Period = "2008-2017"
  # 
  # }
  # 
  # colnames(x) = c("value", "period")
  # colnames(y) = c("value", "period")
  # 
  # xy = rbind(x,y)
  # 
  # ggplot2.histogram(data=xy, xName='value',
  #                   groupName='period',
  #                   legendPosition = "top",
  #                   alpha = 0.5,
  #                   addDensity=TRUE,
  #                   addMeanLine=TRUE,  meanLineSize=1.5) +
  #   ylim(0,100) + xlim(0.15,0.3) + xlab("Presence Probability") +
  #   theme_classic()
  
  res[i,] = c(D, ts, lt, gt)
  
}

df = cbind(latlon, res)

D = df[,c("x","y","D")] #D statistic
ks_ts = df[,c("x","y","ks.ts")] #two-tail
ks_lt = df[,c("x","y","ks.lt")] #less than
ks_gt = df[,c("x","y","ks.gt")] #greater than

sp_map = function(x, y, z, df){
  
  coordinates(df) = ~x + y
  
  auto = autofitVariogram(z ~ 1, df)
  g = gstat(formula = z ~ 1, model = auto$var_model, data = df, maxdist = 0.05)
  
  xrange = range(df$x)
  yrange = range(df$y)
  grid= expand.grid(x = seq(from = xrange[1], to = xrange[2], by = 0.25), 
                    y = seq(from = yrange[1], to = yrange[2], by = 0.25))
  
  gridded(grid) = ~x + y
  
  p = predict(g, newdata = grid)
  
  return(p)
}

D = sp_map(D$x, D$y, D$D, D); plot(D)
ks_ts = sp_map(ks_ts$x, ks_ts$y, ks_ts$ks.ts, ks_ts); plot(ks_ts)
ks_lt = sp_map(ks_lt$x, ks_lt$y, ks_lt$ks.lt, ks_lt); plot(ks_lt)
ks_gt = sp_map(ks_gt$x, ks_gt$y, ks_gt$ks.gt, ks_gt); plot(ks_gt)

D$D = D$var1.pred
D$Two_Sample_KS_Pval = ks_ts$var1.pred
D$P_0.05 = ifelse(ks_ts$var1.pred < 0.05, ks_ts$var1.pred, NA)
D$KS_Less = ks_lt$var1.pred
D$KS_Greater = ks_gt$var1.pred

jpeg(paste0("Two-Sample_KS_test_", test, "_", Sys.Date(), ".jpg"), res = 500, height = 5, width = 7, units = "in")
spplot(D, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:1000)/1000,
       col.regions = parula(1000),
       zcol = c("Two_Sample_KS_Pval"),
       scales=list(draw=T),
       colorkey = T)
dev.off()

jpeg(paste0("Two-Sample_KS_test_P0.05_", test, "_", Sys.Date(), ".jpg"), res = 500, height = 5, width = 7, units = "in")
spplot(D, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:50)/1000,
       col.regions = parula(1000),
       zcol = c("P_0.05"),
       scales=list(draw=T),
       colorkey = T)
dev.off()

jpeg(paste0("Two-sample_KS_Tests_D_", test, "_", Sys.Date(), ".jpg"), res = 500, height = 5, width = 7, units = "in")
spplot(D, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:35)/100,
       col.regions = parula(length(na.omit(D$var1.pred))),
       zcol = c("D"),
       scales=list(draw=T),
       colorkey = T)
dev.off()

# The significance of difference between two REC curves can be assessed by examining the maximum deviation between the two curves across all values of E. This corresponds to the Kolmogorov-Smirnov (KS) two-sample test (DeGroot, 1986) for judging the hypothesis that the error e generated by two models f and g follows the same distribution. The KS-test requires no assumptions about the underlying CDF of z. 