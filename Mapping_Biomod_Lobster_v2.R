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
# setwd("/Users/Kisei/Desktop/lobster")

source("/Users/Kisei/Google Drive/R/misc/color palette function.R")
steps = c("blue", "white", "red")
col = color.palette(steps, space="rgb")

# map average probability first -------------------------------------------
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

png("/Users/Kisei/Desktop/Lobster_Fall.png", width = 5, height = 4, res = 500, units = "in")
spplot(spring, 
       main=list(label="American lobster - fall",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:100)/100, #for slope
       col.regions = matlab.like(100),
       zcol = c("Fall"),
       scales=list(draw=T),
       colorkey = F) 
dev.off()

png("/Users/Kisei/Desktop/Lobster_Spring.png", width = 5, height = 4, res = 500, units = "in")
spplot(spring, 
       main=list(label="American lobster - spring",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:100)/100, #for slope
       col.regions = matlab.like(100),
       zcol = c("Spring"),
       scales=list(draw=T),
       colorkey = F) 
dev.off()

# Map trend ---------------------------------------------------------------
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

max = max(spring@data, na.rm = T)*12000
min = max*-1

# png("C:/Users/Kisei/Google Drive/R/Biomod/Habitat_Change.png", width = 4000, height = 5000, res = 500)
png("/Users/Kisei/Desktop/Lobster_Habitat_Change_Season.png", width = 2000, height = 2200, res = 500)
spplot(spring, 
       main=list(label="American lobster",cex=1.5),
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/10000, #for slope
       col.regions = col,
       zcol = c("Spring", "Fall"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()

# ks test -----------------------------------------------------------------
#load biomod output (cm2.6 1-80 years presence-absence)
season = list("fall", "spring")

#choose analysis
test = "first_last_10_years" #compare first and last 10 years
# test = "first_last_20_years" #compare first and last 40 years

for (j in 1:length(season)){
  
  # j = 2
  
  if (j == 1) df <- read_csv("Biomod_1_80_fall.csv") #lobster fall biomod output
  if (j == 2) df <- read_csv("Biomod_1_80_spring.csv") #lobster spring biomod output
  
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
    #   # ylim(0,100) + xlim(0.15,0.3) + 
    #   xlab("Habitat Suitability") +
    #   theme_classic()
    
    res[i,] = c(D, ts, lt, gt)
    
  }
  
  ks = cbind(latlon, res)
  
  #add GOMGBK and SNE area
  names(ks)
  latlon = ks[,c(1,2)]; plot(latlon)
  library(sp)
  library(maptools)
  library(colorRamps)
  coordinates(latlon)=~x+y
  area = rgdal::readOGR("/Users/Kisei/Google Drive/Research/GIS/NOAA_Statistical_Area/Statistical_Areas.shp")
  area = area[which(area$Id %in% c(464:465, 511:515, 521:526, 537:539, 541:543, 561:562, 611:616, 621:626)),]
  CRS.new = CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") #EPSG:102003
  proj4string(latlon) <- CRS.new 
  proj4string(area) <- CRS.new
  area <- over(latlon, area)
  colnames(area)[1] = "GOMGBK_SNE"
  ks = cbind(ks, area)
  ks$GOMGBK_SNE  = ifelse(ks$GOMGBK_SNE %in% c(464:465, 511:515, 521:526, 561:562), "GOM_GBK", 
                          ifelse(ks$GOMGBK_SNE %in% c(537:539, 541:543, 561:562, 611:616, 621:626), "SNE", "NA"))
  # qplot(ks$x, ks$y, colour = ks$GOMGBK_SNE)
  
  #add nearshore areas
  area = rgdal::readOGR("/Users/Kisei/Google Drive/Research/GIS/Lobster_Management_Areas/Lobster_Management_Areas.shp")
  area = area[which(area$AREANAME %in% c("EEZ Nearshore Management Area 1",
                                         "EEZ Nearshore Management Area 2",
                                         "EEZ Nearshore Management Area 4",
                                         "EEZ Nearshore Outer Cape Lobster Management Area")),]
  proj4string(latlon) <- CRS.new 
  proj4string(area) <- CRS.new
  area <- over(latlon, area)
  colnames(area)[1] = "Nearshore_Area"
  ks = cbind(ks, area)
  ks$Nearshore_Area <- ifelse(ks$Nearshore_Area %in% c("EEZ Nearshore Management Area 1",
                                                       "EEZ Nearshore Outer Cape Lobster Management Area"), 
                              "GOM_GBK Neashore Management Area", 
                              ifelse(ks$Nearshore_Area %in% c("EEZ Nearshore Management Area 2", 
                                                              "EEZ Nearshore Management Area 4"), 
                                     "SNE Neashore Management Area",
                                     "NA"))
  # qplot(ks$x, ks$y, colour = ks$Nearshore_Area)
  
  D = ks[,c("x","y","D", "Nearshore_Area", "GOMGBK_SNE")] #D statistic
  ks_ts = ks[,c("x","y","ks.ts", "Nearshore_Area", "GOMGBK_SNE")] #two-tail
  ks_lt = ks[,c("x","y","ks.lt", "Nearshore_Area", "GOMGBK_SNE")] #less than
  ks_gt = ks[,c("x","y","ks.gt", "Nearshore_Area", "GOMGBK_SNE")] #greater than
  
  xlims <- range(pretty(ks$x));ylims <- range(pretty(ks$y))
  
  ks$P_0.05 = ifelse(ks$ks.ts < 0.05, ks$ks.ts, NA)
  
  ks1 = ks[,c("x","y","GOMGBK_SNE","D", "P_0.05")]; colnames(ks1)[3] = "Area"
  ks2 = ks[,c("x","y","Nearshore_Area","D", "P_0.05")]; colnames(ks2)[3] = "Area"
  ks_df = rbind(ks1, ks2)
  
  ks_df$Area = gsub('GBK', 'GB', ks_df$Area)
  ks_df$Area <- gsub(' Neashore Management Area', ' Neashore', ks_df$Area)
  
  scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
    xbreaks <- seq(xmin,xmax,step)
    xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
    return(scale_x_continuous("", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
  }
  scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
    ybreaks <- seq(ymin,ymax,step)
    ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
    return(scale_y_continuous("", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
  }  
  
  # p1 = ggplot() +
  #   geom_point(data = ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], 
  #              aes(x = x, y = y, color = D),
  #              size = 1) +
  #   # xlab("") +
  #   # ylab("") +
  #   facet_wrap(~Area, scales = "free", ncol = 1) + 
  #   # borders(xlim = xlims,ylim = ylims, fill = "gray") +
  #   # coord_quickmap(xlim = xlims,ylim = ylims) +
  #   coord_quickmap() +
  #   scale_colour_gradientn(colours = col(length(unique(ks_df$D)))) + #parura(100)
  #   theme_pubr() + 
  #   theme(legend.position="right") + 
  #   # theme(legend.justification = c(0,1), legend.position = c(0,1))
  #   scale_x_longitude(xmin=-180, xmax=180, step=2) +
  #   scale_y_latitude(ymin=-180, ymax=180, step=2) +
  #   ggtitle(paste0("American_lobster_", season[j]))

  # p2 = ggplot() +
  #   geom_point(data = ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], 
  #              aes(x = x, y = y, color = P_0.05),
  #              size = 1) +
  #   # xlab("") +
  #   # ylab("") +
  #   facet_wrap(~Area, scales = "free", ncol = 1) + 
  #   # borders(xlim = xlims,ylim = ylims, fill = "gray") +
  #   # coord_quickmap(xlim = xlims,ylim = ylims) +
  #   coord_quickmap() +
  #   scale_colour_gradientn(colours = col(length(unique(ks_df$P_0.05)))) + #parura(100)
  #   theme_pubr() + 
  #   theme(legend.position="right") + 
  #   # theme(legend.justification = c(0,1), legend.position = c(0,1))
  #   scale_x_longitude(xmin=-180, xmax=180, step=2) +
  #   scale_y_latitude(ymin=-180, ymax=180, step=2) +
  #   ggtitle(paste0("American_lobster_", season[j]))
  
  p1 = ggplot(data = ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], 
              aes(x = x, y = y, fill = D)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colours = col(length(unique(ks_df$D)))) + #parura(100)
    facet_wrap(~Area, ncol = 1) + 
    coord_quickmap(xlim = xlims,ylim = ylims) +
    borders(xlim = xlims,ylim = ylims, fill = "gray") +
    theme_pubr(I(18)) + 
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_x_longitude(xmin=-180, xmax=180, step=2) +
    scale_y_latitude(ymin=-180, ymax=180, step=2) +
    ggtitle(paste0("American_lobster_", season[j]))
  
  p2 = ggplot(data = ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], 
              aes(x = x, y = y, fill = P_0.05)) +
    geom_raster(interpolate = F) +
    scale_fill_gradientn(colours = col(length(unique(ks_df$P_0.05)))) + #parura(100)
    facet_wrap(~Area, ncol = 1) + 
    coord_quickmap(xlim = xlims,ylim = ylims) +
    borders(xlim = xlims,ylim = ylims, fill = "gray") +
    theme_pubr(I(18)) + 
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_x_longitude(xmin=-180, xmax=180, step=2) +
    scale_y_latitude(ymin=-180, ymax=180, step=2) +
    ggtitle(paste0("American_lobster_", season[j]))
  
  # # grid.arrange(p1, p2, nrow = 1)
  # 
  # png(paste0("/Users/Kisei/Desktop/Lobster_D_P0.05_", season[[j]], ".png"), 
  #     units = "in", res = 500, width = 8, height = 10) #res = 500, width=6000,height=6000
  # grid.arrange(p1, p2, nrow = 1)
  # dev.off()
  # 
  # png(paste0("/Users/Kisei/Desktop/Lobster_D_", season[[j]], ".png"), 
  #     units = "in", res = 500, width = 5, height = 16) #res = 500, width=6000,height=6000
  # p1
  # dev.off()
  # 
  # png(paste0("/Users/Kisei/Desktop/Lobster_P_", season[[j]], ".png"), 
  #     units = "in", res = 500, width = 5, height = 16) #res = 500, width=6000,height=6000
  # p2
  # dev.off()
  
  d1 = summarySE(ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], measurevar = "D", groupvars = c("Area"))
  d2 = summarySE(ks_df[which(ks_df$Area %in% c("GOM_GB", "SNE", "GOM_GB Neashore", "SNE Neashore")),], measurevar = "P_0.05", groupvars = c("Area"), na.rm = T)
  
  # write.csv(d1, paste0("/Users/Kisei/Desktop/", paste0("Lobster_", season[j]), "_D.csv"))
  # write.csv(d2, paste0("/Users/Kisei/Desktop/", paste0("Lobster_", season[j]), "_P.csv"))
 
  sp_map = function(x, y, z, df){
    
    coordinates(df) = ~x + y
    
    auto = autofitVariogram(z ~ 1, df)
    g = gstat(formula = z ~ 1, model = auto$var_model, data = df, maxdist = 0.05)
    
    xrange = range(df$x)
    yrange = range(df$y)
    grid= expand.grid(x = seq(from = xrange[1], to = xrange[2], by = 0.1), 
                      y = seq(from = yrange[1], to = yrange[2], by = 0.1))
    
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
  
  # jpeg(paste0("/Users/Kisei/Desktop/Lobster_Two-Sample_KS_test_", test, "_", season[[j]], ".jpg"), res = 500, height = 5, width = 7, units = "in")
  # spplot(D, 
  #        sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
  #        at = (0:1000)/1000,
  #        # col.regions = rev(pals::parula(1000)),
  #        col.regions = col(1000),
  #        zcol = c("Two_Sample_KS_Pval"),
  #        scales=list(draw=T),
  #        colorkey = T)
  # dev.off()
  # 
  # jpeg(paste0("/Users/Kisei/Desktop/Lobster_Two-Sample_KS_test_P0.05_", test, "_", season[[j]], ".jpg"), res = 500, height = 5, width = 7, units = "in")
  # spplot(D, 
  #        sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
  #        at = (0:50)/1000,
  #        # col.regions = rev(pals::parula(1000)),
  #        col.regions = col(1000),
  #        zcol = c("P_0.05"),
  #        scales=list(draw=T),
  #        colorkey = T)
  # dev.off()
  # 
  # jpeg(paste0("/Users/Kisei/Desktop/Lobster_Two-sample_KS_Tests_D_", test, "_", season[[j]], ".jpg"), res = 500, height = 5, width = 7, units = "in")
  # spplot(D, 
  #        sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
  #        at = (0:100)/100,
  #        # col.regions = pals::parula(length(na.omit(D$var1.pred))),
  #        col.regions = col(length(na.omit(D$var1.pred))),
  #        zcol = c("D"),
  #        scales=list(draw=T),
  #        colorkey = T)
  # dev.off()
  # 
  # jpeg(paste0("/Users/Kisei/Desktop/Lobster_Two-sample_KS_Tests_P0.05_D_", test, "_", season[[j]], ".jpg"), res = 500, height = 10, width = 7, units = "in")
  # spplot(D, 
  #        main=list(label="H. americanus",cex=1.5),
  #        sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
  #        at = (0:100)/100,
  #        # col.regions = rev(matlab.like(1000)),
  #        col.regions = col(1000),
  #        zcol = c("D", "P_0.05"),
  #        scales=list(draw=T),
  #        colorkey = T)
  # dev.off()
  
  #summarize TS test by mgmt areas
  
  df = cbind(ks, df*0.001)
  
  prob = NULL
  for (i in 1:4) {
    
    if(i == 1) {
      dd = df[which(df$GOMGBK_SNE == "GOM_GBK"),]
      Area = "GOM_GB"
    }
    
    if(i == 2) {
      dd = df[which(df$GOMGBK_SNE == "SNE"),]
      Area = "SNE"
    }
    
    if(i == 3) {
      dd = df[which(df$Nearshore_Area == "GOM_GBK Neashore Management Area"),]
      Area = "GOM_GB Neashore"}
    
    if(i == 4) {
      dd = df[which(df$Nearshore_Area == "SNE Neashore Management Area"),]
      Area = "SNE Neashore"
    }
    
    x = as.vector(t(dd[,10:29])) #first 10 years
    y = as.vector(t(dd[,70:length(dd)])) #last 10 years
    
    #make cdf comparison plots
    x = data.frame(x)
    y = data.frame(y)
    
    x$Period = "First_10_years"
    y$Period = "Last_10_years"
    
    colnames(x) = c("value", "period")
    colnames(y) = c("value", "period")
    
    xy = rbind(x,y)
    
    stat = tapply(xy$value, xy$period, summary)
    
    prob1 = stat[1]$First_10_years[3]
    prob2 = stat[2]$Last_10_years[3]
    
    prob_area = as.data.frame(cbind(prob1, prob2))
    prob_area$Area = Area
    prob_area$change = round(prob_area$prob2 - prob_area$prob1, 2)
    
    prob = rbind(prob, prob_area)
    
  }
  
  dd = df[which(df$GOMGBK_SNE == "GOM_GBK"),]
  Area = "GOM_GB"
  x = as.vector(t(dd[,10:29])) #first 10 years
  y = as.vector(t(dd[,70:length(dd)])) #last 10 years
  x = data.frame(x)
  y = data.frame(y)
  x$Period = "First_10_years"
  y$Period = "Last_10_years"
  colnames(x) = c("value", "period")
  colnames(y) = c("value", "period")
  xy_1 = rbind(x,y)
  xy_1$Area = Area
  
  dd = df[which(df$GOMGBK_SNE == "SNE"),]
  Area = "SNE"
  x = as.vector(t(dd[,10:29])) #first 10 years
  y = as.vector(t(dd[,70:length(dd)])) #last 10 years
  x = data.frame(x)
  y = data.frame(y)
  x$Period = "First_10_years"
  y$Period = "Last_10_years"
  colnames(x) = c("value", "period")
  colnames(y) = c("value", "period")
  xy_2 = rbind(x,y)
  xy_2$Area = Area
  
  dd = df[which(df$Nearshore_Area == "GOM_GBK Neashore Management Area"),]
  Area = "GOM_GB Neashore"
  x = as.vector(t(dd[,10:29])) #first 10 years
  y = as.vector(t(dd[,70:length(dd)])) #last 10 years
  x = data.frame(x)
  y = data.frame(y)
  x$Period = "First_10_years"
  y$Period = "Last_10_years"
  colnames(x) = c("value", "period")
  colnames(y) = c("value", "period")
  xy_3 = rbind(x,y)
  xy_3$Area = Area
  
  dd = df[which(df$Nearshore_Area == "SNE Neashore Management Area"),]
  Area = "SNE Neashore"
  x = as.vector(t(dd[,10:29])) #first 10 years
  y = as.vector(t(dd[,70:length(dd)])) #last 10 years
  x = data.frame(x)
  y = data.frame(y)
  x$Period = "First_10_years"
  y$Period = "Last_10_years"
  colnames(x) = c("value", "period")
  colnames(y) = c("value", "period")
  xy_4 = rbind(x,y)
  xy_4$Area = Area
  
  xy = rbind(xy_1, xy_2, xy_3, xy_4)
  
  library(dplyr)
  mu = xy %>%
    group_by(period, Area) %>% 
    summarise(median = median(value))
  
  # pdf(paste0("/Users/Kisei/Desktop/Lobster_Area_", season[[j]], ".pdf"), width = 3, height = 6) 
  png(paste0("/Users/Kisei/Desktop/Lobster_Area_", season[[j]], ".png"), width = 3, height = 6, res = 500, units = "in") 
  p = ggplot(xy, aes(x = value, fill = period, color = period)) +
    # geom_histogram(position="identity", alpha = 0.5, bins = 50) +
    geom_histogram(aes(y=..count..), position="identity", alpha = 0.5, bins = 50) +
    # geom_histogram(aes(y=..count../sum(..count..)), position="identity", alpha = 0.5, bins = 50) +
    # geom_histogram(aes(y=..ncount..), position="identity", alpha = 0.5, bins = 50) +
    # geom_density(alpha = 0.3) +
    facet_wrap(~Area, scales = "free_y", ncol = 1) +
    geom_vline(data=mu, aes(xintercept=median, color=period),
               linetype="dashed",
               size = 0.5) +
    xlab("Habitat Suitability") +
    ylab(" ") +
    xlim(0,1)+ 
    theme_pubr(I(10)) + ggtitle(paste0("Lobster_", toupper(season[[j]])))  + 
    theme(legend.position="none") + 
    # geom_label(data = prob, aes(label= paste0("change=", change)), 
    #            x = Inf, y = Inf, hjust = 1, vjust = 1, label.size = NA,
    #            inherit.aes = FALSE)
    geom_label(data = prob, aes(label= paste0("First_10_yrs=", round(prob1, 2), "\nLast_10_yrs=", round(prob2, 2))),
               x = Inf, y = Inf, hjust = 1, vjust = 1, label.size = NA, size = 3, alpha = 0.5,
               inherit.aes = FALSE)
  print(p)
  dev.off()
  
  save(xy, file = paste0("/Users/Kisei/Desktop/Lobster_Area_", season[[j]], ".RData"))

}

# The significance of difference between two REC curves can be assessed by examining the maximum deviation between the two curves across all values of E. This corresponds to the Kolmogorov-Smirnov (KS) two-sample test (DeGroot, 1986) for judging the hypothesis that the error e generated by two models f and g follows the same distribution. The KS-test requires no assumptions about the underlying CDF of z. 