library(sp)
library(maptools)
library(colorRamps)
library(ggplot2)
library(readr)
library(Rmisc)
library(ggpmisc)
library(ggpubr)

# setwd("/Users/Kisei/Google Drive/R/Biomod/scallop")
setwd("/Users/Kisei/Desktop/scallop")

season = list("fall", "spring", "annual")[3]

for (j in 1:length(season)){
  
  j = 1
  
  df <- read_csv(paste0("Biomod_1_80_", season[j],".csv")) #biomod output

  ks = df[,c(1,2)]
  df = df[3:82]
  
  #add GOMGBK and SNE area
  colnames(ks)[1:2] = c("x","y")
  latlon = ks[,c(1,2)]
  coordinates(latlon)=~x+y
  area = rgdal::readOGR("/Users/Kisei/Google Drive/Research/GIS/NOAA_Statistical_Area/Statistical_Areas.shp")
  area = area[which(area$Id %in% c(464:465, 511:515, 521:526, 537:539, 541:543, 561:562, 611:616, 621:626)),]
  CRS.new = CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
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
  
  df = cbind(ks, df)
  
  df$ID = paste(df$x, df$y, sep = "_")
  
  # data = NULL
  # 
  # for (i in 1:80) {
  #   
  #   # i = 80
  #   dd = df[,c(3:4, 85, i+4)]
  #   dd$Year = i
  #   colnames(dd) = c("Management_Areas", "Nearshore_Areas", "ID", "HS", "Year")
  #   data <- rbind(data, dd)
  #   
  # }
  
  data = vector("list")
  
  for (i in 1:80) {
    
    # i = 80
    dd = df[,c(3:4, 85, i+4)]
    dd$Year = i
    colnames(dd) = c("Management_Areas", "Nearshore_Areas", "ID", "HS", "Year")
    data[[i]] = dd
    print(i)
    
  }
  
  # data = rbind(data, do.call(rbind, data.list))
  data = rbindlist(data)
  
  data$Year = as.numeric(data$Year)
  data$HS = data$HS/1000
  # data = subset(data, Management_Areas %in% c("GOM_GBK", "SNE"))
  # data = subset(data, Nearshore_Areas %in% c("GOM_GBK Neashore Management Area", "SNE Neashore Management Area"))
  d1 = summarySE(data[which(data$Management_Areas == "GOM_GBK"),], measurevar = "HS", groupvars = c("Year"))
  d2 = summarySE(data[which(data$Management_Areas == "SNE"),], measurevar = "HS", groupvars = c("Year"))
  d3 = summarySE(data[which(data$Nearshore_Areas == "GOM_GBK Neashore Management Area"),], measurevar = "HS", groupvars = c("Year"))
  d4 = summarySE(data[which(data$Nearshore_Areas == "SNE Neashore Management Area"),], measurevar = "HS", groupvars = c("Year"))
  
  d1$Management_Area = "GOM_GB"
  d2$Management_Area = "SNE"
  d3$Management_Area = "GOM_GB Neashore"
  d4$Management_Area = "SNE Neashore"
  
  data = rbind(d1, d2, d3, d4)
  data = data[c("Year", "HS", "Management_Area")]
  colnames(data) = c("x", "y", "Management_Area")

  r_df <- ddply(data, .(Management_Area), summarise, r=round(summary(lm(scale(y)~x))$adj.r.squared, 2))
  b_df <- ddply(data, .(Management_Area), summarise, b=round(summary(lm(scale(y)~x))$coefficients[2], 2))
  p_df <- ddply(data, .(Management_Area), summarise, p=round(summary(lm(scale(y)~x))$coefficients[8], 2))
  p_df$p = ifelse(p_df$p < 0.05, "<0.05", paste("=", p_df$p))
  r_df = merge(r_df, p_df)
  r_df = merge(r_df, b_df)
  r_df = merge(r_df, b_df)
  
  p = ggplot(data, aes(x = x, y = y, color = Management_Area)) +
    geom_text(data = r_df, 
              # aes(label = paste0("r^2: ", r)), 
              # aes(label = paste0(" r?=", r, ", p", p)), 
              # aes(label = paste0("slope=", b, "\nr?=", r, "\np", p)), 
              aes(label = paste0("??=", b, ", r?=", r, ", p", p)), 
              x = Inf, y = -Inf, 
              hjust = 1,
              vjust = -0.4,
              size = 4) +
    geom_point(alpha = 0.2) + 
    # geom_line(alpha = 0.2) + 
    geom_smooth(method = "lm")  +
    theme_pubr(base_size = I(10)) +
    theme(legend.position="none") + 
    facet_wrap( ~ Management_Area, scales="free_y", ncol = 1) +
    labs(x = "Model Year", y = "Relative Habitat Suitability") + ggtitle(paste0("Sea scallop_", season[j]))
  
  p
  
  png(paste0("/Users/Kisei/Desktop/scallop_change_", season[j], ".png"), res = 500, height = 10, width = 2.5, units = "in")
  print(p)
  dev.off()
  
  write_csv(r_df, paste0("/Users/kisei/Desktop/scallop_", season[j], "_ensemble_area_trend.csv"))
  
  
}
