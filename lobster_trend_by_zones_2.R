library(sp)
library(maptools)
library(colorRamps)
library(ggplot2)
library(readr)
library(Rmisc)
library(ggpmisc)
library(ggpubr)
library(data.table)

setwd("/Users/Kisei/Google Drive/R/Biomod/lobster/")

season = list("fall", "spring")

load("final_sdms.RData")

for (j in 1:length(season)){
  
  # j = 2
  
  df <- read_csv(paste0("Biomod_1_80_", season[j],".csv"))
  ks = df[,c(1,2)]
  df_ensemble = df[3:82]
  
  df_individual = list()
  
  for (m in 1:length(final_sdms)) {
    
    df = read_csv(paste0("Biomod_1_80_", final_sdms[m], "_", season[j], ".csv"))
    df = df[4:83]
    df_individual[[m]] = df
    
  }
  
  #add GOMGBK and SNE area
  colnames(ks)[1:2] = c("x","y")
  latlon = ks[,c(1,2)]; plot(latlon)
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
                              "GOM_GBK Nearshore Management Area", 
                              ifelse(ks$Nearshore_Area %in% c("EEZ Nearshore Management Area 2", 
                                                              "EEZ Nearshore Management Area 4"), 
                                     "SNE Nearshore Management Area",
                                     "NA"))
  
  df_ensemble = cbind(ks, df_ensemble); df_ensemble$Type = "Ensemble"
  
  df_individual_model = NULL
  
  for (m in 1:length(final_sdms)) {
    
    model = df_individual[[m]]
    model = cbind(ks, model)
    model$Type = final_sdms[[m]]
    df_individual_model = rbind(df_individual_model, model)
    
  }
  
  df = rbind(df_ensemble, df_individual_model)
  
  df$ID = paste(df$x, df$y, sep = "_")
  
  data.list = vector("list")
  
  for (i in 1:80) {
    
    # i = 80
    dd = df[,c(3:4, 85, 86, i+4)]
    dd$Year = i
    colnames(dd) = c("Management_Areas", "Nearshore_Areas", "Type", "ID", "HS", "Year")
    data.list[[i]] = dd
    print(i)
    
  }
  
  # data = rbind(data, do.call(rbind, data.list))
  data = rbindlist(data.list)
  
  data$Year = as.numeric(data$Year)
  data$HS = data$HS/1000
  
  data_ensemble = subset(data, Type == "Ensemble")
  d1_e = summarySE(data_ensemble[which(data_ensemble$Management_Areas == "GOM_GBK"),], measurevar = "HS", groupvars = c("Year"))
  d2_e = summarySE(data_ensemble[which(data_ensemble$Management_Areas == "SNE"),], measurevar = "HS", groupvars = c("Year"))
  d3_e = summarySE(data_ensemble[which(data_ensemble$Nearshore_Areas == "GOM_GBK Nearshore Management Area"),], measurevar = "HS", groupvars = c("Year"))
  d4_e = summarySE(data_ensemble[which(data_ensemble$Nearshore_Areas == "SNE Nearshore Management Area"),], measurevar = "HS", groupvars = c("Year"))
  d1_e$Management_Area = "GOM_GB"
  d2_e$Management_Area = "SNE"
  d3_e$Management_Area = "GOM_GB Nearshore"
  d4_e$Management_Area = "SNE Nearshore"
  d1_e$Type = "Ensemble"
  d2_e$Type = "Ensemble"
  d3_e$Type = "Ensemble"
  d4_e$Type = "Ensemble"
  
  ensemble_model = rbind(d1_e, d2_e, d3_e, d4_e)
  
  individual_model = NULL
  
  for (m in 1:length(final_sdms)) {
    
    data_model = subset(data, Type == final_sdms[m])
    d1_b = summarySE(data_model[which(data_model$Management_Areas == "GOM_GBK"),], measurevar = "HS", groupvars = c("Year"))
    d2_b = summarySE(data_model[which(data_model$Management_Areas == "SNE"),], measurevar = "HS", groupvars = c("Year"))
    d3_b = summarySE(data_model[which(data_model$Nearshore_Areas == "GOM_GBK Nearshore Management Area"),], measurevar = "HS", groupvars = c("Year"))
    d4_b = summarySE(data_model[which(data_model$Nearshore_Areas == "SNE Nearshore Management Area"),], measurevar = "HS", groupvars = c("Year"))
    d1_b$Management_Area = "GOM_GB"
    d2_b$Management_Area = "SNE"
    d3_b$Management_Area = "GOM_GB Nearshore"
    d4_b$Management_Area = "SNE Nearshore"
    d1_b$Type = final_sdms[m]
    d2_b$Type = final_sdms[m]
    d3_b$Type = final_sdms[m]
    d4_b$Type = final_sdms[m]
    
    individual_model = rbind(individual_model, d1_b, d2_b, d3_b, d4_b)
  }
  
  rm(d1_b, d2_b, d3_b, d4_b)
  
  data = rbind(ensemble_model, individual_model)
  data = data[c("Year", "HS", "Management_Area", "Type")]
  colnames(data) = c("x", "y", "Management_Area", "Model")
  
  r <- ddply(data[which(data$Model == "Ensemble"),], .(Management_Area), summarise, r = round(summary(lm(scale(y)~x))$adj.r.squared, 2))
  b <- ddply(data[which(data$Model == "Ensemble"),], .(Management_Area), summarise, b = round(summary(lm(scale(y)~x))$coefficients[2], 2))
  p <- ddply(data[which(data$Model == "Ensemble"),], .(Management_Area), summarise, p = round(summary(lm(scale(y)~x))$coefficients[8], 2))
  p$p = ifelse(p$p < 0.05, "<0.05", paste("=", p$p))
  ensemble_coef = merge(r, p)
  ensemble_coef = merge(ensemble_coef, b)
  ensemble_coef$Model = "Ensemble"
  
  individual_coef = NULL
  
  for (m in 1:length(final_sdms)) {
    
    r <- ddply(data[which(data$Model == final_sdms[m]),], .(Management_Area), summarise, r = round(summary(lm(scale(y)~x))$adj.r.squared, 2))
    b <- ddply(data[which(data$Model == final_sdms[m]),], .(Management_Area), summarise, b = round(summary(lm(scale(y)~x))$coefficients[2], 2))
    p <- ddply(data[which(data$Model == final_sdms[m]),], .(Management_Area), summarise, p = round(summary(lm(scale(y)~x))$coefficients[8], 2))
    p$p = ifelse(p$p < 0.05, "<0.05", paste("=", p$p))
    ind_coef = merge(r, p)
    ind_coef = merge(ind_coef, b)
    ind_coef$Model = final_sdms[m]
    individual_coef = rbind(individual_coef, ind_coef)
  }
  
  
  data$Model = gsub("lobster_AllData_RUN", "", data$Model)
  data$Model = paste0(substr(data$Model, 3,5), substr(data$Model, 1,1))
  data$Model = gsub("semE", "Ensemble", data$Model)
  
  coef = rbind(ensemble_coef, individual_coef)
  coef$Model = gsub("lobster_AllData_RUN", "", coef$Model)
  coef$Model = paste0(substr(coef$Model, 3,5), "_", substr(coef$Model, 1,1))
  coef$Model = gsub("sem_E", "Ensemble", coef$Model)
  
  data$Lobster = data$Model
  
  jpeg(paste0("/Users/Kisei/Desktop/lobster_change_legend.jpg"), res = 500, height = 1, width = 4, units = "in")
  legend = ggplot(data, aes(x = x, y = y, color = Lobster)) +
    geom_line( data = subset( data, !(Model %in% 'Ensemble') ),
               # color = I( 'black' ),
               # show.legend = FALSE,
               size = I( 2 )) +  
    theme(legend.position = "bottom") 
  legend = cowplot::get_legend(legend)
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
  
  p = ggplot(data, aes(x = x, y = y, color = Model)) +
    geom_line( data = subset( data, !(Model %in% 'Ensemble') ),
               # method = "lm", se = F, 
               # show.legend = FALSE,
               # alpha = 0.2, 
               size = I( 0.2 )) + 
    geom_line( data = subset( data, Model %in% 'Ensemble' ),
               color = I( 'black' ),
               show.legend = FALSE,
               size = I( 1 )) +
    theme_pubr(base_size = I(10)) +
    theme(legend.position="none") + 
    facet_wrap( ~ Management_Area, scales="free_y", ncol = 4) +
    labs(x = "Model Year", y = "Relative Habitat Suitability") + 
    ggtitle(paste0("American_lobster_", season[j]))
  
  p
  
  png(paste0("/Users/Kisei/Desktop/lobster_change_", season[j], ".png"), 
      res = 500, height = 3, width = 8, units = "in")
  print(p)
  dev.off()
  
  prob = coef[,c("p", "b", "Model", "Management_Area")]
  prob = subset(prob, Model != "Ensemble")
  
  prob$outcome <- with(prob, ifelse(p == "<0.05" & b > 0, "+", "-"))
  prob$outcome <- with(prob, ifelse(p != "<0.05", "NA", outcome))
  prob$increase = ifelse(prob$outcome == "+", 1, 0)
  prob$decrease = ifelse(prob$outcome == "-", 1, 0)
  prob$no_change = ifelse(prob$outcome == "NA", 1, 0)
  
  x = prob[,c("Management_Area", "increase", "decrease", "no_change")]
  
  x = aggregate(. ~ Management_Area, data = x, FUN = sum)
  x$total = rowSums(x[2:4])
  x$increase = x$increase/x$total
  x$decrease = x$decrease/x$total
  x$no_change = x$no_change/x$total
  x = x[,c("Management_Area", "increase", "decrease", "no_change")]
  x
  rowSums(x[2:4])
  
  y = prob[,c("Management_Area", "Model", "increase", "decrease", "no_change")]
  colnames(y)[2] = "SDM"
  y = merge(y, mw)
  y$increase = y$increase*y$weight
  y$decrease = y$decrease*y$weight
  y$no_change = y$no_change*y$weight
  y = y[,c("Management_Area", "increase", "decrease", "no_change")]
  y = aggregate(. ~ Management_Area, data = y, FUN = sum)
  y$total = rowSums(y[2:4])
  y$increase = y$increase/y$total
  y$decrease = y$decrease/y$total
  y$no_change = y$no_change/y$total
  y = y[,c("Management_Area", "increase", "decrease", "no_change")]
  y
  rowSums(y[2:4])
  
  colnames(x)[2:4] = c("Increase_Unweighted", "Decrease_Unweighted", "No_Change_Unweighted")
  colnames(y)[2:4] = c("Increase_Weighted", "Decrease_Weighted", "No_Change_Weighted")
  
  sdm = merge(x, y, all = T)
  sdm = sdm[,c(1,2,5,3,6,4,7)]
  sdm
  
  write_csv(sdm, paste0("/Users/kisei/Desktop/lobster_", season[j], "_sdm_area_agreement.csv"))
  
  
}
