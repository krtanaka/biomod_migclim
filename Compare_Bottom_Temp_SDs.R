d = data.table::fread(file = "/Users/Kisei/Desktop/Bottom_Temperature.csv")

names(d)
latlon = d[,c(1:2)]; plot(latlon)
coordinates(latlon)=~Var1+Var2
lme<-rgdal::readOGR("/Users/Kisei/Google Drive/Research/GIS/NOAA_Statistical_Area/Statistical_Areas.shp")
CRS.new<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
proj4string(latlon) <- CRS.new 
proj4string(lme) <- CRS.new
area <- over(latlon,lme)
colnames(area)[1] = "stat"
d = cbind(d, area)
d = subset(d, stat %in% c("464", "465", "511", "512", "513", "514", "515", "521", "522", "525", "526", "537", "538", "539", "551", "552", "561", "562", "611", "612", "613", "614", "615", "616", "621", "622", "623", "625", "626"))

cm = colMeans(d[,3:962])
  
cm_month_sd = NULL

for (i in 1:12) {

  interval = seq(i, 960, 12)
sd = sd(cm[interval])

cm_month_sd = rbind(cm_month_sd, sd)

}

cm_month_sd = as.data.frame(cbind(c(1:12), cm_month_sd))

load("/Users/Kisei/biomod_migclim/lobster/lobster_survey_data_spring_fall_combined_1984-2016.RData") #load survey data
d = lobster

d = d[,c("Month", "Obs_Bottom_Temperature")]

d = d %>% group_by(Month) %>% summarise_each(funs(sd))

colnames(d) = c("Month", "Observed Bottom Temp SD")
colnames(cm_month_sd) = c("Month", "Modeled Bottom Temp SD")

d = merge(d, cm_month_sd)

setwd("C:/Users/Kisei/Desktop")
readr::write_csv(d, "sd_comparison.csv")


