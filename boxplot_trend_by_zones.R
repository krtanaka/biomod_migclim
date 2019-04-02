library(sp)
library(maptools)
library(colorRamps)
library(ggplot2)
library(readr)
library(Rmisc)

setwd("/Users/Kisei/Google Drive/R/Biomod/scallop")

df <- read_csv("Biomod_1_80_annual.csv") #lobster fall biomod output
names(df) 

ks = df[,c(85,84)]
df = df[1:80]

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
                            "GOM_GBK Neashore Management Area", 
                            ifelse(ks$Nearshore_Area %in% c("EEZ Nearshore Management Area 2", 
                                                            "EEZ Nearshore Management Area 4"), 
                                   "SNE Neashore Management Area",
                                   "NA"))

df = cbind(ks, df)

df$ID = paste(df$x, df$y, sep = "_")

data = NULL

for (i in 1:80) {
  
  # i = 80
  dd = df[,c(3:4, 85, i+4)]
  dd$Year = i
  colnames(dd) = c("Management_Areas", "Nearshore_Areas", "ID", "HS", "Year")
  data <- rbind(data, dd)
  
}

data$Year = as.numeric(data$Year)
data$HS = data$HS/1000
# data = subset(data, Management_Areas %in% c("GOM_GBK", "SNE"))
# data = subset(data, Nearshore_Areas %in% c("GOM_GBK Neashore Management Area", "SNE Neashore Management Area"))
d1 = summarySE(data[which(data$Management_Areas == "GOM_GBK"),], measurevar = "HS", groupvars = c("Year"))
d2 = summarySE(data[which(data$Management_Areas == "SNE"),], measurevar = "HS", groupvars = c("Year"))
d3 = summarySE(data[which(data$Nearshore_Areas == "GOM_GBK Neashore Management Area"),], measurevar = "HS", groupvars = c("Year"))
d4 = summarySE(data[which(data$Nearshore_Areas == "SNE Neashore Management Area"),], measurevar = "HS", groupvars = c("Year"))

d1$Management_Area = "GOM_GBK"
d2$Management_Area = "SNE"
d3$Management_Area = "GOM_GBK Neashore"
d4$Management_Area = "SNE Neashore"

data = rbind(d1, d2, d3, d4)

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


  
p = ggplot(data, aes(x = Year, y = HS, color = Management_Area)) + 
  geom_point( aes(group = Year)) +
  # geom_line() +
  geom_smooth(method = "lm")  + 
  theme_classic(base_size = I(10))+ 
  theme(legend.position="none") + 
  facet_wrap( ~ Management_Area, scales="free", ncol = 1) +
  labs(x = "Year", y = "") + ggtitle("Atlantic scallop")

png("scallop_change_1_80.png", res = 500, height = 8, width = 2, units = "in")
print(p)
dev.off()
