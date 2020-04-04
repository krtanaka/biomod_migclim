library(ggpmisc)
library(ggplot2)
library(sp)
library(Rmisc)
library(ggpubr)
library(gridExtra)

# attach necessary GIS information to cm2.6 data, skip after first --------

bt = Bottom_Temperature

ks = bt[,c(1,2)]
df = bt[3:962]

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

area = c('GOM_GBK', 'SNE', 'GOM_GBK Neashore Management Area', 'SNE Neashore Management Area')

get_area_average = function(area, df){
  
  if (area %in% c('GOM_GBK', 'SNE')) df = df[which(df$GOMGBK_SNE == area),]
  if (area %in% c('GOM_GBK Neashore Management Area', 'SNE Neashore Management Area')) df = df[which(df$Nearshore_Area== area),]
  
  data1 = NULL
  
  for (i in 1:160) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data1 <- rbind(data1, dd)
    
  }
  
  data2 = NULL
  
  for (i in 161:320) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data2 <- rbind(data2, dd)
    
  }
  
  data3 = NULL
  
  for (i in 321:480) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data3 <- rbind(data3, dd)
    
  }
  
  data4 = NULL
  
  for (i in 481:640) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data4 <- rbind(data4, dd)
    
  }
  
  data5 = NULL
  
  for (i in 641:800) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data3 <- rbind(data3, dd)
    
  }
  
  data6 = NULL
  
  for (i in 801:960) {
    
    print(i)
    dd = df[,c(1:4, i+4)]
    dd$Year = i
    colnames(dd) = c("x", "y", "GOMGBK_SNE", "Nearshore_Area", "bt", "Year")
    data4 <- rbind(data4, dd)
    
  }
  
  data = rbind(data1, data2, data3, data4, data5, data6)
  
  data$Year = as.numeric(data$Year)
  
  data = data[,c("Year", "bt"),]
  
  d = summarySE(data, measurevar = "bt", groupvars = c("Year"))
  
  d$Management_Area = area
  
  return(d)
  
}

d1 = get_area_average("GOM_GBK", df)
d2 = get_area_average("SNE", df)
d3 = get_area_average("GOM_GBK Neashore Management Area", df)
d4 = get_area_average("SNE Neashore Management Area", df)

d = rbind(d1, d2, d3, d4)
d$Management_Area <- gsub('GBK', 'GB', d$Management_Area)
d$Management_Area <- gsub(' Neashore Management Area', ' Neashore', d$Management_Area)

colnames(d)[1] = "x"; colnames(d)[3] = "y"

# compute trend in each area within NEUS-LME ------------------------------
load("/Users/Kisei/Google Drive/Research/NW_Clim/Data/cm2.6_w_NWA_GIS_bt.RData")
load("/Users/Kisei/Google Drive/Research/NW_Clim/Data/cm2.6_w_NWA_GIS_bs.RData")

bt$Property = "Bottom Temperature Anomalies (deg C)"
bs$Property = "Bottom Salinity Anomalies (ppt)"

d = rbind(bt, bs)

d$Property = factor(d$Property, levels=c("Bottom Temperature Anomalies (deg C)","Bottom Salinity Anomalies (ppt)" ))

b_df <- ddply(d, .(Management_Area, Property), summarise, slope=round(summary(lm(scale(y)~x))$coefficients[2], 3))
p_df <- ddply(d, .(Management_Area, Property), summarise, p=summary(lm(scale(y)~x))$coefficients[8])
p_df$p = ifelse(p_df$p < 0.001, "<0.001", paste0("", p_df$p))
summary = merge(b_df, p_df)
summary

p = ggplot(d, aes(x, y, color = Management_Area)) +
  geom_text(data = summary,
            aes(label = paste0("\n β=", slope, "\n p", p)),
            x = -Inf, y = -Inf,
            hjust = -0.1,
            vjust = -0.2,
            size = 3) +
  # geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 1) +
  facet_grid(Property~Management_Area, scales = "free") +
  scale_x_continuous(breaks = seq(0, 960, by = 240)) + 
  theme_pubr(base_size = I(10)) + 
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.7)) + 
  labs(x = "Model Month", y = "")

png("/Users/Kisei/Desktop/bt_bs.png", res = 500, height = 5, width = 6, units = "in")
p
dev.off()

formula = y ~ x

ggplot(d, aes(x, y)) +
  geom_point(aes(color = Management_Area), size = 1, alpha = 0.1) +
  facet_grid(Property~Management_Area, scales = "free") + 
  theme_pubr(base_size = I(10)) +
  theme(legend.position="none") +
  stat_smooth( aes(color = Management_Area, fill = Management_Area), method = "lm") +
  # stat_cor(aes(color = Management_Area), label.y = 4.6) +
  stat_cor(aes(color = Management_Area, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = c(2, 5)) + 
  stat_poly_eq(
    aes(color = Management_Area, label = ..eq.label..),
    formula = formula, label.y = c(1, 4), coef.digits = 2, parse = TRUE) + 
  labs(x = "Model Month")


d = bt

formula = y ~ x


p1 <- ggplot(d, aes(x, y)) +
  geom_point(aes(color = Management_Area), size = 1, alpha = 0.2) +
  facet_wrap(~Management_Area, ncol = 1)+ 
  theme_pubr(base_size = I(10)) +
  theme(legend.position="none") +
  stat_smooth( aes(color = Management_Area, fill = Management_Area), method = "lm") +
  # stat_cor(aes(color = Management_Area), label.y = 4.6) +
  stat_cor(aes(color = Management_Area, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = 5) + 
  stat_poly_eq(
    aes(color = Management_Area, label = ..eq.label..),
    formula = formula, label.y = 4, coef.digits = 2, parse = TRUE) + 
  labs(x = "Model Month", y = "Bottom temperature anomaly (deg C)") + ggtitle("")

d = bs

formula = y ~ x


p2 <- ggplot(d, aes(x, y)) +
  geom_point(aes(color = Management_Area), size = 1, alpha = 0.2) +
  facet_wrap(~Management_Area, ncol = 1)+ 
  theme_pubr(base_size = I(10)) +
  theme(legend.position="none") +
  stat_smooth( aes(color = Management_Area, fill = Management_Area), method = "lm") +
  # stat_cor(aes(color = Management_Area), label.y = 4.6) +
  stat_cor(aes(color = Management_Area, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = 1.1) + 
  stat_poly_eq(
    aes(color = Management_Area, label = ..eq.label..),
    formula = formula, label.y = 0.9, coef.digits = 2, parse = TRUE) + 
  labs(x = "Model Month", y = "Bottom salinity anomaly (ppt)") + ggtitle("")

grid.arrange(p1, p2, nrow = 1)

png("/Users/Kisei/Google Drive/Research/Manuscripts/Biomod/bt.png", res = 500, height = 10, width = 4, units = "in")
p1
dev.off()

png("/Users/Kisei/Google Drive/Research/Manuscripts/Biomod/b2.png", res = 500, height = 10, width = 4, units = "in")
p2
dev.off()
