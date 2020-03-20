library(sp)

rm(list = ls())

scallop = read.csv("/Users/kisei/Google Drive/R/Biomod/scallop/Biomod_1_80_annual.csv")
Lob_fall = read.csv("/Users/kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_fall.csv")
Lob_spring = read.csv("/Users/kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_spring.csv")

ks = scallop[,c(1,2)]
df_ensemble = scallop[3:82]

#add GOMGBK and SNE area
colnames(ks)[1:2] = c("x","y")
latlon = ks[,c(1,2)]; plot(latlon)
coordinates(latlon)=~x+y
area = rgdal::readOGR("/Users/kisei/Google Drive/Research/GIS/NOAA_Statistical_Area/Statistical_Areas.shp")
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
area = rgdal::readOGR("/Users/kisei/Google Drive/Research/GIS/Lobster_Management_Areas/Lobster_Management_Areas.shp")
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

Scallop_annual = cbind(scallop, ks[,c(3:4)])
Lobster_spring = cbind(Lob_spring, ks[,c(3:4)])
Lobster_fall = cbind(Lob_fall, ks[,c(3:4)])

# load('habitat_change.RData')


L = list(Scallop_annual = Scallop_annual, Lobster_fall = Lobster_fall, Lobster_spring = Lobster_spring)

## High = 666-1000
## Med = 333 - 666
## low= < 333

df = NULL

Area = c("GOM_GB", "SNE", "SNE Nearshore", "GOM_GB Nearshore")

for (a in 1:length(Area)) {
  
  area = Area[[a]]
  
  for (k in 1:length(L)){
    
    x = L[[k]] 
    
    ct = data.frame(Start_Quant = c('High','High','High', 'Med','Med','Med', 'Low','Low','Low'),
                    End_Quant = c('High', 'Med', 'Low','High', 'Med', 'Low','High', 'Med', 'Low'),
                    Frequency = rep(NA,9 ))
    
    if (area == "GOM_GB") x = subset(x, GOMGBK_SNE == "GOM_GBK")
    if (area == "SNE") x = subset(x, GOMGBK_SNE == "SNE")
    if (area == "GOM_GB Nearshore") x = subset(x, Nearshore_Area == "GOM_GBK Nearshore Management Area")
    if (area == "SNE Nearshore") x = subset(x, Nearshore_Area == "SNE Nearshore Management Area")
    
    x$mean1 = apply(x[,3:12], 1 , mean )
    x$mean2 = apply(x[,72:82], 1 , mean )
    
    highx = which(x$mean1 > 666)
    medx = which(x$mean1 <= 666 & x$mean1 >= 333 )
    lowx = which(x$mean1 < 333 )
    
    ct[1,3] = length(which(x$mean2[highx] > 666))
    ct[2,3] = length(which(x$mean2[highx] <= 666 & x$mean2[highx] >= 333))
    ct[3,3] = length(which(x$mean2[highx] < 333))
    
    ct[4,3] = length(which(x$mean2[medx] > 666))
    ct[5,3] = length(which(x$mean2[medx] <= 666 & x$mean2[medx] >= 333))
    ct[6,3] = length(which(x$mean2[medx] < 333))
    
    ct[7,3] = length(which(x$mean2[lowx] > 666))
    ct[8,3] = length(which(x$mean2[lowx] <= 666 & x$mean2[lowx] >= 333))
    ct[9,3] = length(which(x$mean2[lowx] < 333))
    
    ct$Sp = names(L)[k]
    ct$Area = area
    
    ct$Start_Quant <- factor(ct$Start_Quant, levels = c("High", "Med", "Low"))
    ct$End_Quant <- factor(ct$End_Quant, levels = c("High", "Med", "Low"))
    
    df = rbind(df, ct) 
    
  }
  
}

df

df[is.na(df)] <- 0

df = df %>% 
  group_by(Sp) %>%
  mutate(a = Frequency/sum(Frequency))

df = as.data.frame(df)

df1 = df[,c("Start_Quant", "Frequency", "Sp", "a", "Area")]
df2 = df[,c("End_Quant", "Frequency", "Sp", "a", "Area")]

df1$Subject = 1:dim(df1)[1]
df2$Subject = 1:dim(df2)[1]

colnames(df1) = c("response", "freq", "sp", "a", "area", "subject")
colnames(df2) = c("response", "freq", "sp", "a", "area", "subject")

df1$survey = "p1"
df2$survey = "p2"

df1 = df1 %>% 
  group_by(sp, area) %>%
  mutate(a = freq/sum(freq))

df2 = df2 %>% 
  group_by(sp, area) %>%
  mutate(a = freq/sum(freq))

df = rbind(df1, df2)

# df = df %>% 
#   group_by(sp) %>%
#   mutate(a = freq/sum(freq))

png('/Users/kisei/Desktop/Habitat_quantiles.png', width = 10, height = 10, units = "in", res = 500)

ggplot(df,
       aes(x = survey, stratum = response, alluvium = subject,
           y = freq, fill = response, label = round(a, 2))) +
  geom_lode() + geom_flow() +
  geom_stratum(alpha = 0) +
  geom_text(stat = "stratum", size = 4) + 
  theme_pubr(I(15)) + 
  # facet_wrap(~sp + area, scales = "free_y", dir = "v", ncol = 3) +
  facet_grid(area ~ sp, scales = "free") +
  
  scale_y_continuous("Freqency")+
  scale_x_discrete("")+
  scale_fill_manual(values = rev(matlab.like(4)), "Habitat Suitability Category") + 
  theme(legend.position = c("bottom"))

dev.off()



for (kk in 1:length(unique(df$Area))){
  
  dff = subset(df, Area==unique(df$Area)[kk])
  
  
  pp = ggplot(data = dff, aes(y = Frequency, axis1 = Start_Quant, axis2 = End_Quant)) + 
    geom_alluvium(aes(fill = Start_Quant), width = 1/12) + 
    #facet_wrap(~Area+Sp ) + 
    facet_wrap(~Sp ) + 
    geom_stratum(width = 1/12,color = "black", lwd=1.5) +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits = c("1-10 years", "70-80 years"), expand = c(.05, .05))+
    scale_fill_viridis_d("Habitat_Category")
  assign(paste("Area",kk, sep="_"),pp)
  
}

png('GOM_GB.png', width =12, height = 10,units = "in", res = 300)
plot(Area_1)
dev.off()

png('SNE.png', width =12, height = 10,units = "in", res = 300)
plot(Area_2)
dev.off()

png('SNE Nearshore.png', width =12, height = 10,units = "in", res = 300)
plot(Area_3)
dev.off()

png('GOM_GB Nearshore.png', width =12, height = 10,units = "in", res = 300)
plot(Area_4)
dev.off()










