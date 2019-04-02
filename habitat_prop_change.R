library(sp)
library(ggalluvial)

rm(list = ls())

scallop = read.csv("/Users/Kisei/Google Drive/R/Biomod/scallop/Biomod_1_80_annual.csv")
Lob_fall = read.csv("/Users/Kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_fall.csv")
Lob_spring = read.csv("/Users/Kisei/Google Drive/R/Biomod/lobster/Biomod_1_80_spring.csv")

ks = scallop[,c(1,2)]
df_ensemble = scallop[3:82]

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

Scallop_annual = cbind(scallop, ks[,c(3:4)])
Lobster_spring = cbind(Lob_spring, ks[,c(3:4)])
Lobster_fall = cbind(Lob_fall, ks[,c(3:4)])

L = list(Scallop_annual = Scallop_annual, Lobster_fall = Lobster_fall, Lobster_spring = Lobster_spring)

## High = 666-1000
## Med = 333 - 666
## low= < 333

df = NULL
cc = NULL

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
    
    x1 = as.data.frame(x[,c("mean1")])
    x2 = as.data.frame(x[,c("mean2")])
    
    x1$P = "P1"
    x2$P = "P2"
    
    colnames(x1) = c("value", "period")
    colnames(x2) = c("value", "period")
    
    h1 = which(x1$value > 666)
    m1 = which(x1$value <= 666 & x1$value >= 333 )
    l1 = which(x1$value < 333 )
    
    h2 = which(x2$value > 666)
    m2 = which(x2$value <= 666 & x2$value >= 333 )
    l2 = which(x2$value < 333 )
    
    d1 = as.data.frame(rbind(c("High", length(h1)), 
                             c("Med", length(m1)), 
                             c("Low", length(l1))))
    
    d2 = as.data.frame(rbind(c("High", length(h2)), 
                             c("Med", length(m2)), 
                             c("Low", length(l2))))
    
    d1$p = "1-10 yrs"
    d2$p = "70-80 yrs"
    
    colnames(d1) = c("Habitat_Category", "Freq", "Period")
    colnames(d2) = c("Habitat_Category", "Freq", "Period")
    
    d1$Prop = as.numeric(d1$Freq)/sum(as.numeric(d1$Freq))
    d2$Prop = as.numeric(d2$Freq)/sum(as.numeric(d2$Freq))
    
    d = rbind(d1, d2)
    
    d$Sp = names(L)[k]
    d$Area = area
    d$Habitat_Category = as.factor(d$Habitat_Category)
    d$Period = as.factor(d$Period)

    df = rbind(df, d) 
    
    highx = which(x$mean1 > 666)
    medx = which(x$mean1 <= 666 & x$mean1 >= 333 )
    lowx = which(x$mean1 < 333 )

    ## percentage change
    ct[1,3] = (length(which(x$mean2[highx] > 666))/length(highx))*100
    ct[2,3] = (length(which(x$mean2[highx] <= 666 & x$mean2[highx] >= 333))/length(highx))*100
    ct[3,3] = (length(which(x$mean2[highx] < 333))/length(highx))*100

    ct[4,3] = (length(which(x$mean2[medx] > 666))/length(medx))*100
    ct[5,3] = (length(which(x$mean2[medx] <= 666 & x$mean2[medx] >= 333))/length(medx))*100
    ct[6,3] = (length(which(x$mean2[medx] < 333))/length(medx))*100

    ct[7,3] = (length(which(x$mean2[lowx] > 666))/length(lowx))*100
    ct[8,3] = (length(which(x$mean2[lowx] <= 666 & x$mean2[lowx] >= 333))/length(lowx))*100
    ct[9,3] = (length(which(x$mean2[lowx] < 333))/length(lowx))*100

    ct$Sp = names(L)[k]
    ct$Area = area

    cc = rbind(cc, ct)
    
  }
  
}

df = df %>% group_by(Period) %>% mutate(Subject = row_number())

levels(df$Habitat_Category) <- rev(levels(df$Habitat_Category))

df$Habitat_Category <- factor(df$Habitat_Category, levels = c("High", "Med", "Low"))

png('/Users/Kisei/Desktop/habitat_prop.png', height = 6, width = 6, res = 100, units = "in")
ggplot(df,
       aes(x = Period, stratum = Habitat_Category, alluvium = Subject,
           y = Prop,
           fill = Habitat_Category, label = Habitat_Category)) +
  scale_x_discrete(expand = c(.2, .2)) +
  geom_flow() +
  geom_stratum(alpha = .5)+
  geom_text(stat = "stratum", size = 3) +
  # theme_pubr() + 
  facet_wrap(~Sp + Area, dir = "v")+ 
  theme(legend.position = "none")
dev.off()
cc
