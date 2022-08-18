

require(raster)
require(maptools)
require(ggplot2)
require(pals)

setwd("C:/Users/Mike/Documents/NY Bight Project") 

COAST = readOGR("Environmental Data/US_Coastlines/US_Coast.shp")

POS = read.csv(file = 'POS_scallop_tuned4.csv', header = TRUE, row.names=NULL) 
names(POS)[2:3] = c('longitude','latitude')
load('cPOS.RData')
CS = read.csv(file = 'CurrentScallop.csv', header = TRUE, row.names=NULL) 

POS = POS[,c(2:length(POS))]
CS = CS[,c(2,3,4)]

rm(list= ls()[!(ls() %in% c('POS','CS', 'COAST','avg' ))])

### Make a grid ###

xrange = range(CS$longitude)
yrange = range(CS$latitude)

grid= expand.grid(lon = seq(from = xrange[1], to = xrange[2], by = 0.01), 
                  lat = seq(from = yrange[1], to = yrange[2], by = 0.01))

coordinates(grid)= ~ lon+ lat

CS <- rasterFromXYZ(CS[, c("longitude", "latitude", "var1.pred")])

CSval =extract(CS, grid)

POS1 = cbind(as.data.frame(grid),CSval)

avgval =extract(avg, grid)

POS1 = cbind(POS1,avgval)


P1 <- rasterFromXYZ(POS1[, c('lon','lat', 'CSval')])

P2 = rasterFromXYZ(POS1[, c('lon','lat', 'avgval')])

 plot(P2)
 plot(P1)

for (i in 1:(length(POS)-2)){
   
   P1 <- rasterFromXYZ(POS[, c('longitude','latitude', paste("hsmap",i,sep=""))])
   P2 = extract(P1, grid)
   
   POS1 = cbind(POS1,P2)
   names(POS1)[i+4] = paste("hsmap",i,sep="")
   
   print(i)
}

 for (i in 5:length(POS1)){
    POS1[,i] = POS1[,i]/ POS1[,4]
 }
 
 
 

 ## Projection on current distribution plot
  
 for (i in 1:(length(POS1)-2)){

    x = as.data.frame(cbind(POS1$lon, POS1$lat, (POS1$CSval*POS1[,i+4])))
    
     CSRas = rasterFromXYZ(x[, c('V1','V2', 'V3')])
 
 
 

 
  png("C:/Users/Mike/Documents/NY Bight Project/Projections/",file = paste("hsmap",i,".png",sep=""),
      width = 6, height = 6, units = "in", res = 200)
  brk=seq(0,35,length.out=100)
  mypal <- colorRampPalette(c("black", "blue", "purple","yellow"), bias=1)
  plot(CSRas)
  plot(COAST, col='grey', add=TRUE)
  dev.off()
 
  print(i)
 }
 

 ## projection  plots 
 for (i in 1:(length(POS)-2)){
    
    
    hsi = rasterFromXYZ(POS[, c('longitude','latitude', paste("hsmap",i,sep=""))])
    
    
    
    png("C:/Users/Mike/Documents/NY Bight Project/Projections/",file = paste("hsi",i,".png",sep=""),
        width = 6, height = 6, units = "in", res = 200)
    brk=seq(0,1000,length.out=100)
    mypal <- colorRampPalette(c("black", "blue", "purple","yellow"), bias=1)
    plot((hsi), breaks=brk,col= parula(100))
    plot(COAST, col='grey', add=TRUE)
    dev.off()
    
    print(i)
 }
 
 
 
 ## Environment  plots 
 for (i in 1:(length(POS)-2)){
    
    
    hsi = rasterFromXYZ(POS[, c('longitude','latitude', paste("hsmap",i,sep=""))])
    
    
    
    png("C:/Users/Mike/Documents/NY Bight Project/Projections/",file = paste("hsi",i,".png",sep=""),
        width = 6, height = 6, units = "in", res = 200)
    brk=seq(0,1000,length.out=100)
    mypal <- colorRampPalette(c("black", "blue", "purple","yellow"), bias=1)
    plot((hsi), breaks=brk,col= parula(100))
    plot(COAST, col='grey', add=TRUE)
    dev.off()
    
    print(i)
 }
 
 
 
 
 
 ##Change Plot
 
 POS80 = rowMeans(POS[, 45:49])
 POS10 = rowMeans(POS[, 3:7])
 
 
 POSD = POS80 - POS10
 
 
 
 f = as.data.frame(cbind(POS$longitude, POS$latitude, POSD))
 
 CSRas = rasterFromXYZ(f[, c('V1','V2', 'POSD')])
 
 
 
 mypal <- colorRampPalette(c("blue",
                             MixColor('blue', 'white', amount1 = 0.5),
                             "white",
                             MixColor('red', 'white', amount1 = 0.5),
                             "red"), 
                           bias=0.75)
 
 
 png("C:/Users/Mike/Documents/NY Bight Project",file = 'scaldelta.png',
     width = 6, height = 6, units = "in", res = 200)
 plot(CSRas, col=mypal(100))
 plot(COAST, col='grey', add=TRUE)
 dev.off()
 
 
 plot(hsi)
 
 
avg = resample(avg, CS)

POS1 = POS1[,c(2,3,4)]
POS1 = cbind(POS1, as.data.frame(avg))

# png("C:/Users/Mike/Documents/NY Bight Project/shapes.png",
#     width = 6, height = 6, units = "in", res = 600)
# plot(avg, col='green')
# plot(CS, col='blue', add=TRUE)
# plot(COAST, col='grey', add=TRUE)
# dev.off




for (i in 5:length(POS1)){
  POS1[,i] = POS1[,i]/ POS1[,4]
}

POS2 = POS1
for (i in 5:length(POS1)){
  POS2[,i] = POS1[,i]* POS1[,3]
}


p3 = rasterFromXYZ(POS1[, c('longitude','latitude', paste("hsmap",i,sep=""))])
p3 <- flip(p3, direction='y')

plot(p3)


POS2 = POS1[!is.na(POS1$var1.pred) & !is.na(POS1$layer)
            & !is.na(POS1$hsmap33) & !is.na(POS1$hsmap49), ]




P2 <- rasterFromXYZ(POS1[, c("longitude", "latitude", "hsmap1")])



CSc = crop(CS,P1)








cLat = resample(cLat, CS)
cLon = resample(cLon)

avg = as.data.frame(avg)
POS1 = as.data.frame(POS1)


r3 <- overlay(r1, r2, fun=function(x,y){return(x+y)})








   
   
   
   
   r3 <- overlay(r1, r2, fun=sum)

