#open mapping biomod lobster and scallop, load spring, fall, annual SpatialPixelsDataFrame first

annual$Sea_Scallop_Annual = annual$var1.pred
annual$American_Lobster_Fall = fall$var1.pred
annual$American_Lobster_Spring = spring$var1.pred

png("/Users/Kisei/Desktop/habitat_mean.png", width = 18, height = 6, res = 500, units = "in")
spplot(annual, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (0:100)/100, #for slope
       col.regions = matlab.like(100),
       par.settings=list(fontsize=list(text=25)),
       zcol = c("American_Lobster_Fall", "American_Lobster_Spring", "Sea_Scallop_Annual"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()

max = max(abs(annual@data), na.rm = T)*12000
min = max*-1

png("/Users/Kisei/Desktop/habitat_change.png", width = 18, height = 6, res = 100, units = "in")
spplot(annual, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/10000, #for slope
       col.regions = col,
       par.settings=list(fontsize=list(text=25)),
       zcol = c("American_Lobster_Fall", "American_Lobster_Spring", "Sea_Scallop_Annual"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()
