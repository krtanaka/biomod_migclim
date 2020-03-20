#open mapping biomod lobster and scallop, load spring, fall, annual SpatialPixelsDataFrame first

<<<<<<< HEAD
annual$Sea_Scallop_Annual = annual$var1.pred*80
annual$American_Lobster_Fall = fall$var1.pred*80
annual$American_Lobster_Spring = spring$var1.pred*80
=======
annual$Sea_Scallop_Annual = annual$var1.pred
annual$American_Lobster_Fall = fall$var1.pred
annual$American_Lobster_Spring = spring$var1.pred
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e

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

<<<<<<< HEAD
max = max(abs(annual@data), na.rm = T)* 120
min = max*-1

png("/Users/Kisei/Desktop/habitat_change.png", width = 18, height = 6, res = 100, units = "in")
=======
max = max(abs(annual@data), na.rm = T)*100
min = max*-1

png("/Users/Kisei/Desktop/habitat_change.png", width = 18, height = 6, res = 500, units = "in")
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
spplot(annual, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/100, #for slope
       col.regions = col,
       par.settings=list(fontsize=list(text=25)),
       zcol = c("American_Lobster_Fall", "American_Lobster_Spring", "Sea_Scallop_Annual"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()
