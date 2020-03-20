#open mapping biomod lobster and scallop, load spring, fall, annual SpatialPixelsDataFrame first

annual$Sea_Scallop_Annual = annual$var1.pred*80
annual$American_Lobster_Fall = fall$var1.pred*80
annual$American_Lobster_Spring = spring$var1.pred*80

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

png("/Users/Kisei/Desktop/habitat_change.png", width = 18, height = 6, res = 100, units = "in")
max = max(abs(annual@data), na.rm = T)*100
min = max*-1

png("/Users/Kisei/Desktop/habitat_change.png", width = 18, height = 6, res = 500, units = "in")
spplot(annual, 
       sp.layout = list(list("sp.polygons", countriesHigh, lwd=0.1, fill="grey")),
       at = (min:max)/100, #for slope
       col.regions = col,
       par.settings=list(fontsize=list(text=25)),
       zcol = c("American_Lobster_Fall", "American_Lobster_Spring", "Sea_Scallop_Annual"),
       scales=list(draw=T),
       colorkey = T) 
dev.off()

load("C:/Users/Kisei/Desktop/Lobster_spring.RData"); df1 = ks_df
load("C:/Users/Kisei/Desktop/Lobster_fall.RData"); df2 = ks_df
load("C:/Users/Kisei/Desktop/Scallop_annual.RData"); df3 = ks_df
df = rbind(df1, df2, df3)

source("/Users/Kisei/Google Drive/R/misc/color palette function.R")
steps = c("blue", "white", "red")
col = color.palette(steps, space="rgb")

xlims <- range(pretty(df$x));ylims <- range(pretty(df$y))

p1 = ggplot(data = df[which(df$Area %in% c("GOM_GB", "SNE", "GOM_GB Nearshore", "SNE Nearshore")),], 
            aes(x = x, y = y, fill = D)) +
        geom_raster(interpolate = TRUE) +
        scale_fill_gradientn(colours = col(length(unique(df$D)))) + #parura(100)
        facet_grid(Area~Species, scales = "free") + 
        coord_quickmap(xlim = xlims,ylim = ylims) +
        borders(xlim = xlims,ylim = ylims, fill = "gray") +
        theme_pubr(I(10)) + 
        theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1)) + 
        scale_x_longitude(xmin=-180, xmax=180, step=2) +
        scale_y_latitude(ymin=-180, ymax=180, step=2)

p2 = ggplot(data = df[which(df$Area %in% c("GOM_GB", "SNE", "GOM_GB Nearshore", "SNE Nearshore")),], 
            aes(x = x, y = y, fill = P_0.05)) +
        geom_raster(interpolate = TRUE) +
        scale_fill_gradientn(colours = col(length(unique(ks_df$P_0.05)))) + #parura(100)
        facet_grid(Area~Species, scales = "free") + 
        coord_quickmap(xlim = xlims,ylim = ylims) +
        borders(xlim = xlims,ylim = ylims, fill = "gray") +
        theme_pubr(I(10)) + 
        theme(legend.position="right", axis.text.x = element_text(angle = 45, hjust = 1)) + 
        scale_x_longitude(xmin=-180, xmax=180, step=2) +
        scale_y_latitude(ymin=-180, ymax=180, step=2)

png("/Users/Kisei/Desktop/D_P0.05.png", units = "in", res = 500, width = 7, height = 7) #res = 500, width=6000,height=6000
p2
dev.off()
