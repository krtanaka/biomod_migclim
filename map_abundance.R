load("/Users/Kisei/biomod_migclim/lobster/lobster_survey_data_spring_fall_combined_1984-2016.RData") #load survey data

xlims <- range(pretty(lobster$Lon))
ylims <- range(pretty(lobster$Lat))

lobster$CPUE = lobster$expcatnum

l = ggplot()+
  geom_point(data=lobster, aes(x=Lon, y=Lat, color=log10(CPUE+1)), size = 1, alpha = 0.5) +
  borders(xlim=xlims,ylim=ylims, fill = "gray") +
  coord_quickmap(xlim = xlims,ylim=ylims) +
  theme_pubr() + 
  scale_colour_gradientn(colours=matlab.like(100)) + 
  theme(legend.position = "bottom", 
        legend.justification = c(1, 0))

lobster=subset(lobster, Vessel == "CT")
xlims <- range(pretty(lobster$Lon))
ylims <- range(pretty(lobster$Lat))
xlims <- c( -74, -71.5)
ylims <- c(40.75,41.35)==+=

l = ggplot()+
  geom_point(data = lobster, aes(x=Lon, y=Lat, color=log10(CPUE+1)), size = 5, alpha = 0.1) +
  borders(xlim=xlims,ylim=ylims, fill = "gray") +
  coord_quickmap(xlim = xlims,ylim=ylims) +
  theme_pubr(I(20)) +
  scale_colour_gradientn(colours=matlab.like(100)) + 
  theme(legend.position = c(0.88, 0.32))

load(paste0(dir, "biomod_migclim/scallop/scallop_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data
scallop$year = substr(scallop$cruise6, 1, 4); scallop = subset(scallop, year %in% c(1984:2016))
scallop = subset(scallop, lat >= 37)

scallop$n = scallop$expcatnum

s = ggplot()+
  geom_point(data=scallop, aes(x=lon, y=lat, color=log10(n+1)), size = 1, alpha = 0.5) +
  borders(xlim=xlims,ylim=ylims, fill = "gray") +
  coord_quickmap(xlim = xlims,ylim=ylims) +
  scale_colour_gradientn(colours=matlab.like(100)) + 
  theme(legend.position = c(0.7,0.2))


setwd("/Users/ktanaka/Desktop")
jpeg("lobster_abundance.jpg", res = 500, height = 5, width = 11, units = "in")
l
dev.off()

setwd("/Users/Kisei/Desktop")
jpeg("scallop_abundance.jpg", res = 500, height = 5, width = 5, units = "in")
s
dev.off()
