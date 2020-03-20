load("C:/Users/Kisei/Desktop/lobster_change_fall.RData"); ls = data; ls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Desktop/lobster_change_spring.RData"); lf = data; lf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Desktop/scallop_change_annual.RData"); sa = data; sa$Species = "Scallop_Annual"

colnames(ls) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')
colnames(lf) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')
colnames(sa) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')

data = rbind(lf, ls, sa)

png("/Users/Kisei/Desktop/Area.png", width = 7, height = 7, res = 500, units = "in") 
ggplot(data, aes(x = x, y = y, color = Model)) +
  geom_line( data = subset( data, !(Model %in% 'Ensemble') ),
             size = I( 0.2 )) +
  geom_line( data = subset( data, Model %in% 'Ensemble' ),
             color = I( 'black' ),
             show.legend = FALSE,
             size = I( 1 )) +
  theme_pubr(base_size = I(10)) +
  scale_color_discrete("") + 
  theme(legend.position = "right", 
        legend.justification = c(1,0)) + 
  facet_wrap(.~ Species + Management_Area, scales = "free_y") +
  # facet_grid(Management_Area ~ Species, scales = "free_y") +
  labs(x = "Model Year", y = "Relative Habitat Suitability")
dev.off()