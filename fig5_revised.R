load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/lobster_change_fall.RData"); ls = data; ls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/lobster_change_spring.RData"); lf = data; lf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/scallop_change_annual.RData"); sa = data; sa$Species = "Scallop_Annual"

colnames(ls) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')
colnames(lf) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')
colnames(sa) = c('x', 'y', 'Management_Area', 'Model', 'Lobster_Scallop', 'Species')

data = rbind(lf, ls, sa)

pdf("/Users/Kisei/Desktop/Figure_5.pdf", width = 6, height = 6) 
ggplot(data, aes(x = x, y = y, color = Model)) +
  geom_line( data = subset( data, !(Model %in% 'Ensemble') ),
             size = I( 0.2 )) +
  geom_line( data = subset( data, Model %in% 'Ensemble' ),
             color = I( 'black' ),
             show.legend = FALSE,
             size = I( 1 )) +
  theme_pubr(base_size = I(10)) +
  scale_color_discrete("") + 
  theme(legend.position = "none", 
        legend.justification = c(1,0)) + 
  facet_wrap(.~ Species + Management_Area, scales = "free_y", dir = "v") +
  # facet_grid(Management_Area ~ Species, scales = "free") +
  labs(x = "Model Year", y = "Habitat Suitability")
dev.off()
