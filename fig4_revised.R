load("C:/Users/Kisei/Desktop/Lobster_Area_fall.RData"); ls = xy; ls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Desktop/Lobster_Area_spring.RData"); lf = xy; lf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Desktop/Scallop_Area_annual.RData"); sa = xy; sa$Species = "Scallop_Annual"

load("C:/Users/Kisei/Desktop/Lobster_Area_Mu_fall.RData"); mls = mu; mls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Desktop/Lobster_Area_Mu_spring.RData"); mlf = mu; mlf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Desktop/Scallop_Area_Mu_annual.RData"); msa = mu; msa$Species = "Scallop_Annual"

load("C:/Users/Kisei/Desktop/Lobster_Area_Prob_fall.RData"); pls = prob; pls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Desktop/Lobster_Area_Prob_spring.RData"); plf = prob; plf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Desktop/Scallop_Area_Prob_annual.RData"); psa = prob; psa$Species = "Scallop_Annual"

xy = rbind(lf, ls, sa)
mu = rbind(mlf, mls, msa)
prob = rbind(plf, pls, psa)

png("/Users/Kisei/Desktop/Area.png", width = 7, height = 7, res = 500, units = "in") 
ggplot(xy, aes(x = value, fill = period, color = period)) +
  geom_histogram(aes(y =..count..), position = "identity", alpha = 0.5, bins = 50) +
  facet_grid(Area ~ Species, scales = "free_y") +
  geom_vline(data = mu, aes(xintercept = median, color = period),
             linetype = "dashed",
             size = 0.5) +
  scale_color_discrete("") + 
  scale_fill_discrete("") +
  scale_x_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  xlab("Habitat Suitability") +
  ylab(" ") +
  theme_pubr(I(10)) + 
  theme(legend.position = "bottom", 
        legend.justification = c(1,0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_label(data = prob, 
             aes(label= paste0("First_10_yrs=", round(prob1, 2), "\nLast_10_yrs=", round(prob2, 2))),
             x = Inf, y = Inf, 
             hjust = 1, vjust = 1, 
             label.size = NA, size = 3, alpha = 0.5,
             inherit.aes = FALSE)
dev.off()