load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_fall.RData"); ls = xy; ls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_spring.RData"); lf = xy; lf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Scallop_Area_annual.RData"); sa = xy; sa$Species = "Scallop_Annual"

load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_Mu_fall.RData"); mls = mu; mls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_Mu_spring.RData"); mlf = mu; mlf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Scallop_Area_Mu_annual.RData"); msa = mu; msa$Species = "Scallop_Annual"

load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_Prob_fall.RData"); pls = prob; pls$Species = "Lobster_Fall"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Lobster_Area_Prob_spring.RData"); plf = prob; plf$Species = "Lobster_Spring"
load("C:/Users/Kisei/Dropbox/lobster_scallop_manuscript/3rd submission/Scallop_Area_Prob_annual.RData"); psa = prob; psa$Species = "Scallop_Annual"

xy = rbind(lf, ls, sa)
mu = rbind(mlf, mls, msa)
prob = rbind(plf, pls, psa)

xy$period <- relevel(xy$period, "First_10_years")


pdf("/Users/Kisei/Desktop/Figure_4.pdf", width = 7, height = 7) 
ggplot(xy, aes(x = value, fill = period, color = period)) +
  geom_histogram(aes(y =..count..), position = "identity", alpha = 0.5, bins = 50) +
  # geom_histogram(aes(y=..count../sum(..count..)), position = "identity", alpha = 0.5, bins = 50) +
  facet_grid(Area~ Species, scales = "free") +
  geom_vline(data = mu, aes(xintercept = median, color = period),
             linetype = "dashed",
             size = 0.5) +
  scale_color_discrete("") + 
  scale_fill_discrete("") +
  scale_x_continuous(limits = c(0,1), expand = c(0.05, 0.05)) + 
  xlab("Habitat Suitability") +
  ylab("Count") +
  theme_pubr(I(10)) + 
  theme(legend.position = "bottom", 
        legend.justification = c(1,0)) +
  geom_label(data = prob, 
             aes(label= paste0("First_10_yrs=", round(prob1, 2), "\nLast_10_yrs=", round(prob2, 2))),
             x = Inf, y = Inf, 
             hjust = 1, vjust = 1, 
             label.size = NA, size = 3, alpha = 0.5,
             inherit.aes = FALSE) + 
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)))
dev.off()
