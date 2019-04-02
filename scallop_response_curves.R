library(gridGraphics)
library(lemon)
library(cowplot)

setwd("/Users/Kisei/Google Drive/R/Biomod/scallop/Fancy Forecast")
load("scallop_final_sdms_tuned.RData")

models = unique(mw$models)
models <- gsub(x = models, pattern = "MAXENT1", replacement = "MAXENT.Phillips")
models = models[c(1:2, 4:5, 7)]

png("/Users/Kisei/Desktop/response_ggplotS2.png", width = 15, height = 10, units = "in", res = 100)
# 
# models = c(
#   'GLM',
#   'GAM',
#   'GBM',
#   'CTA',
#   'ANN',
#   # 'SRE',
#   # 'FDA',
#   'MARS',
#   'RF',
#   'MAXENT.Tsuruoka',
#   'MAXENT.Phillips')

# temp --------------------------------------------------------------------
#rows = c(2:5) [,rows]

df_total = data.frame()

for (i in 1:length(models)){
  
  sdm = read.csv(paste0(models[i], "L1.csv"))
  sdm = cbind(sdm,apply(sdm[,3:5],1, median))
  names(sdm)[6] = "median"
  sdm$Model = models[i]
  sdm = sdm[,c(2,6:7)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  
  df_total = rbind(df_total, sdm)
}

library(ggplot2)
p1 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Bottom Temperature (deg C)") +
  ylab("Prob of Presence") +
  theme_classic() + ggtitle("Sea Scallop") + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(size=20), 
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1.5,1,1),"cm")) #top, right, bottom, left

# salinity --------------------------------------------------------------------
#rows = c(6:9)

df_total = data.frame()

for (i in 1:length(models)){
  
  sdm = read.csv(paste0(models[i], "L2.csv"))
  sdm = cbind(sdm,apply(sdm[,3:5],1, median))
  names(sdm)[6] = "median"
  sdm$Model = models[i]
  sdm = sdm[,c(2,6:7)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  
  df_total = rbind(df_total, sdm)
}

p2 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Bottom Salinity") +
  ylab("Prob of Presence") +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(size=20), 
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1.5,1,1),"cm")) #top, right, bottom, left

# depth --------------------------------------------------------------------
rows = c(10:13)

df_total = data.frame()

for (i in 1:length(models)){
  
  sdm = read.csv(paste0(models[i], "L3.csv"))
  sdm = cbind(sdm,apply(sdm[,3:5],1, median))
  names(sdm)[6] = "median"
  sdm$Model = models[i]
  sdm = sdm[,c(2,6:7)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  
  df_total = rbind(df_total, sdm)
}

p3 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Depth (m)") +
  ylab("Prob of Presence") +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(size=20), 
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1.5,1,1),"cm")) #top, right, bottom, left

# lat --------------------------------------------------------------------
rows = c(14:17)

df_total = data.frame()

for (i in 1:length(models)){
  
  sdm = read.csv(paste0(models[i], "L4.csv"))
  sdm = cbind(sdm,apply(sdm[,3:5],1, median))
  names(sdm)[6] = "median"
  sdm$Model = models[i]
  sdm = sdm[,c(2,6:7)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  
  df_total = rbind(df_total, sdm)
}

p4 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Latitude") +
  ylab("Prob of Presence") +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(size=20), 
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1.5,1,1),"cm")) #top, right, bottom, left

# lon --------------------------------------------------------------------


df_total = data.frame()

for (i in 1:length(models)){
  
  sdm = read.csv(paste0(models[i], "L5.csv"))
  sdm = cbind(sdm,apply(sdm[,3:5],1, median))
  names(sdm)[6] = "median"
  sdm$Model = models[i]
  sdm = sdm[,c(2,6:7)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  
  df_total = rbind(df_total, sdm)
}

p5 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Longitude") +
  ylab("Prob of Presence") +
  theme_classic() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(size=20), 
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=unit(c(1,1.5,1,1),"cm")) #top, right, bottom, left

p6 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + theme(legend.text=element_text(size=25))

legend <- g_legend(p6) 


p = cowplot::plot_grid(p1, p2, p3, p4, p5, legend, align = "h", ncol = 3, nrow = 2, hjust = -1)



print(p)

dev.off()
