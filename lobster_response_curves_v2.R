dir = "/Users/Kisei/"
library(ggplot2)

# setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
<<<<<<< HEAD
setwd("/Users/ktanaka//Desktop/lobster/")
=======
setwd("/Users/kisei//Desktop/scallop"))
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e

load("final_sdms.RData")

model = unique(mw$models)
# model <- gsub(x = model, pattern = "MAXENT1", replacement = "MAXENT.Phillips")
model <- gsub(x = model, pattern = "MAXENT1", replacement = "MAXENT")

png(paste0(dir, "Desktop/response_ggplot_", Sys.Date(), ".png"), width = 15, height = 10, units = "in", res = 100)

# temp --------------------------------------------------------------------
rows = c(2:5)

df_total = data.frame()

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:3],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  df_total = rbind(df_total, sdm)
}

p1 = ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Bottom Temperature (deg C)") +
  ylab("Prob of Presence") +
  theme_classic() + ggtitle("American lobster") + 
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
rows = c(6:9)

df_total = data.frame()

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:3],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
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

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:3],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
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

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:3],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
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
rows = c(18:21)

df_total = data.frame()

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:3],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
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

legend <- lemon::g_legend(p6) 


p = cowplot::plot_grid(p1, p2, p3, p4, p5, legend, align = "h", ncol = 3, nrow = 2, hjust = -1)

print(p)

dev.off()


<<<<<<< HEAD

# steepness --------------------------------------------------------------------
rows = c(22:25)

df_total = data.frame()

for (i in 1:length(model)){
  
  sdm = read.csv(paste0("Biomod_", model[i], ".csv"))[,rows]
  sdm = cbind(sdm,apply(sdm[,2:4],1, median))
  names(sdm)[5] = "median"
  sdm$Model = model[i]
  sdm = sdm[,c(1,5:6)]
  colnames(sdm) = c("Var", "Value", "SDM")
  
  df_total = rbind(df_total, sdm)
}

ggplot(df_total, aes(x = Var, y = Value, group = SDM, colour=SDM)) + 
  geom_line(size = 2) + 
  xlab("Steepness") +
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

p = cowplot::plot_grid(p1, p2, p3, p4, p5, legend, align = "h", ncol = 3, nrow = 2, hjust = -1)

print(p)

dev.off()


=======
>>>>>>> c6a54766214084fa3f5a0b3ca60e276f5475b34e
