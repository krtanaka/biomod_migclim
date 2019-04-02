# Load the saved Biomod results, compare results

rm(list = ls())

library(biomod2)
library(raster)
library(readr)
library(colorRamps)
library(automap)
library(beepr)
library(gstat)
library(cowplot)
library(maps)
library(Rmisc)
library(svMisc)
library(pals)
library(doParallel)
library(data.table)
library(grid)
library(ENMeval)
library(ggpubr)

dir = paste0("/Users/", Sys.info()[7], "/")

sp = c("lobster", "scallop")[2]

setwd(paste0(dir, "/Desktop/", sp)) #with slope
setwd(paste0(dir, "/biomod_migclim/", sp)) #without slope

op = c("tuned", "default")[1]

if (op == "default") {
  
  load("Biomod_Results_with_Default_Options.RData")
  
}  else {
  
  load("Biomod_Results.RData")
  
}

# Get TSS and ROC summaries
# summarize TSS and ROC scores of all models
tss = myBiomodModelEval["TSS","Testing.data",,,]; tss = as.data.frame(tss)
roc = myBiomodModelEval["ROC","Testing.data",,,]; roc = as.data.frame(roc)

# roc$mean = apply(roc, 1, mean)
# tss$mean = apply(tss, 1, mean)

get_variables_importance(myBiomodModelOut)

models = c('GLM',
           'GAM',
           'GBM',
           'CTA',
           'ANN',
           'SRE',
           'FDA',
           'MARS',
           'RF',
           'MAXENT.Tsuruoka',
           'MAXENT.Phillips')

tss = cbind(models, tss)
colnames(tss) = c("SDM", "TSS_RUN1","TSS_RUN2","TSS_RUN3")
roc = cbind(models, roc)
colnames(roc) = c("SDM", "ROC_RUN1","ROC_RUN2","ROC_RUN3")
tss_roc = merge(tss, roc)
tss_roc = subset(tss_roc, SDM!="MAXENT.Tsuruoka")
tss_roc$SDM = gsub("MAXENT.Phillips", "MAXENT", tss_roc$SDM)
tss_roc
# write_csv(tss_roc, paste0(dir, "Desktop/", sp, "_ROC_TSS.csv"))

# Compare model performance and select final models, both lobster and scallop
# organize all SDMs by TSS and ROC
model_eval = t(as.data.frame(myBiomodModelEval))
seq = seq(1, nrow(model_eval), by = 4)
model_eval = model_eval[seq,]
model_eval = as.data.frame(model_eval)
model_eval$combo = model_eval$TSS + model_eval$ROC
model_eval=model_eval[order(-model_eval$combo),]
row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Phillips", replacement = "MAXENT1")  
row.names(model_eval) <- gsub(x = row.names(model_eval), pattern = "MAXENT.Tsuruoka", replacement = "MAXENT2")  

sdm = data.frame(model = rep(NA,nrow(model_eval)),run = rep(NA,nrow(model_eval)))

for (j in 1:nrow(model_eval)){
  ff = strsplit(row.names(model_eval[j,]), '[._]')[[1]]
  sdm[j,1] = ff[3]
  sdm[j,2] = gsub("RUN", "", ff[length(ff)-1]) 
} 

evals = cbind(sdm,model_eval);rownames(evals) <- c()

d1 = summarySE(evals, measurevar = "TSS", groupvars = c("model"))
d2 = summarySE(evals, measurevar = "ROC", groupvars = c("model"))

eval = merge(d1[c("model","TSS","se")],d2[c("model","ROC","se")], by = "model")
names(eval) = c("Model", "TSS", "SE2", "ROC", "SE1")
eval = subset(eval, Model!="MAXENT2")
eval$Model = gsub("MAXENT1", "MAXENT", eval$Model)

# jpeg(paste0("/Users/Kisei/Desktop/model_evaluation_legend.jpg"), res = 500, height = 2.5, width = 1, units = "in")
# legend = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
#   geom_point(size=4) +
#   scale_shape_manual(values = rep(18, 11))+
#   geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
#   geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) +
#   ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95)) +
#   theme(legend.position = "right", legend.title = element_blank())
# 
# legend = cowplot::get_legend(legend)
# grid.newpage()
# grid.draw(legend)
# dev.off()

png(paste0(dir, "/Desktop/model_evaluation_", sp, ".png"), res = 500, height = 5, width = 3, units = "in")

if (sp == "lobster") {
  
  p = ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
    geom_point(size=4) +
    scale_shape_manual(values = rep(18, 11))+
    geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
    geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
    ggtitle("American lobster")+
    ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+ 
    theme_classic() + 
    theme(legend.position="none")    
}

if (sp == "scallop") {
  
  p =  ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, shape = Model)) +
    geom_point(size=4) +
    scale_shape_manual(values = rep(18, 11))+
    geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
    geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) + 
    ggtitle("Sea scallop")+
    ylim(c(0.2,0.7)) + xlim(c(0.6, 0.95))+
    theme_classic() + 
    theme(legend.position="none")
  
}

print(p)

# ggplot(data = eval,aes(x = ROC, y = TSS, colour = Model, label = Model)) +
#   geom_point(size=2) + geom_text(aes(label = Model), hjust = -0.1, vjust = -0.8, size = 5) +
#   scale_shape_manual(values = rep(18, 11))+
#   geom_errorbar(aes(ymin = TSS-SE2, ymax = TSS + SE2)) +
#   geom_errorbarh(aes(xmin = ROC-SE1, xmax = ROC + SE1)) +
#   ggtitle("") +
#   xlim(0.5, 1) +
#   ylim(0.3, 0.8) +
#   theme_classic(base_size = I(20)) +
#   theme(legend.position="none")

dev.off()

#remove SDMs with TSS < 0.5
evals = subset(evals, TSS > 0.5)
evals = subset(evals, ROC > 0.8)

# weighting all models by TSS, ROC, and combined weight
weight = evals[c(1:nrow(evals)),5] #4 is ROC, 5 is combo
scale <- function(x){(x-min(x))/(max(x)-min(x))} #scale weight to 0-1
weight = scale(weight)
models = evals$model
run = evals$run
mw = data.frame(models, run, weight)
mw$SDM = paste0(mw$models, "_", mw$run)

# jpeg(paste0("/Users/Kisei/Desktop/model_weights_legend.jpg"), res = 500, height = 5, width = 1, units = "in")
# legend = ggplot(mw, aes(x=weight, y=SDM, color=models)) + 
#   geom_point(aes(size = weight)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# legend = cowplot::get_legend(legend)
# grid::grid.newpage()
# grid::grid.draw(legend)
# dev.off()

jpeg(paste0(dir, "Desktop/model_weights_", sp, ".jpg"), res = 500, height = 4, width = 8, units = "in")

colnames(mw) = c("SDM", "run", "Weight", "SDM_runs")

if (sp == "lobster") {
  p = ggplot(mw, aes(x=Weight, y=SDM_runs, color=SDM)) + 
    # geom_point(size = 5) + 
    geom_point(aes(size = Weight)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("American lobster")+ 
    theme_classic() + 
    theme(legend.position="right")  
}

if (sp == "scallop") {
  p = ggplot(mw, aes(x=Weight, y=SDM_runs, color=SDM)) + 
    # geom_point(size = 5) + 
    geom_point(aes(size = Weight)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Sea scallop")+
    theme_classic() + 
    theme(legend.position="right")  
}

print(p)
dev.off()

# good_sdms = models with TSS > 0.5
a = paste0(sp, "_AllData_RUN")
good_sdms = paste0(a, mw$run, "_", mw$models)
good_sdms = gsub("MAXENT1", "MAXENT.Phillips", good_sdms)
good_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", good_sdms)

d = myBiomodProj@proj@val@layers
dd = stack(d)
ddd = subset(dd,good_sdms)

#all sdm > 0.5 TSS
# jpeg(paste0("/Users/Kisei/Desktop/good_sdms_", sp, ".jpg"), res = 500, height = 10, width = 20, units = "in")
plot(ddd, col = matlab.like(100), nc = 6)
# dev.off()

# create final_sdms = visually check each model output and remove any unrealistic models

if (sp == "lobster") {
  mw = mw[mw$models!="CTA",]
  # mw = mw[mw$models!="MAXENT1",]
  mw = mw[mw$models!="GBM",]
}

if (sp == "scallop") {
  mw = mw[mw$models!="CTA",]
  # mw = mw[mw$models!="MAXENT1",]
  mw = mw[mw$models!="GBM",]
}

final_sdms = paste0(a, mw$run, "_", mw$models)
final_sdms = gsub("MAXENT1", "MAXENT.Phillips", final_sdms)
final_sdms = gsub("MAXENT2", "MAXENT.Tsuruoka", final_sdms)
d = myBiomodProj@proj@val@layers
dd = stack(d)
ddd = subset(dd,final_sdms)

#do weighted mean
avg = weighted.mean(ddd, mw$weight); dev.off(); 
avg@data@values = (avg@data@values - avg@data@min)/(avg@data@max - avg@data@min)

# df = as.data.frame(rasterToPoints(avg))
# load(paste0(dir, "biomod_migclim/lobster/lobster_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data
# obs = lobster[,c("Lon", "Lat", "expcatnum")]
# 
# load(paste0(dir, "biomod_migclim/scallop/scallop_survey_data_spring_fall_combined_1984-2016.RData")) #load survey data
# scallop$year = substr(scallop$cruise6, 1, 4); scallop = subset(scallop, year %in% c(1984:2016))
# 
# obs = scallop[,c(19, 20, 16)] 
# 
# obs$expcatnum = obs$expcatnum/(max(obs$expcatnum)-min(obs$expcatnum))
# 
# par(mfrow = c(2,1), mar = c(3,3,3,5))
# plot(df$x, df$layer, col=4, pch = 20)
# points(obs$lon, obs$expcatnum, col = 2, pch = 20)
# 
# plot(df$y, df$layer, col=4, pch = 20)
# points(obs$lat, obs$expcatnum, col = 2, pch = 20)
# 
# df$x = round(df$x, 2)
# df$y = round(df$y, 2)
# obs$Lon = round(obs$Lon, 2)
# obs$Lat = round(obs$Lat, 2)
# 
# colnames(df) = c("lon", "lat", "pred")
# colnames(obs) = c("lon", "lat", "obs")
# 
# df = merge(obs, df)
# cor(df$obs, df$pred)


# jpeg(paste0(dir, "Desktop/Biomod_Parsimonious_Model_Prediction", sp, ".jpg"), res = 500, height = 5, width = 10, units = "in")

par(mfrow = c(1,2), mar = c(3,3,3,5))
plot(avg, zlim = c(0,1), col = matlab.like(100), 
     main = ifelse(sp == "lobster", "American lobster", "Sea scallop"),
     xaxt='n', yaxt='n', 
     legend = F)
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
box()
degAxis(1, cex.axis = 1)
degAxis(2, cex.axis = 1, las = 1)
legend("bottomright", "Mean", bty = "n", cex = 1.4)

#plot standard error from final SDMs
plot(calc(ddd, fun=sd)*0.001, 
     col = matlab.like(100), 
     zlim = c(0,1),
     axes = F)
map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
box()
degAxis(1, cex.axis = 1)
degAxis(2, cex.axis = 1, las = 2)
legend("bottomright", "Standard deviation", bty = "n", cex = 1.4)

# dev.off()

#save list of final SDMs and weights.
# if (op == "default") {
#   
#   save(final_sdms, mw, file = "final_sdms_default.RData")
#   
# }else{
#   
#   save(final_sdms, mw, file = "final_sdms.RData")
#   
# }



