library(reshape2)
library(ggpubr)
library(ggplot2)

dir = "/Users/Kisei/Google Drive/R/Biomod/"

species = c("lobster", "scallop")
op = c("tuned", "default")[1]

var = NULL

for (i in 1:length(species)) {
  
  setwd(paste0("/Users/Kisei/Google Drive/R/Biomod/", species[[i]]))
  load("final_sdms.RData")
  
  if (op == "default") {
    
    load("Biomod_Results_with_Default_Options.RData")
    
  }  else {
    
    load("Biomod_Results.RData")
    
  }
  
  df = myBiomodModelOut@variables.importances@val
  df = as.data.frame(df)
  names(df) <- sub(".AllData", "", names(df))
  rownames(df) = c("Bottom Temperature", "Bottom Salinity", "Depth", "Latitude", "Longitude")
  
  if (i == 1) {
    final_sdms = gsub("lobster_AllData_", "", final_sdms)
    
  } else{
    final_sdms = gsub("scallop_AllData_", "", final_sdms)
    
  }
  final_sdm_1 = sapply(strsplit(final_sdms,"_"), `[`, 1)
  final_sdm_2 = sapply(strsplit(final_sdms,"_"), `[`, 2)
  final_sdm = paste0(final_sdm_2, ".", final_sdm_1)
  
  # df <- df[, final_sdm] #variable importance based on final or all models
  
  df = as.data.frame(t(df))
  df$Species = species[i]
  
  var = rbind(var, df)
  
  
}

var <- melt(var,id.vars='Species', measure.vars = c("Bottom Temperature", "Bottom Salinity", "Depth", "Latitude", "Longitude"))

# png("/Users/kisei/Desktop/biomod_var_contribution.png", height = 8, width = 10, units = "in", res = 500)
p = ggplot(var) +
  geom_boxplot(aes(x=variable, y=value, fill=Species)) + 
  xlab("") + ylab("Variable Importancece Score") + 
  theme_pubr(I(20)) +  
  theme(
    legend.position = c(0,1),
    legend.justification = c(0, 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
# dev.off()
