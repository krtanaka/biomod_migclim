setwd("C:/Users/Kisei/Google Drive/Research/Charlie/VAST/NWA_american_lobster_VAST_output")
load("Save.RData")
Sdreport=Save$Opt[["SD"]]
TmbData = Save$TmbData

Year_Set = c(1981:2017)

if ("ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
  CogName = "mean_Z_tm"
  EffectiveName = "effective_area_tl"
  TmbData[["n_c"]] = 1
}
if ("ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
  CogName = "mean_Z_ctm"
  EffectiveName = "effective_area_ctl"
}
if ("ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport))) {
  CogName = "mean_Z_cym"
  EffectiveName = "effective_area_cyl"
  TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
}

category_names = 1:TmbData$n_c

SD = TMB::summary.sdreport(Sdreport)

df = data.frame(SD[which(rownames(SD) == CogName), 
                     c("Estimate", "Std. Error")])
df1 = df[1:length(Year_Set),]
df2 = df[length(Year_Set):max(length(Year_Set)),]

plot()