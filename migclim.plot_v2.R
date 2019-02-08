migclim.plot = function () {
  
  # sim = paste0(par1, "_", par2, "_", par3, "_", par4, "_", par5, "_", par6) 
  sim = paste0(par1, "_", par2, "_", par3) 
  
  
  all_raster = list.files(getwd(), pattern =".asc$", full.names=TRUE)
  
  all_raster = stack(all_raster)
  
  Rst = mean(all_raster)
  
  rstVals <- sort(raster::unique(Rst))
  negativeNb <- length(which(rstVals < 0))
  positiveNb <- length(which(rstVals > 1 & rstVals < 30000))
  zeroExists <- any(rstVals == 0)
  oneExists <- any(rstVals == 1)
  unilimtedExists <- any(rstVals == 30000)
  
  # Colors <- rep("yellow", negativeNb)
  # if (zeroExists) Colors <- c(Colors, "grey94")
  # if (oneExists) Colors <- c(Colors, "black")
  # Colors <- c(Colors, rainbow(positiveNb, start = 0, end = 0.4))
  # if (unilimtedExists) Colors <- c(Colors, "pink")
  
  Colors = parula(length(rstVals))
  # Colors[1:negativeNb] <- rep("yellow", negativeNb) #replace first color scheme to represent negativeNB
  if (zeroExists) Colors[negativeNb+1] <- "white" #2nd color scheme
  if (oneExists) Colors[negativeNb+2] <- "black" #3rd color scheme
  if (unilimtedExists) Colors[length(rstVals)] <- "white"
  
  setwd(paste0(dir, "/Google Drive/R/Biomod/lobster"))
  
  image(Rst, col = Colors, breaks = c(min(rstVals) - 1, rstVals), xaxt='n', yaxt='n', xlab = "", ylab = "")
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1)
  
  legend("bottomright",
         c("Immediate colonization","","","","Late colonization"),
         fill = Colors[c(seq(3, length(Colors), by = length(Colors)/5))], 
         # xpd = NA,
         # cex = 1.5,
         bty = "n")
  
  legend("topleft",
         c(paste0("barrierType = ", strsplit(sim,"_")[[1]][1]),
           paste0("dispKernel = ", strsplit(sim,"_")[[1]][2]),
           paste0("iniMatAge = ", strsplit(sim,"_")[[1]][3]),
           paste0("lddFreq = ", strsplit(sim,"_")[[1]][4]),
           paste0("lddMinDist = ", strsplit(sim,"_")[[1]][5]),
           paste0("lddMaxDist = ", strsplit(sim,"_")[[1]][6])),
         # xpd = NA,
         cex = 1.5,
         bty = "n")
  
  jpeg(filename = paste0("Migclim_", sim, Sys.Date(), ".jpg"), units = "in", width = 10, height = 10, res = 500)
  
  image(Rst, col = Colors, breaks = c(min(rstVals) - 1, rstVals), xaxt='n', yaxt='n', xlab = "", ylab = "")
  map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
  box()
  degAxis(1, cex.axis = 1)
  degAxis(2, cex.axis = 1)
  
  legend("bottomright",
         c("Immidiete colonization","","","","Late colonization"),
         fill = Colors[c(seq(3, length(Colors), by = length(Colors)/5))], 
         # xpd = NA,
         # cex = 1.5,
         bty = "n")
  
  legend("topleft",
         c(paste0("barrierType = ", strsplit(sim,"_")[[1]][1]),
           paste0("dispKernel = ", strsplit(sim,"_")[[1]][2]),
           paste0("iniMatAge = ", strsplit(sim,"_")[[1]][3]),
           paste0("lddFreq = ", strsplit(sim,"_")[[1]][4]),
           paste0("lddMinDist = ", strsplit(sim,"_")[[1]][5]),
           paste0("lddMaxDist = ", strsplit(sim,"_")[[1]][6])),
         # xpd = NA,
         cex = 1.5,
         bty = "n")

  dev.off()
  
  rm(Rst, rstVals, negativeNb, positiveNb, zeroExists, oneExists, unilimtedExists, Colors)
}
