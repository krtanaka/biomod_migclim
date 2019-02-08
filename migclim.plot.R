migclim.plot = function (asciiFile, 
                         outDir = "", 
                         fileFormat = "jpeg", 
                         fullOutput = FALSE) {
  


  if (substr(asciiFile, nchar(asciiFile) - 3, nchar(asciiFile)) != ".asc") 
    asciiFile <- paste(asciiFile, ".asc", sep = "")
  
  if (fileFormat != "jpeg" & fileFormat != "png" & fileFormat != "inR") 
    stop("Input error: 'fileFormat' must be one of 'jpeg' or 'png'.\n")
  
  if (!file.exists(asciiFile)) 
    stop("Input error: ", asciiFile, " could not be found.\n")
  
  if (outDir != "") 
    if (!file.exists(outDir)) 
      stop("Input error: 'outDir' directory could not be found.\n")
  
  baseName <- substr(basename(asciiFile), 1, nchar(basename(asciiFile)) - 11)
  
  if (outDir == "" & dirname(asciiFile) != ".") 
    outDir <- dirname(asciiFile)
  
  outDir <- paste(outDir, "/", sep = "")
  
  inDir <- paste(dirname(asciiFile), "/", sep = "")
  
  if (inDir == "./") 
    inDir <- ""
  
  if (fullOutput) {
    
    stepName <- paste(baseName, "_step_", sep = "")
    fileList <- list.files(dirname(asciiFile))
    fileList <- fileList[which(substr(fileList, 1, nchar(stepName)) == stepName)]
    fileList <- c(basename(asciiFile), fileList)
    rm(stepName)
    
  } else {
    
    fileList <- basename(asciiFile)
  }
  
  fileList <- substr(fileList, 1, nchar(fileList) - 4)
  
  for (fileName in fileList) {
    
    cat("plotting data for", paste(fileName, ".asc", sep = ""), "\n")
    
    Rst <- raster(paste(inDir, fileName, ".asc", sep = ""))
    
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
    
    if (fileFormat == "jpeg") {
      # jpeg(filename = paste(outDir, fileName, ".jpg", sep = ""), units = "in", width = 10, height = 10, res = 500)
      dev.off()
      image(Rst, col = Colors, breaks = c(min(rstVals) - 1, rstVals), xaxt='n', yaxt='n', xlab = "", ylab = "")
      map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
      box()
      degAxis(1, cex.axis = 1)
      degAxis(2, cex.axis = 1)
      # dev.off()
    }
    
    if (fileFormat == "png") {
      png(filename = paste(outDir, fileName, ".png", sep = ""), 
          width = 2000, height = 3000 * ((ymax(Rst) - ymin(Rst))/(xmax(Rst) - xmin(Rst))), res = 300)
      image(Rst, col = Colors, breaks = c(min(rstVals) - 1, rstVals))
      dev.off()
    }
    
    if (fileFormat == "inR") {
      dev.new(width = 10, height = 10 * ((ymax(Rst) - ymin(Rst))/(xmax(Rst) - xmin(Rst))))
      plot(Rst, col = Colors, breaks = c(min(rstVals) - 1, rstVals), legend = FALSE)
      map("world", type = "b", col = "gray", add = T, fill = T, resolution = 0, border=FALSE)
    }
    
    rm(Rst, rstVals, negativeNb, positiveNb, zeroExists, oneExists, unilimtedExists, Colors)
  }
}