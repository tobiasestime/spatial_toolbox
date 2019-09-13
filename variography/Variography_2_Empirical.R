## See Variography_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/variography/")

library(spdep)
library(ggplot2)
source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")

ARGS <- commandArgs(TRUE)
AGGREGATES <- pathArgument(ARGS[1])
DISTANCES <- pathArgument(ARGS[2])
BINS <- pathArgument(ARGS[3])
EMPIRICALS <- pathArgument(ARGS[4])
VARIABLE <- as.character(ARGS[5])
BIVARIABLE <- as.character(ARGS[6])
STANDARD <- as.numeric(ARGS[7])
AUTO <- as.numeric(ARGS[8])
AUTO <- ifelse(VARIABLE == BIVARIABLE, AUTO, 0)
POINTS <- as.numeric(ARGS[9])
PRECISION <- 0.1

aggregateFiles <- list.files(AGGREGATES, pattern = ".csv")
info <- strsplit(aggregateFiles, "_")
tids <- unique(sapply(info, "[", 1))

for (tid in tids) {
  print(tid)
  
  leftHemisphere <- read.csv(sprintf("%s/%s_lh_final.csv", AGGREGATES, tid))
  rightHemisphere <- read.csv(sprintf("%s/%s_rh_final.csv", AGGREGATES, tid))
  
  structureType <- ifelse(STANDARD, "standard", tid)
  
  leftDistances <- readRDS(sprintf("%s/%s_lh_distances.rds", DISTANCES, structureType))
  rightDistances <- readRDS(sprintf("%s/%s_rh_distances.rds", DISTANCES, structureType))
  
  binSize <- floor((max(nrow(leftHemisphere), nrow(rightHemisphere)) ^ 2) / POINTS)
  
  binFile <- sprintf("%s/%s_%s_bins.rds", BINS, structureType, POINTS)
  
  if (!file.exists(binFile)) {
    saveRDS(NULL, binFile)
    print("binning...")
    
    maxExtent <- floor(max(c(unlist(leftDistances$d), unlist(rightDistances$d))))
    
    leftBins <- getBins(leftHemisphere, leftDistances, maxExtent, PRECISION, binSize)
    rightBins <- getBins(rightHemisphere, rightDistances, maxExtent, PRECISION, binSize)
    
    bins <- list("lh" = leftBins, "rh" = rightBins)
    saveRDS(bins, binFile)
  } else {
    bins <- readRDS(binFile)
  }
  
  leftBins <- bins$lh
  rightBins <- bins$rh
  
  columns <- colnames(leftHemisphere)
  columns <- columns[!grepl("_d", columns)]
  variablePoints <- columns[grepl(VARIABLE, columns)]
  bivariablePoints <- columns[grepl(BIVARIABLE, columns)]
  
  if (length(bivariablePoints) > 0) {
    for (variablePoint in variablePoints) {
      
      if (AUTO) {
        bivariablePoints <- variablePoint
      }
      
      for (bivariablePoint in bivariablePoints) {
        
        variogramFile <- sprintf("%s/%s_%s_%s.rds", EMPIRICALS, tid, variablePoint, bivariablePoint)
        redundantFile <- sprintf("%s/%s_%s_%s.rds", EMPIRICALS, tid, bivariablePoint, variablePoint)
        plotFile <- sprintf("%s/plots/%s_%s_%s.png", EMPIRICALS, tid, variablePoint, bivariablePoint)
        
        if (!file.exists(variogramFile) & !file.exists(redundantFile)) {
          saveRDS(NULL, variogramFile)
          
          print(paste(tid, ":", variablePoint, "|", bivariablePoint, ":", Sys.time()))
          
          print("variographing...")
          
          esvgLH <- empiricalVariogram(leftHemisphere, leftDistances, variablePoint, bivariablePoint, leftBins)
          esvgLH$hemisphere <- c("lh")
          rownames(esvgLH) <- NULL
          
          esvgRH <- empiricalVariogram(rightHemisphere, rightDistances, variablePoint, bivariablePoint, rightBins)
          esvgRH$hemisphere <- c("rh")
          rownames(esvgRH) <- NULL
          
          esvgLHD <- empiricalVariogram(leftHemisphere, leftDistances, sprintf("%s_d", variablePoint), sprintf("%s_d", bivariablePoint), leftBins)
          esvgLHD$hemisphere <- c("lh_d")
          rownames(esvgLHD) <- NULL
          
          esvgRHD <- empiricalVariogram(rightHemisphere, rightDistances, sprintf("%s_d", variablePoint), sprintf("%s_d", bivariablePoint), rightBins)
          esvgRHD$hemisphere <- c("rh_d")
          rownames(esvgRHD) <- NULL
          
          esvg <- rbind(esvgLH, esvgLHD, esvgRH, esvgRHD)
          rownames(esvg) <- NULL
          
          print("saving...")
          
          saveRDS(esvg, variogramFile)
          
          ggplot(esvg, aes(x = distance, y = gamma, group = hemisphere, color = hemisphere)) +
            geom_point() + scale_color_manual(values = c("#109876", "#36ddb3", "#ab1970", "#e51d94")) +
            ggtitle(sprintf("%s\n%s | %s", tid, variablePoint, bivariablePoint)) +
            xlab("Distance (mm)") + ylab("Semivariance") +
            theme(plot.title = element_text(size = 14, hjust = 0.5))
          ggsave(plotFile, width = 9.6, height = 6.4, units = "in", dpi = 150)
        }
      }
    }
  }
}

