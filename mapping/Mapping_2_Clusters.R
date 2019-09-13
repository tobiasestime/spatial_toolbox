## See Mapping_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/mapping/")

library(spdep)
library(rgeos)
library(RColorBrewer)
library(classInt)
source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")
source("../vertex_sequence/Vertex_Definitions_Regions.R")

ARGS <- commandArgs(TRUE)

AGGREGATES <- pathArgument(ARGS[1])
POLYGONS <- pathArgument(ARGS[2])
MODELS <- pathArgument(ARGS[3])
MATRICES <- pathArgument(ARGS[4])
MAPS <- pathArgument(ARGS[5])
REGIONS <- as.character(ARGS[6])
VARIABLE <- as.character(ARGS[7])
SIGNIFICANCE <- as.numeric(ARGS[8])
STANDARD <- as.numeric(ARGS[9])

aggregateFiles <- list.files(AGGREGATES, pattern = ".csv")
modelSummary <- readRDS(MODELS)

if (STANDARD) {
  standardPolygons <- list(
    "lh" = readRDS(sprintf("%s/standard_lh_polygon.rds", POLYGONS)),
    "rh" = readRDS(sprintf("%s/standard_rh_polygon.rds", POLYGONS))
  )
  standardMatrices <- list(
    "lh" = readRDS(sprintf("%s/standard_lh_matrix.rds", MATRICES)),
    "rh" = readRDS(sprintf("%s/standard_rh_matrix.rds", MATRICES))
  )
}

for (aggregateFile in aggregateFiles) {
  info <- strsplit(aggregateFile, "_")[[1]]
  tid <- info[1]
  hemisphere <- info[2]
  
  aggregateTable <- read.csv(sprintf("%s/%s", AGGREGATES, aggregateFile), stringsAsFactors = FALSE)
  
  variables <- colnames(aggregateTable)
  variables <- variables[!grepl("_d", variables)]
  variables <- variables[grepl(VARIABLE, variables)]
  
  for (variable in variables) {
    
    clusterPlot <- sprintf("%s/%s_%s_%s.png", MAPS, tid, hemisphere, variable)
    
    if (!file.exists(clusterPlot)) {
      print(paste(tid, hemisphere, variable, collapse = " : "))
      saveRDS(NULL, clusterPlot)
      
      if (STANDARD) {
        aggregatePolygon <- standardPolygons[[hemisphere]]
        aggregateMatrix <- standardMatrices[[hemisphere]]
      } else {
        aggregatePolygon <- readRDS(sprintf("%s/%s_%s_polygon.rds", POLYGONS, tid, hemisphere))
        aggregateMatrix <- readRDS(sprintf("%s/%s_%s_matrix.rds", MATRICES, tid, hemisphere))
      }
      
      aggregatePolygon@data[, variable] <- aggregateTable[, variable]
      
      ## Parameters for geostatistical weighting function
      selectModel <- modelSummary[modelSummary$tid == tid & modelSummary$hemisphere == paste(hemisphere, "d", sep = "_") & modelSummary$variable == variable & modelSummary$bivariable == variable, ]

      if (selectModel$r != 0) {
        lagParameters <- list(r = selectModel$r, s = selectModel$s, a = selectModel$a, n = selectModel$n)
        lagMatrix <- lag_stable(aggregateMatrix, lagParameters)
                
        iResult <- localI(aggregateTable[, variable], aggregateTable[, variable], lagMatrix)
                
        aggregateTable$localI <- iResult$localI
        aggregateTable$iCategory <- iResult$iCategory
        
        ## Correction for multiple comparisons
        pValue <- 1 - ((1 - SIGNIFICANCE) ^ (1 / nrow(aggregateTable)))
        nPermutations <- round((pValue ^ -1), 0)
        
        aggregateTable$iProbability <- inferenceLocalI(aggregateTable[, variable], aggregateTable[, variable], lagMatrix, aggregateTable$localI, nPermutations)
      
        ## One-sided test
        aggregateTable$iCategory[!(aggregateTable$iProbability >= (1 - pValue) | aggregateTable$iProbability <= pValue)] <- 0
        aggregatePolygon@data$iCategory <- aggregateTable$iCategory
  
        ## Output map
        plotTitle <- sprintf("%s %s\n%s\nGlobal I = %s, p = %s", tid, hemisphere, variable, round(mean(aggregateTable$localI), 2), SIGNIFICANCE)
        
        png(clusterPlot, width = 960, height = 960, units = "px", res = 150)
        print(spplot(aggregatePolygon, "iCategory", col = "transparent", main = list(label = plotTitle, cex = 1.5),
                     col.regions = c("#dddddd", "#dd2626", "#f6a4a4", "#9dbdfd", "#2660dd"),
                     at = seq(-0.5, 4.5, by = 1),
                     colorkey = list(labels = list(labels = c("Not Significant", "High Cluster", "High Outlier", "Low Outlier", "Low Cluster"))),
        par.settings = list(axis.line = list(col = "transparent"))))
        dev.off()
      }
    }
  }
}

