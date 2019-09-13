## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

library(concaveman)
library(dismo)
library(raster)
library(rgeos)
library(sf)
library(sp)
source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")
source("Vertex_Definitions_Regions.R")

ARGS <- commandArgs(TRUE)

FINALS <- pathArgument(ARGS[1])
POLYGONS <- pathArgument(ARGS[2])
REGIONS <- as.character(ARGS[3])
STANDARD <- as.numeric(ARGS[4])

aggregateFiles <- list.files(FINALS, pattern = ".csv")

## Assumes that first two files are lh and rh standard space vertex tables
if (STANDARD) {
  aggregateFiles <- aggregateFiles[1:2]
}

for (aggregateFile in aggregateFiles) {
  
  aggregateName <- aggregateFile
  if (STANDARD) {
    aggregateName <- strsplit(aggregateName, "_")[[1]]
    aggregateName[1] <- "standard"
    aggregateName <- paste(aggregateName, collapse = "_")
  }
  
  aggregateTable <- read.csv(sprintf("%s/%s", FINALS, aggregateFile))
  aggregateTable[, c("x3", "y3", "z3")] <- aggregateTable[, c("x", "y", "z")]
  aggregateTable[, c("x", "y", "z")] <- NULL
  
  aggregatePoints <- aggregateTable
  coordinates(aggregatePoints) <- ~xp+yp
  
  concaveHull <- concaveman(aggregatePoints)
  concaveHull <- gBuffer(concaveHull, width = 1, joinStyle = "BEVEL")
  
  polygonFile <- sprintf("%s/%s", POLYGONS, gsub("final.csv", "polygon.rds", aggregateName))
  
  hemisphere <- strsplit(aggregateFile, "_")[[1]][2]
  
  if (!file.exists(polygonFile)) {
    saveRDS(NULL, polygonFile)
    print(polygonFile)
    
    polygonData <- aggregateTable[, c("xp", "yp")]
    colnames(polygonData) <- c("x", "y")
    polygonData <- voronoi(polygonData)
    
    polygonData <- gIntersection(polygonData, concaveHull, byid = TRUE, drop_lower_td = FALSE) 
    
    for (i in 1:length(slot(polygonData, "polygons"))) {
      slot(slot(polygonData, "polygons")[[i]], "ID") = row.names(aggregateTable)[i]
    }
    
    polygonData <- SpatialPolygonsDataFrame(polygonData, aggregateTable)
    
    if (nrow(aggregateTable) != length(polygonData)) {
      print(paste("correcting:", hemisphere))
      aggregateTable[, c("x", "y", "z")] <- aggregateTable[, c("x3", "y3", "z3")]
      aggregateTable[, c("x3", "y3", "z3")] <- NULL
      removes <- which(!(rownames(aggregateTable) %in% polygonData@data$id))
      aggregateTable <- aggregateTable[-removes, ]
      rownames(aggregateTable) <- NULL
      rownames(polygonData@data) <- NULL
      polygonData@data$id <- as.numeric(rownames(aggregateTable))
      write.csv(aggregateTable, sprintf("%s/%s", FINALS, aggregateFile), row.names = FALSE)
    }
    
    saveRDS(polygonData, polygonFile)
  }
}

