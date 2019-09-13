## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

library(spdep)
library(rgdal)
source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")

ARGS <- commandArgs(TRUE)

FINALS <- pathArgument(ARGS[1])
POLYGONS <- pathArgument(ARGS[2])
DISTANCES <- pathArgument(ARGS[3])
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

  polygonFile <- sprintf("%s/%s", POLYGONS, gsub("final.csv", "polygon.rds", aggregateName))
  distanceFile <- sprintf("%s/%s", DISTANCES, gsub("final.csv", "distances.rds", aggregateName))
  aggregateFile <- sprintf("%s/%s", FINALS, aggregateFile)
  
  if (!file.exists(distanceFile)) {
    print(aggregateFile)
    saveRDS(NULL, distanceFile)
    aggregateTable <- read.csv(aggregateFile)
    aggregatePolygon <- readRDS(polygonFile)
  
    vertexCount <- nrow(aggregateTable)
    vertexList <- 1:vertexCount
    
    allNeighbors <- poly2nb(aggregatePolygon)
    allDistances <- vector("list", vertexCount)
    
    ## Check for polygons without any adjacent polygons (neighborless)
    removes <- c()
    for (i in vertexList) {
      if (allNeighbors[[i]][1] != 0) {
        ## Three-space distance
        allDistances[[i]] <- spaceDistance(aggregateTable, i, allNeighbors[[i]])
      }
      else {
        removes <- c(removes, i)
      }
    }
    
    ## Remove neighborless polygons
    if (length(removes) != 0) {
      aggregateTable <- aggregateTable[-removes, ]
      aggregatePolygon <- aggregatePolygon[-removes, ]
      
      rownames(aggregateTable) <- NULL
      rownames(aggregatePolygon@data) <- NULL
      
      vertexCount <- nrow(aggregateTable)
      vertexList <- 1:vertexCount
      
      allNeighbors <- poly2nb(aggregatePolygon)
      allDistances <- vector("list", vertexCount)
      
      for (i in vertexList) {
        allDistances[[i]] <- spaceDistance(aggregateTable, i, allNeighbors[[i]])
      }
    }
    
    lags <- vector("list", vertexCount)
    lagDistances <- vector("list", vertexCount)
    
    ## Produce distance list
    for (i in vertexList) {
      print(paste(round((i / vertexCount) * 100, 2), "%", sep = ""))
      neighborTable <- list(neighbors = i, distances = 0, limits = Inf)
      neighborTable <- getDistanceLag(allNeighbors, allDistances, neighborTable)
      lags[[i]] <- neighborTable$neighbors
      lagDistances[[i]] <- neighborTable$distances
    }
    
    write.csv(aggregateTable, aggregateFile, row.names = FALSE)
    saveRDS(aggregatePolygon, polygonFile)
    
    check <- sum(unlist(lapply(lags, length)) == 0)
    
    if (check <= 0) {
      distanceList <- list("j" = lags, "d" = lagDistances)
      saveRDS(distanceList, file = distanceFile)
      saveRDS(allNeighbors, file = gsub("distances.rds", "neighbors.rds", distanceFile))
    }
    else {
      ## This step should not be reached
      ## This should only occur if there is a structural problem relating to neighborless polygons
      print(paste("Error:", print(aggregateFile)))
    }
  }
}

