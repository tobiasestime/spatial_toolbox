## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

source("../utility/Functions_General.R")

library(compiler)
library(data.table)
library(dplyr)

ARGS <- commandArgs(TRUE)

ORIGINALS <- pathArgument(ARGS[1])
AGGREGATES <- pathArgument(ARGS[2])
RESOLUTION <- as.numeric(ARGS[3])
STANDARD <- as.numeric(ARGS[4])

originalFiles <- list.files(ORIGINALS, pattern = ".csv")

## Assumes that first two files are lh and rh standard space vertex tables
if (STANDARD) {
  originalFiles <- originalFiles[1:2]
}

## Aggregation function (is run in loop below)
## The table provided must have x, y, z coordinates, and ivx - a unique identifier for each vertex
aggregateVertices <- function(newTable) {
  distanceBreak <- 1

  while (distanceBreak != 0) {
    pointList <- 1:nrow(newTable)
    pointLength <- length(pointList)
    
    to <- numeric((((pointLength + 1) * pointLength) / 2) - pointLength)
    from <- to
    
    j <- 1
    k <- pointLength - 1
    for (i in pointList) {
      if (k > 0) {
        to[j:(j + k - 1)] <- rep(i, k)
        pointList <- pointList[pointList != i]
        from[j:(j + k - 1)] <- pointList
        j <- j + k
        k <- k - 1
      }
    }
    
    pairs <- data.table("from" = from, "to" = to)

    pairs[, distance := sqrt((newTable[from, x] - newTable[to, x]) ^ 2 + (newTable[from, y] - newTable[to, y]) ^ 2 + (newTable[from, z] - newTable[to, z]) ^ 2)]
    pairs[, x := (newTable[to, x] + newTable[from, x]) / 2]
    pairs[, y := (newTable[to, y] + newTable[from, y]) / 2]
    pairs[, z := (newTable[to, z] + newTable[from, z]) / 2]

    pairs <- setorder(pairs, distance)
    pointList <- unique(pairs$from)
 
    mins <- pairs[ , .SD[which.min(distance)], by = from]
    mins <- mins[ , .SD[which.min(distance)], by = to]
    mins <- subset(mins, mins$distance < RESOLUTION)
        
    for (i in mins$from) {
      mins <- mins[mins$to != i, ]
    }
    
    pairs <- pairs[!(pairs$from %in% mins$from), ]
    pairs <- pairs[!(pairs$to %in% mins$to), ]
    pairs <- pairs[!(pairs$from %in% mins$to), ]
    pairs <- pairs[!(pairs$to %in% mins$from), ]
    
    while (sum(pairs$distance < RESOLUTION) != 0) {
      newMins <- pairs[ , .SD[which.min(distance)], by = from]
      newMins <- newMins[ , .SD[which.min(distance)], by = to]
      newMins <- subset(newMins, newMins$distance < RESOLUTION)
      
      for (i in newMins$from) {
        newMins <- newMins[newMins$to != i, ]
      }
      
      mins <- rbind(mins, newMins)
      
      pairs <- pairs[!(pairs$from %in% newMins$from), ]
      pairs <- pairs[!(pairs$to %in% newMins$to), ]
      pairs <- pairs[!(pairs$from %in% newMins$to), ]
      pairs <- pairs[!(pairs$to %in% newMins$from), ]
    }

    distanceBreak <- nrow(mins)
    
    if (distanceBreak > 0) {
      pointList <- mins$from
      for (point in pointList) {
        subpairs <- subset(mins, from == point)
        toPoint <- subpairs$to
        newTable[point, c("x", "y", "z")] <- subpairs[, c("x", "y", "z")]
        ## Track the original vertices that comprise each aggregated point
        newTable[point, "ivx"] <- paste(newTable[point, "ivx"], newTable[toPoint, "ivx"], sep = "|")
        newTable[toPoint, "ivx"] <- NA
      }
    }
    newTable <- subset(newTable, !is.na(newTable$ivx))
  }
  return(newTable)
}

aggregateVertices <- cmpfun(aggregateVertices)

for (originalFile in originalFiles) {

  originalName <- originalFile
  if (STANDARD) {
    originalName <- strsplit(originalFile, "_")[[1]]
    originalName[1] <- "standard"
    originalName <- paste(originalName, collapse = "_")
  }
  
  aggregateFile <- sprintf("%s/%s", AGGREGATES, gsub("original", "aggregate", originalName))
  
  if (!file.exists(aggregateFile)) {
    write.csv(NULL, file = aggregateFile)
    print(paste(which(originalFile == originalFiles), ":", originalFile, ":", Sys.time()))
    
    originalTable <- read.csv(file = sprintf("%s/%s", ORIGINALS, originalFile))
    originalTable$ivx <- as.character(originalTable$ivx)
    originalTable <- data.table(originalTable)
    
    regions <- sort(unique(originalTable$labels))
    
    ## Aggregate one region
    aggregateTable <- aggregateVertices(originalTable[originalTable$labels == regions[1], ])
    
    regions <- regions[2:length(regions)]
    
    ## Aggregate the other regions individually; combine to existing aggregated regions
    for (region in regions) {
      regionAggregate <- aggregateVertices(originalTable[originalTable$labels == region, ])
      aggregateTable <- rbind(aggregateTable, regionAggregate)
    }
    
    ## Finally, aggregate the combined aggregated regions
    aggregateTable <- aggregateVertices(aggregateTable)
    
    write.csv(aggregateTable, file = aggregateFile, row.names = FALSE)
  }
}

