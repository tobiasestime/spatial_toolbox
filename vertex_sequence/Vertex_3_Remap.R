## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")

ARGS <- commandArgs(TRUE)

AGGREGATES <- pathArgument(ARGS[1])
ORIGINALS <- pathArgument(ARGS[2])
FINALS <- pathArgument(ARGS[3])
STANDARD <- as.numeric(ARGS[4])

originalFiles <- list.files(ORIGINALS, pattern = ".csv")

for (originalFile in originalFiles) {
  finalFile <- sprintf("%s/%s", FINALS, gsub("original", "final", originalFile))
  if (!file.exists(finalFile)) {
    write.csv(NULL, finalFile)
    print(originalFile)
    
    aggregateFile <- gsub("original", "aggregate", originalFile)
    if (STANDARD) {
      aggregateFile <- strsplit(aggregateFile, "_")[[1]]
      aggregateFile[1] <- "standard"
      aggregateFile <- paste(aggregateFile, collapse = "_")
    }
    
    aggregateTable <- read.csv(sprintf("%s/%s", AGGREGATES, aggregateFile), stringsAsFactors = FALSE)
    originalTable <- read.csv(sprintf("%s/%s", ORIGINALS, originalFile), stringsAsFactors = FALSE)
    
    excludeColumns <- c("ipvx", "ivx", "ibnd", "xp", "yp", "x", "y", "z", "labels")
    allColumns <- colnames(aggregateTable)
    aggregateColumns <- allColumns[which(!(allColumns %in% excludeColumns))]
    
    for (vertex in 1:nrow(aggregateTable)) {
      ## Obtain original vertices that comprise an aggregated vertex
      vertices <- as.numeric(strsplit(aggregateTable[vertex, "ivx"], split = "\\|")[[1]])
      originalVertices <- originalTable[which(originalTable$ivx %in% vertices), ]
  
      distances <- sqrt((originalVertices$x - aggregateTable[vertex, "x"]) ^ 2 + (originalVertices$y - aggregateTable[vertex, "y"]) ^ 2 + (originalVertices$z - aggregateTable[vertex, "z"]) ^ 2)
      minIndex <- which(distances == min(distances))[1]
  
      ## Re-assign location to the closest original location
      aggregateTable[vertex, excludeColumns] <- originalVertices[minIndex, excludeColumns]
      ## Average the original vertices to obtain a value for the aggregate
      aggregateTable[vertex, aggregateColumns] <- colMeans(originalVertices[, aggregateColumns])
    }
    
    ## Check that should be passed; if not, there may be something strange about the vertex table
    if (sum(duplicated(aggregateTable$ivx)) != 0) {
      print(paste("review:", aggregateFile))
    }
    
    write.csv(aggregateTable, finalFile, row.names = FALSE)
  }
}

