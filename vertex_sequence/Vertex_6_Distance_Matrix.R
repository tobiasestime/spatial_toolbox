## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

source("../utility/Functions_General.R")

ARGS <- commandArgs(TRUE)

DISTANCES <- pathArgument(ARGS[1])
MATRICES <- pathArgument(ARGS[2])
STANDARD <- as.numeric(ARGS[3])

distanceFiles <- list.files(DISTANCES)
distanceFiles <- distanceFiles[grepl("distances.rds", distanceFiles)]

## Assumes that first two files are lh and rh standard space distance lists
if (STANDARD) {
  distanceFiles <- distanceFiles[1:2]
}

for (distanceFile in distanceFiles) {
  
  if (STANDARD) {
    distanceFile <- strsplit(distanceFile, "_")[[1]]
    distanceFile[1] <- "standard"
    distanceFile <- paste(distanceFile, collapse = "_")
  }

  matrixFile <- sprintf("%s/%s", MATRICES, gsub("_distances", "_matrix", distanceFile))
  if (!file.exists(matrixFile)) {
    saveRDS(NULL, matrixFile)
    print(distanceFile)  
    distanceList <- readRDS(sprintf("%s/%s", DISTANCES, distanceFile))
    nPoints <- length(distanceList$d)
    distanceMatrix <- matrix(0, nrow = 0, ncol = nPoints)
    for (i in 1:nPoints) {
      distanceRow <- numeric(nPoints)
      for (j in distanceList$j[[i]]) {
        distanceIndex <- which(distanceList$j[[i]] == j)
        distanceRow[j] <- distanceList$d[[i]][distanceIndex]
      }
      distanceMatrix <- rbind(distanceMatrix, as.numeric(distanceRow))
    }
    saveRDS(distanceMatrix, matrixFile)
  }
}

