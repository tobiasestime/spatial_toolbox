## Various functions used across scripts
## September 2019

library(raster)

## Variography

getBins <- function(fromTable, weight, extendTo, precision, perBin) {
  vertexList <- as.numeric(row.names(fromTable))
  maxLags <- extendTo / precision
  distanceBins <- data.frame(from = numeric(maxLags), to = numeric(maxLags), nPairs = numeric(maxLags))
  colNames <- colnames(distanceBins)
  
  ii <- 1
  nPairs <- 0
  startRange <- as.numeric(0)
  saveStart <- 0
  endRange <- precision
  while(endRange <= extendTo) {
    for (vertex in vertexList) {
      allDistances <- weight$d[[vertex]]
      distanceSet <- which(allDistances > startRange & allDistances <= endRange)
      nPairs <- nPairs + length(distanceSet)
    }
    if (nPairs >= perBin) {
      distanceBins[ii, colNames] <- c(saveStart, endRange, nPairs)
      nPairs <- 0
      saveStart <- endRange
      ii <- ii + 1
    }
    startRange <- endRange
    endRange <- endRange + precision
  }
  
  return(c(0, distanceBins[distanceBins$to > 0, "to"]))
}

variance <- function(tableFrom, variable, bivariable, atVertex, toVertices) {
  z1i <- tableFrom[as.character(atVertex), variable]
  z1j <- tableFrom[as.character(toVertices), variable]
  z2i <- tableFrom[as.character(atVertex), bivariable]
  z2j <- tableFrom[as.character(toVertices), bivariable]
  ((z1i - z1j) * (z2i - z2j))
}

empiricalVariogram <- function(fromTable, weight, variable, bivariable, distanceBins) {
  vertexList <- as.numeric(row.names(fromTable))
  nBins <- length(distanceBins)
  empirical <- data.frame(distance = distanceBins, gamma = numeric(nBins), gammaOutliers = numeric(nBins), nPairs = numeric(nBins), outliers = numeric(nBins))
  
  cloud <- vector("list", length(distanceBins))
  
  for (vertex in vertexList) {
    allDistances <- weight$d[[vertex]]
    allGammas <- variance(fromTable, variable, bivariable, vertex, weight$j[[vertex]])
    
    j <- 1
    while (j < nBins) {
      selectGamma <- allGammas[which(allDistances > distanceBins[j] & allDistances <= distanceBins[j + 1])]
      selector <- empirical$distance == distanceBins[j]
      cloud[[j]] <- c(cloud[[j]], selectGamma)
      empirical[selector, "gammaOutliers"] <- empirical[selector, "gammaOutliers"] + sum(selectGamma)
      empirical[selector, "nPairs"] <- empirical[selector, "nPairs"] + length(selectGamma)
      j <- j + 1
    }
  }
  
  cloud[empirical$nPairs <= 0] <- NULL
  empirical <- empirical[empirical$nPairs > 0, ]

  index <- 1
  for (subCloud in cloud) {
    lowerBound <- quantile(subCloud)[[2]]
    upperBound <- quantile(subCloud)[[4]]
    iqr <- upperBound - lowerBound
    upperBound <- (iqr * 1.5) + upperBound
    lowerBound <- lowerBound - (iqr * 1.5)
    trimmedCloud <- subCloud <= upperBound & subCloud >= lowerBound
    empirical$gamma[index] <- mean(subCloud[trimmedCloud]) / 2
    empirical$outliers[index] <- sum(!trimmedCloud)
    index <- index + 1
  }
  
  empirical$distance <- (empirical$distance + distanceBins[2:nBins]) / 2
  empirical$gammaOutliers <- empirical$gammaOutliers / (2 * empirical$nPairs)
  
  return(empirical)
}

semivariogramCloud <- function(fromTable, variable) {
  vertexList <- 1:nrow(fromTable)
  nPairs <- nrow(fromTable) - 1
  svgLength <- (nPairs) * nrow(fromTable)
  cloud <- data.frame(distance = numeric(svgLength), gamma = numeric(svgLength))
  ii <- 1
  for (vertex in vertexList) {
    cloud[ii:(ii - 1 + nPairs), "distance"] <- spaceDistance(fromTable, vertex, vertexList[-vertex])
    cloud[, "gamma"] <- variance(fromTable, variable, vertex)
    ii <- ii + nPairs
  }
  return(cloud)
}

## Variogram models

model_stable <- function(h, r, s, a, n) {
  gamma <- n + (s * (1 - exp(-1 * (((3 * h) ^ a) / (r ^ a)))))
  return(gamma)
}

## Lag functions

lag_stable <- function(h, parameters) {
  r <- parameters[["r"]]
  a <- parameters[["a"]]
  lag <- 1 - model_stable(h, r, 1, a, 0)
  diag(lag) <- 0
  lag <- lag / rowSums(lag)
  return(lag)
}

## Distance

spaceDistance <- function(pointsTable, from, to) {
  sqrt((pointsTable[from, "x"] - pointsTable[to, "x"]) ^ 2 + (pointsTable[from, "y"] - pointsTable[to, "y"]) ^ 2 + (pointsTable[from, "z"] - pointsTable[to, "z"]) ^ 2)
}

getDistances <- function(pointsTable, startPoly) {
  distances <- c()
  for (toPoly in 1:nrow(pointsTable)) {
    distances <- c(distances, spaceDistance(pointsTable, startPoly, toPoly))
  }
  return(distances)
}

## Get neighbors by distance and contiguity

getDistanceLag <- function(neighborList, distanceList, selectList) {
  proceedCheck <- TRUE
  while (proceedCheck) {
    features <- selectList$neighbors
    previousDistances <- selectList$distances
    for (feature in features) {
      neighbors <- neighborList[[feature]]
      if (length(neighbors) > 0) {
        neighborDistances <- distanceList[[feature]]
        
        limit <- max(selectList$limits[selectList$neighbors == feature])[1]
        distance <- min(selectList$distances[selectList$neighbors == feature])[1]
        selector <- neighborDistances <= limit
        
        neighbors <- neighbors[selector]
        neighborDistances <- neighborDistances[selector]
        
        setDistance <- distance + neighborDistances
        setLimits <- limit - neighborDistances
        sortOrder <- order(setDistance, decreasing = FALSE)
        
        selectList$distances <- c(selectList$distances, setDistance[sortOrder])
        selectList$neighbors <- c(selectList$neighbors, neighbors[sortOrder])
        selectList$limits <- c(selectList$limits, setLimits[sortOrder])
      }
    }
    
    keeps <- unique(selectList$neighbors)
    keepLength <- length(keeps)
    keepNeighbors <- numeric(keepLength)
    keepDistances <- keepNeighbors
    keepLimits <- keepNeighbors
    for (keepIndex in 1:keepLength) {
      indexList <- which(selectList$neighbors == keeps[keepIndex])
      distanceSet <- selectList$distances[indexList]
      selection <- which(distanceSet == min(distanceSet))[1]
      selection <- indexList[selection]
      keepNeighbors[keepIndex] <- keeps[keepIndex]
      keepDistances[keepIndex] <- selectList$distances[selection]
      keepLimits[keepIndex] <- selectList$limits[selection]
    }
    
    selectList$neighbors <- keepNeighbors
    selectList$distances <- keepDistances
    selectList$limits <- keepLimits
    
    if (length(previousDistances) == keepLength) {
      if (sum(previousDistances %in% keepDistances) == keepLength) {
        proceedCheck <- FALSE
      }
    }
  }
  
  if (keepLength >= 2) {
    selectList$neighbors <- selectList$neighbors[2:keepLength]
    selectList$distances <- selectList$distances[2:keepLength]
    selectList$limits <- NULL
  }
  else {
    selectList$neighbors <- c()
    selectList$distances <- c()
    selectList$limits <- NULL
  }
  
  return(selectList)
}

## Clustering

## Local I calculation
localI <- function(variable, bivariable, lagMatrix) {
  zi <- zScore(variable)
  zj <- zScore(bivariable)
  wij_zj <- lagMatrix %*% as.matrix(zj)
  localI <- (zi * wij_zj)
  iCategory <- ifelse(zi > 0, ifelse(wij_zj > 0, 1, 2), ifelse(wij_zj > 0, 3, 4))
  return(list("localI" = as.numeric(localI), "zi" = zi, "lag" = zScore(as.numeric(wij_zj)), "iCategory" = iCategory))
}

## Permutation (bootstrap) to assess probability
inferenceLocalI <- function(variable, bivariable, lagMatrix, calculatedLocalI, permute) {
  nSamples <- length(variable)
  pValue <- numeric(nSamples)
  print("permuting...")
  for (i in 1:permute) {
    bivariable <- sample(bivariable, nSamples, replace = FALSE)
    permutedI <- localI(variable, bivariable, lagMatrix)
    pValue <- pValue + (permutedI$localI >= calculatedLocalI)
  }
  return((pValue + 1) / (permute + 1))
}

## Raster

fillRaster <- function(x) {
  len <- length(x)
  i <- (len / 2) + 0.5
  if (is.na(x)[i]) {
    if (sum(is.na(x)) != len) {
      return(1)
    }
  }
  return(x[i])
}

