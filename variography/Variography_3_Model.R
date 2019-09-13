## See Variography_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/variography/")

source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")

ARGS <- commandArgs(TRUE)
EMPIRICALS <- pathArgument(ARGS[1])
MODELS <- pathArgument(ARGS[2])
HEMISPHERES <- ARGS[3:length(ARGS)]

empiricalFiles <- list.files(EMPIRICALS, pattern = ".rds")

fitError <- function(e) { 0 }
aicError <- function(e) { Inf }

fitSummary <- function(summaryStats, defaultNugget) {
  modelStats <- list()
  modelStats[["r"]] <- c(0, 0, 0, 0)
  modelStats[["s"]] <- c(0, 0, 0, 0)
  modelStats[["a"]] <- c(0, 0, 0, 0)
  modelStats[["n"]] <- c(0, 0, 0, 0)
  modelStats[["sigma"]] <- Inf
  modelStats[["df"]] <- c(0, 0)
  
  if (typeof(summaryStats) == "list") {
    parameterNames <- row.names(summaryStats$coefficients)
    for (parameterName in parameterNames) {
      modelStats[[parameterName]] <- as.vector(summaryStats$coefficients[parameterName, ])
    }
    modelStats[["sigma"]] <- summaryStats$sigma
    modelStats[["df"]] <- summaryStats$df
  }
  else {
    modelStats[["n"]] <- c(defaultNugget, 0, 0, 0)
    modelStats[["a"]] <- c(1, 0, 0, 0)
  }
  return(modelStats)
}

## If true, fix the variogram sill to a certain value such as the maximum gamma in the empirical variogram
sillTypes <- c("free", "fixed")

for (empiricalFile in empiricalFiles) {
  print(empiricalFile)

  modelFile <- sprintf("%s/%s", MODELS, gsub(".rds", "_model.rds", empiricalFile))
  
  modelSummary <- list()
  esvg <- readRDS(sprintf("%s/%s", EMPIRICALS, empiricalFile))
  
  for (hemisphere in HEMISPHERES) {
    esvgHemisphere <- esvg[esvg$hemisphere == hemisphere, ]
    rownames(esvgHemisphere) <- NULL
    hemispherePoints <- nrow(esvgHemisphere)
    ## Number of points from variogram to be included in the model
    ## Modeling will use a minimum of 20% to a maximum of 60% of distance covered by the empirical variogram
    distanceFractions <- floor(hemispherePoints * 0.2):floor(hemispherePoints * 0.6)
    
    for (sillType in sillTypes) {
      for (distanceFraction in distanceFractions) {
        esvgSegment <- esvgHemisphere[1:distanceFraction, ]

        ## Parameter bounds and initial parameter estimation
        
        maxRange <- max(esvgHemisphere$distance)
        
        ## Max sill may also be set to variance of original data to fix that outcome in models
        maxSill <- max(esvgHemisphere$gamma)
        minSill <- min(esvgSegment$gamma)
        
        initialGamma <- esvgHemisphere$gamma[1]
        trendGamma <- esvgHemisphere$gamma[5]
        
        ## Positive or negative variograms considered for cross-variograms
        direction <- (trendGamma - initialGamma) >= 0
        maxSill <- ifelse(direction, maxSill, 0)
        minSill <- ifelse(direction, 0, minSill)
              
        ## Fixed sill option
        if (sillType == "fixed") {
          minSill <- maxSill
        }

        maxNugget <- ifelse(direction, initialGamma, initialGamma)
        minNugget <- ifelse(direction, 0, initialGamma)
        
        ## Specify what fraction of the distance covered by the empirical variogram was used for the model (see comment on distanceFraction)
        fractionModeled <- sprintf("d,%s", round(distanceFraction / hemispherePoints, 1))

        ## Upper and lower bounds for parameters
        ## r - range, s - sill, a - parameter specific to stable model, n - nugget
        stableUpper <- c(r = maxRange, s = maxSill, a = 2, n = maxNugget)
        stableLower <- c(r = 0, s = minSill, a = 1, n = minNugget)
        ## Starting parameters
        stableParameters <- list(r = max(esvgSegment$distance), s = (minSill + maxSill) / 2, a = 2, n = (minNugget + maxNugget) / 2)
        
        stableFit <- tryCatch(nls(gamma ~ model_stable(distance, r, s, a, n), data = esvgSegment, start = stableParameters, algorithm = "port", lower = stableLower, upper = stableUpper), error = fitError)
        stableAIC <- tryCatch(AIC(stableFit), error = aicError)
        stableFitSummary <- tryCatch(summary(stableFit), error = fitError)
        modelSummary[[hemisphere]][["stable"]][[sillType]][[fractionModeled]] <- fitSummary(stableFitSummary, mean(esvgHemisphere$gamma))
        modelSummary[[hemisphere]][["stable"]][[sillType]][[fractionModeled]]["AIC"] <- stableAIC
      }
    }
  }
  saveRDS(modelSummary, modelFile)
}

