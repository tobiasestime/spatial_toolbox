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

nRows <- length(HEMISPHERES) * length(empiricalFiles)

## r - range
## s - sill
## a - parameter specific to stable model
## n - nugget
## d - fraction of empirical variogram used for modeling
summaryTable <- data.frame(
  "tid" = character(nRows),
  "hemisphere" = character(nRows),
  "variable" = character(nRows),
  "bivariable" = character(nRows),
  "r" = numeric(nRows),
  "s" = numeric(nRows),
  "a" = numeric(nRows),
  "n" = numeric(nRows),
  "d" = numeric(nRows),
  stringsAsFactors = FALSE
)

summaryRow <- 1

for (empiricalFile in empiricalFiles) {
  print(empiricalFile)
  info <- strsplit(empiricalFile, "_")[[1]]
  tid <- info[1]
  variable <- info[2]
  bivariable <- gsub(".rds", "", info[3])
  
  for (hemisphere in HEMISPHERES) {
    modelFile <- sprintf("%s/%s", MODELS, gsub(".rds", "_model.rds", empiricalFile))
    esvg <- readRDS(sprintf("%s/%s", EMPIRICALS, empiricalFile))
    
    esvg <- esvg[esvg$hemisphere == hemisphere, ]
    
    variogramModel <- readRDS(modelFile)
    variogramModel <- variogramModel[[hemisphere]]
    
    sillTypes <- names(variogramModel[[1]])
    distanceFractions <- names(variogramModel[[1]][[1]])
    
    ## Compare half of empirical data to modeled data
    comparisonFraction <- round(nrow(esvg) * 0.5)
    empiricalData <- esvg$gamma[1:comparisonFraction]
    
    ## Start with pure nugget model
    selectDistance <- "d,0"
    selectSill <- "free"
    selectModel <- "nugget"
    r <- 0
    s <- 0
    a <- 1
    n <- mean(esvg$gamma[1:comparisonFraction])  
    
    ## (*) Substitute other criteria for model selection here
    bestCriteria <- sqrt(mean((n - empiricalData) ^ 2))
    
    ## Replace previous selected model and distance modeled with stronger (based on bestCriteria) model if available
    for (sillType in sillTypes) {
        for (distanceFraction in distanceFractions) {
        rr <- variogramModel[["stable"]][[sillType]][[distanceFraction]][["r"]][1]
        ss <- variogramModel[["stable"]][[sillType]][[distanceFraction]][["s"]][1]
        aa <- variogramModel[["stable"]][[sillType]][[distanceFraction]][["a"]][1]
        nn <- variogramModel[["stable"]][[sillType]][[distanceFraction]][["n"]][1]

        modeledData <- model_stable(esvg$distance[1:comparisonFraction], rr, ss, aa, nn)

        ## Match criteria to previous at line (*)
        criteria <- sqrt(mean((modeledData - empiricalData) ^ 2))
        
        if (criteria < bestCriteria) {
          selectModel <- "stable"
          selectDistance <- distanceFraction
          selectSill <- sillType
          bestCriteria <- criteria
        }
      }
    }
    
    ## Select parameters, round, add to table
    if (selectModel != "nugget") {
      r <- variogramModel[[selectModel]][[selectSill]][[selectDistance]][["r"]][1]
      s <- variogramModel[[selectModel]][[selectSill]][[selectDistance]][["s"]][1]
      a <- variogramModel[[selectModel]][[selectSill]][[selectDistance]][["a"]][1]
      n <- variogramModel[[selectModel]][[selectSill]][[selectDistance]][["n"]][1]
    }
    
    r <- round(ifelse(is.null(r), 0, r), 0)
    s <- round(ifelse(is.null(s), 0, s), 6)
    a <- round(ifelse(is.null(a), 0, a), 3)
    n <- round(ifelse(is.null(n), 0, n), 6)
    
    summaryTable[summaryRow, colnames(summaryTable)] <- c(
      "tid" = tid,
      "hemisphere" = hemisphere,
      "variable" = variable,
      "bivariable" = bivariable,
      "r" = r,
      "s" = s,
      "a" = a,
      "n" = n,
      "d" = as.numeric(strsplit(selectDistance, split = ",")[[1]][2])
    )
    summaryRow <- summaryRow + 1
  }
}

summaryTable$hemisphere <- as.factor(summaryTable$hemisphere)
summaryTable$variable <- as.factor(summaryTable$variable)
summaryTable$bivariable <- as.factor(summaryTable$bivariable)
summaryTable$r <- as.numeric(summaryTable$r)
summaryTable$s <- as.numeric(summaryTable$s)
summaryTable$a <- as.numeric(summaryTable$a)
summaryTable$n <- as.numeric(summaryTable$n)
summaryTable$d <- as.numeric(summaryTable$d)
summaryTable <- summaryTable[!duplicated(summaryTable), ]

saveRDS(summaryTable, sprintf("%s/summary_models.rds", MODELS))

