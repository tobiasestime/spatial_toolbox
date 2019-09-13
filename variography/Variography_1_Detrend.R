## See Variography_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

source("../utility/Functions_General.R")

ARGS <- commandArgs(TRUE)
AGGREGATES <- pathArgument(ARGS[1])
DETREND <- pathArgument(ARGS[2])
VARIABLE <- as.character(ARGS[3])

aggregateFiles <- list.files(AGGREGATES, pattern = ".csv")

## This step is performed because the number of related variables (NC, PVC, multiple time points, etc.) is not assessed in advance
## Rather than constantly changing the size of the regression summary, we assume no more than 100 related variables per vertex table
nVars <- length(aggregateFiles) * 100

regressionSummary <- data.frame(
  tid = rep(NA, nVars),
  hemisphere = rep(NA, nVars),
  variable = rep(NA, nVars),
  rSquared = rep(NA, nVars),
  adjRSquared = rep(NA, nVars),
  rmse = rep(NA, nVars),
  pvalue = rep(NA, nVars),
  intercept = rep(NA, nVars),
  slopeX = rep(NA, nVars),
  slopeY = rep(NA, nVars),
  stringsAsFactors = FALSE
)

## Counter for regression summary row
summaryRow <- 1

for (aggregateFile in aggregateFiles) {
  print(aggregateFile)
  info <- strsplit(aggregateFile, "_")[[1]]
  tid <- info[1]
  hemisphere <- info[2]
  aggregateTable <- read.csv(sprintf("%s/%s", AGGREGATES, aggregateFile))

  columns <- colnames(aggregateTable)
  variableColumns <- columns[grep(VARIABLE, columns)]
  
  for (variableColumn in variableColumns) {
    ## Subsequent scripts will call for detrended variables regardless of whether or not variable was detrended
    ## Summary indicates which variables have actually been detrended
    aggregateTable[, sprintf("%s_d", variableColumn)] <- aggregateTable[, variableColumn]

    regress <- lm(aggregateTable[, variableColumn] ~ aggregateTable$xp + aggregateTable$yp)
    regressSummary <- summary(regress)
    regressCoefficients <- regressSummary$coefficients
    regressR2 <- regressSummary$r.squared
    regressAR2 <- regressSummary$adj.r.squared
    regressSigma <- regressSummary$sigma
    regressP <- lmProb(regress)
    regressInt <- regressCoefficients[1, "Estimate"]
    regressXP <-regressCoefficients[2, "Estimate"]
    regressYP <- regressCoefficients[3, "Estimate"]
    
    regressionSummary[summaryRow, ] <- c(
      tid, hemisphere, variableColumn, regressR2, regressAR2, regressSigma, regressP, regressInt, regressXP, regressYP
    )
    
    if (regressP < 0.05) {
      aggregateTable[, sprintf("%s_d", variableColumn)] <- regress$residuals
    }
    summaryRow <- summaryRow + 1
  }
  ## Overwrite vertex file input
  write.csv(aggregateTable, sprintf("%s/%s", AGGREGATES, aggregateFile), row.names = FALSE)
}

regressionSummary <- regressionSummary[!is.na(regressionSummary$tid), ]
regressionSummary$tid <- as.factor(regressionSummary$tid)
regressionSummary$hemisphere <- as.factor(regressionSummary$hemisphere)
regressionSummary$rSquared <- as.numeric(regressionSummary$rSquared)
regressionSummary$adjRSquared <- as.numeric(regressionSummary$adjRSquared)
regressionSummary$rmse <- as.numeric(regressionSummary$rmse)
regressionSummary$pvalue <- as.numeric(regressionSummary$pvalue)
regressionSummary$intercept <- as.numeric(regressionSummary$intercept)
regressionSummary$slopeX <- as.numeric(regressionSummary$slopeX)
regressionSummary$slopeY <- as.numeric(regressionSummary$slopeY)
## Denote which variables (rows) have been detrended
regressionSummary$detrend <- as.numeric(regressionSummary$pvalue < 0.05)
regressionSummary$detrend <- as.factor(regressionSummary$detrend)

saveRDS(regressionSummary, sprintf("%s/%s", DETREND, "summary_detrend.rds"))

