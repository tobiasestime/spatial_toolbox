## Various functions used across scripts
## September 2019

## General

pathArgument <- function(path) {
  path <- as.character(path)
  pathLength <- nchar(path)
  path <- ifelse(substr(path, pathLength, pathLength) == "/", substr(path, 1, pathLength - 1), path)
  return(path)
}

outlier <- function(original) {
  lowerBound <- quantile(original, na.rm = TRUE)[[2]]
  upperBound <- quantile(original, na.rm = TRUE)[[4]]
  iqr <- upperBound - lowerBound
  upperBound <- (iqr * 1.5) + upperBound
  lowerBound <- lowerBound - (iqr * 1.5)
  trimmed <- original < upperBound & original > lowerBound
  trimmed[is.na(trimmed)] <- FALSE
  return(trimmed)
}

corr.summary <- function(corr) {
  r <- round(corr$estimate, 2)
  p <- ifelse(corr$p.value < 0.001, "< 0.001", paste("=", round(corr$p.value, 3)))
  cl <- round(corr$conf.int[1], 2)
  cu <- round(corr$conf.int[2], 2)
  pm <- abs(cl - cu) / 2
  n <- corr$parameter + 2
  sprintf("r = %s ± %s, p %s", r, pm, p)
  sprintf("r = %s ± %s, p %s", r, pm, p)
}

lmProb <- function(model) {
  fstat <- summary(model)$fstatistic
  pvalue <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  return(pvalue[["value"]])
}

zScore <- function(r) {
  z <- r - mean(r)
  (z / sd(r))
}

