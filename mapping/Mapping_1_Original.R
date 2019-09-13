## See Mapping_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/mapping/")

library(spdep)
library(rgeos)
library(RColorBrewer)
library(classInt)
source("../utility/Functions_General.R")
source("../vertex_sequence/Vertex_Definitions_Regions.R")

ARGS <- commandArgs(TRUE)
AGGREGATES <- pathArgument(ARGS[1])
POLYGONS <- pathArgument(ARGS[2])
MAPS <- pathArgument(ARGS[3])
REGIONS <- as.character(ARGS[4])
VARIABLE <- as.character(ARGS[5])
QUANTILES <- as.numeric(ARGS[6])
STANDARD <- as.numeric(ARGS[7])

aggregateFiles <- list.files(AGGREGATES, pattern = ".csv")

if (STANDARD) {
  standardPolygons <- list(
    "lh" = readRDS(sprintf("%s/standard_lh_polygon.rds", POLYGONS)),
    "rh" = readRDS(sprintf("%s/standard_rh_polygon.rds", POLYGONS))
  )
}

for (aggregateFile in aggregateFiles) {
  info <- strsplit(aggregateFile, "_")[[1]]
  tid <- info[1]
  hemisphere <- info[2]
  
  aggregateTable <- read.csv(sprintf("%s/%s", AGGREGATES, aggregateFile), stringsAsFactors = FALSE)

  variables <- colnames(aggregateTable)
  variables <- variables[grepl(VARIABLE, variables)]
  
  for (variable in variables) {
    variablePlot <- sprintf("%s/%s_%s_%s.png", MAPS, tid, hemisphere, variable)
    
    if (!file.exists(variablePlot)) {
      print(paste(tid, hemisphere, variable, collapse = " : "))
      if (STANDARD) {
        aggregatePolygon <- standardPolygons[[hemisphere]]
      } else {
        aggregatePolygon <- readRDS(sprintf("%s/%s_%s_polygon.rds", POLYGONS, tid, hemisphere))
      }
      
      aggregatePolygon@data[, variable] <- aggregateTable[, variable]
      
      regions <- gUnaryUnion(aggregatePolygon, id = aggregatePolygon@data$labels)
      spRegions <- list("sp.polygons", SpatialPolygons(regions@polygons), fill = "#eeeeee", col = "#aaaaaa", lwd = 2)
      spText <- list("sp.text", coordinates(regions), toupper(REGION_DEFINITIONS[[sprintf("%s_labels", REGIONS)]]), cex = 0.8, col = "#333333", fontface = "bold")
      
      classBreaks <- classIntervals(aggregatePolygon@data[, variable], QUANTILES, style = "quantile")
      classBreaks <- classBreaks$brks
      classBreaks[1] <- classBreaks[1] - 0.002
      classBreaks[length(classBreaks)] <- classBreaks[length(classBreaks)] + 0.002
      colorPalette <- rev(brewer.pal(n = QUANTILES, name = "RdYlBu"))
      colorPalette <- adjustcolor(colorPalette, alpha.f = 0.8)
      
      plotTitle <- sprintf("%s %s\n%s", tid, hemisphere, variable)
      png(variablePlot, width = 1200, height = 1200, units = "px", res = 150)
      print(spplot(aggregatePolygon, variable, col = "transparent", main = list(label = plotTitle, cex = 1.6),
                   sp.layout = list(spRegions, spText),
                   col.regions = colorPalette,
                   at = classBreaks,
                   xlim = c(aggregatePolygon@bbox["x", "min"] - 10, aggregatePolygon@bbox["x", "max"] + 10), ylim = c(aggregatePolygon@bbox["y", "min"] - 10, aggregatePolygon@bbox["y", "max"] + 10),
                   par.settings = list(axis.line = list(col = "transparent"))))
      dev.off()
    } 
  }
}

