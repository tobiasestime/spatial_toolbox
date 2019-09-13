## See Vertex_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/vertex_sequence/")

source("../utility/Functions_General.R")
source("Vertex_Definitions_Regions.R")

ARGS <- commandArgs(TRUE)

PATCHES <- pathArgument(ARGS[1])
ORIGINALS <- pathArgument(ARGS[2])
REGIONS  <- pathArgument(ARGS[3])
HEMISPHERES <- as.character(ARGS[4:length(ARGS)])

originalPatchFiles <- list.files(PATCHES, pattern = ".rds")
vertexFiles <- list.files(ORIGINALS)

for (originalPatch in originalPatchFiles) {
  ## Assumes filename with format: TID_*
  tid <- strsplit(originalPatch, "_")[[1]][1]
  rasterData <- readRDS(sprintf("%s/%s", PATCHES, originalPatch))
  for (hemisphere in HEMISPHERES) {
    vertexFile <- sprintf("%s/%s_%s_original.csv", ORIGINALS, tid, hemisphere)
    if (!(vertexFile %in% vertexFiles)) {
      vertices <- rasterData[[hemisphere]]$patch$MegaPatch2$pdv
      vertices <- vertices[vertices$labels %in% REGION_DEFINITIONS[[REGIONS]], ]
      print(vertexFile)
      write.csv(vertices, file = vertexFile, row.names = FALSE)
    }
  }
}

