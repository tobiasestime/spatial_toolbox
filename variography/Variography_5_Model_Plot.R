## See Variography_0_Notes document for details
## September 2019
## tobias.estime@gmail.com

setwd("spatial_toolbox/variography/")

library(ggplot2)
source("../utility/Functions_General.R")
source("../utility/Functions_Spatial.R")

ARGS <- commandArgs(TRUE)
SUMMARY <- pathArgument(ARGS[1])
EMPIRICALS <- pathArgument(ARGS[2])
PLOTS <- pathArgument(ARGS[3])
HEMISPHERES <- ARGS[4:length(ARGS)]

model_summary <- readRDS(SUMMARY)
empiricalFiles <- list.files(EMPIRICALS, pattern = ".rds")

## Color palette can be expanded if required
colorSet <- c("#109876", "#ab1970", "#36ddb3", "#e51d94")
while (length(colorSet) < length(HEMISPHERES)) {
  colorSet <- c(colorSet, "#222222")
}
hemisphereColors <- as.list(colorSet[1:length(HEMISPHERES)])
names(hemisphereColors) <- HEMISPHERES

for (empiricalFile in empiricalFiles) {
  info <- strsplit(empiricalFile, "_")[[1]]
  tid <- info[1]
  variable <- info[2]
  bivariable <- gsub(".rds", "", info[3])
  print(paste(tid, ":", variable, "|", bivariable, sep = " "))
    
  for (hemisphere in HEMISPHERES) {
    esvg <- readRDS(sprintf("%s/%s", EMPIRICALS, empiricalFile))
    esvg <- esvg[esvg$hemisphere == hemisphere, ]
    maxDistance <- max(esvg$distance)
    
    selectModel <- "stable"
    model <- model_summary[model_summary$tid == tid & model_summary$hemisphere == hemisphere & model_summary$variable == variable & model_summary$bivariable == bivariable, ]
    
    parameters <- list(
      r = model$r,
      s = model$s,
      a = model$a,
      n = model$n
    )
    
    modelType <- ifelse(parameters[["r"]] == 0, "Pure Nugget", "Stable Model")
    
    plotFile <- sprintf("%s/%s_%s.png", PLOTS, gsub(".rds", "", empiricalFile), hemisphere)

    ggplot(esvg, aes(x = distance, y = gamma)) +
      geom_point(color = hemisphereColors[[hemisphere]]) +
      ggtitle(sprintf(
        "%s, %s\n%s | %s\n%s\nrange: %s | partial sill: %s | nugget: %s | alpha: %s\nfraction modeled: %s",
        tid, toupper(hemisphere), variable, bivariable,
        modelType, parameters[["r"]], parameters[["s"]], parameters[["n"]], parameters[["a"]], model$d
      )) +
      geom_hline(yintercept = parameters[["s"]] + parameters[["n"]], color = "#aaaaaa", linetype = "dashed") +
      geom_vline(xintercept = parameters[["r"]], color = "#aaaaaa", linetype = "dashed") +
      stat_function(fun = sprintf("model_%s", selectModel), args = parameters, geom = "line", color = "#222222", lwd = 1.25) +
      xlab("Distance (mm)") + ylab("Gamma") +
      xlim(0, maxDistance) +
      theme_minimal() + theme(plot.title = element_text(size = 12), plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "in"))
    ggsave(plotFile, width = 9.6, height = 6.4, units = "in", dpi = 150)
  }
}

