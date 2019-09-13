## Define regions here to provide arguments for various scripts
## Format:
## region_name = numeric
## region_name_labels = character vector
## region_name is the argument you will use for scripts that require a SUB_REGION argument
## Add and remove as necessary
## ! Subject to change as patch labels change

REGION_DEFINITIONS <- list(
  "temporal" = c(1, 2, 3, 4, 7, 13),
  "temporal_labels" = c("Entorhinal", "Inferior Temp.", "Fusiform", "Parahippocampal", "Lingual", "Middle Temp."),
  "standard_temporal_full" = 1:8,
  "standard_temporal_labels" = c("Entorhinal", "Parahippocampal", "Fusiform", "Inferior Temp.", "Middle Temp.", "Lingual")
)