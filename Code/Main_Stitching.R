# import libraries
library(stringr)
library(tidyverse)

# import stitching functions
source("Simulation_Stitch_Functions.R")

# get subfolders of Results/partial
sbfldrs <- list.files(path = "Results/Partial")

# create paths to subfolders
sbfldr_paths <- map(
  .x = sbfldrs,
  .f = function(fldr) paste0("Results/Partial/", fldr)
)

# create names for summary objects
output_names <- map(
  .x = sbfldrs,
  .f = function(fldr) str_split(string = fldr, pattern = "_partial")
)

# summarize partial results
summary_results <- map(
  .x = sbfldr_paths,
  .f = function(p) try(simulation_stitch(path = p))
)

# save in summaries folder
for (i in 1:length(sbfldrs)) {
  saveRDS(
    object = summary_results[[i]],
    file = paste0("Results/Summaries/", output_names[[i]][[1]][1], ".RDS")
  )
}
