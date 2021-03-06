# import libraries
library(stringr)
library(tidyverse)

# import stitching functions
source("Simulation_Stitch_Functions.R")

# get subfolders of Results/partial
sbfldrs <- list.files(path = "Results/Paper/High_Num_bsplines_fpcr_partial/")

# create paths to subfolders
sbfldr_paths <- map(
  .x = sbfldrs,
  .f = function(fldr) paste0("Results/Paper/High_Num_bsplines_fpcr_partial/", fldr)
)

# create names for summary objects
output_names <- map(
  .x = sbfldrs,
  .f = function(fldr) str_split(string = fldr, pattern = "_partial")[[1]][1]
)

# summarize partial results
summary_results <- map(
  .x = 1:length(sbfldr_paths),
  .f = function(p) {
    try(simulation_stitch(path = sbfldr_paths[[p]]))
  }
)

# save in summaries folder
for (i in 1:length(sbfldrs)) {
  if(length(summary_results[[i]]) != 1){
    saveRDS(
      object = summary_results[[i]],
      file = paste0("Results/Paper/Summaries_additional/", output_names[[i]], ".RDS")
    )
    write_csv(
      x = summary_results[[i]], 
      file = paste0("Results/Paper/Summaries_additional/", output_names[[i]], ".csv")
    )
  }
}
