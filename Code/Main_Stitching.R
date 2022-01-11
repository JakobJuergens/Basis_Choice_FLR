# import stitching functions
source("Simulation_Stitch_Functions.R")

# set paths to partial simulation result folders
bspline_partial_folder <- 'Results/bspline_sim_partial/'
# fourier_partial_folder <- 'Results/fourier_sim_partial'
# pca_bspline_partial_folder <- 'Results/pca_bspline_sim_partial'
# pca_fourier_partial_folder <- 'Results/pca_fourier_sim_partial'

# summarize partial results
bspline_results <- simulation_stitch(path = bspline_partial_folder)
# fourier_results <- simulation_stitch(path = fourier_partial_folder)
# pca_bspline_results <- fpcr_simulation_stitch(path = pca_bspline_partial_folder)
# pca_fourier_results <- fpcr_simulation_stitch(path = pca_fourier_partial_folder)

# save in summaries folder
saveRDS(object = bspline_results, file = 'Results/Summaries_sim/bspline_sim.RDS')
# saveRDS(object = fourier_results, file = 'Results/Summaries_sim/fourier_sim.RDS')
# saveRDS(object = pca_bspline_results, file = 'Results/Summaries_sim/pca_bspline_sim.RDS')
# saveRDS(object = pca_fourier_results, file = 'Results/Summaries_sim/pca_fourier_sim.RDS')