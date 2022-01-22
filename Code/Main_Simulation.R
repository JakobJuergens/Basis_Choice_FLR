##### Set up R session for simulation #####

### load libraries
suppressMessages(library(refund))
suppressMessages(library(MASS))
suppressMessages(library(tidyverse))
suppressMessages(library(fda))
suppressMessages(library(fpca))
suppressMessages(library(caret))

### load data set
data(gasoline)
# octane <- (gasoline$octane)
NIR <- as.matrix(gasoline$NIR)

##### import simulation functions #####
source("Simulation_Functions.R")
### load code from other source files
source("data_generator.R")
source("k_fold_CV_function.R")

##### Perform Simulations #####

# How to save files? --> Please safe them accordingly in their folder and name them with reps, obs, seed!
# conduct 5x500 reps for seed 100-seed 104

# Jakob 5x500
for (i in 100:149) {
  test_bspline_function <- bspline_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_bspline_function,
    file = paste0("Results/Paper/Partial/bspline_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

# Jonghun 5x500
for (i in 100:149) {
  test_fourier_function <- fourier_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, 
    even_basis = FALSE, debug = TRUE
  )
  saveRDS(
    object = test_fourier_function,
    file = paste0("Results/Paper/Partial/fourier_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

for (i in 100:149) {
  test_monomial_function <- monomial_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, 
    even_basis = FALSE, debug = TRUE
  )
  saveRDS(
    object = test_fourier_function,
    file = paste0("Results/Paper/Partial/monomial_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

# Jona 5x500
for (j in 2:4) {
  for (i in 100:149) {
    test_fpcr <- fpcr_function(
      rep = 100, my_data = NULL, n_obs = 200, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr,
      file = paste0("Results/Paper/Partial/pca_bspline_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}

for (j in 2:4) {
  for (i in 100:149) {
    test_fpcr <- fpcr_monomial_function(
      rep = 100, my_data = NULL, n_obs = 200, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr,
      file = paste0("Results/Paper/Partial/pca_monomial_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}

# Jakob 5x500
for (j in 2:4) {
  for (i in 100:149) {
    test_fpcr2 <- fpcr_fourier_function(
      rep = 100, my_data = NULL, nharm = j, n_obs = 200,
      seed = i, even_basis = FALSE, debug = TRUE
    )
    saveRDS(
      object = test_fpcr2,
      file = paste0("Results/Paper/Partial/pca_fourier_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}

################################################################
##### The following simulations use the original NIR data ######
################################################################

# Jonghun 5x500
for (i in 100:149) {
  test_bspline_function_NIR <- bspline_function(
    rep = 100, my_data = NIR, n_obs = 60, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_bspline_function_NIR,
    file = paste0("Results/Paper/Partial/bspline_NIR_partial/NIRrep100_n_obs60_seed", i, ".RDS")
  )
}

# Jona 5x500
for (i in 100:149) {
  test_fourier_function_NIR <- fourier_function(
    rep = 100, my_data = NIR, n_obs = 60, seed = i, debug = TRUE, even_basis = FALSE
  )
  saveRDS(
    object = test_fourier_function_NIR,
    file = paste0("Results/Paper/Partial/fourier_NIR_partial/NIR_fourier_rep100_n_obs60_seed", i, ".RDS")
  )
}

# Jona 5x500 (Jakob)
for (j in 4:4) {
  for (i in 100:149) {
    test_fpcr_NIR <- fpcr_function(
      rep = 100, my_data = NIR, n_obs = 60, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr_NIR,
      file = paste0("Results/Paper/Partial/pca_bspline_nharm", j, "_NIR_partial/NIRrep100_n_obs60_seed", i, ".RDS")
    )
  }
}

# Jonghun 5x500
for (j in 2:4) {
  for (i in 100:149) {
    test_fpcr2_NIR <- fpcr_fourier_function(
      rep = 100, my_data = NIR, n_obs = 60, nharm = j, 
      even_basis = FALSE, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr2_NIR,
      file = paste0("Results/Paper/Partial/pca_fourier_nharm", j, "_NIR_partial/NIRrep100_n_obs60_seed", i, ".RDS")
    )
  }
}
