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

##### Perform Simulations #####

# Jakob
for (i in 100:149) {
  test_bspline_function <- bspline_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_bspline_function,
    file = paste0("Results/Paper/Partial/bspline_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

# Jakob
for (i in 140:149) {
  test_fourier_function <- fourier_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, debug = TRUE,
    even_basis = FALSE
  )
  saveRDS(
    object = test_fourier_function,
    file = paste0("Results/Paper/Partial/fourier_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

# Jonghun
for (i in 100:149) {
  test_monomial_function <- monomial_function(
    rep = 100, my_data = NULL, n_obs = 200, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_monomial_function,
    file = paste0("Results/Paper/Partial/monomial_sim_partial/rep100_n_obs200_seed", i, ".RDS")
  )
}

# Jona
for (j in 2:4) {
  for (i in 100:149) {
    test_bspline_fpcr <- fpcr_function(
      rep = 100, my_data = NULL, n_obs = 200, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_bspline_fpcr,
      file = paste0("Results/Paper/Partial/pca_bspline_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}

# Jona
for (j in 2:4) {
  for (i in 100:149) {
    test_fourier_fpcr <- fpcr_fourier_function(
      rep = 100, my_data = NULL, nharm = j, n_obs = 200,
      seed = i, even_basis = FALSE, debug = TRUE
    )
    saveRDS(
      object = test_fourier_fpcr,
      file = paste0("Results/Paper/Partial/pca_fourier_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}

# Jonghun
for (j in 2:4) {
  for (i in 100:149) {
    test_monomial_fpcr <- fpcr_monomial_function(
      rep = 100, my_data = NULL, n_obs = 200, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_monomial_fpcr,
      file = paste0("Results/Paper/Partial/pca_monomial_nharm", j, "_sim_partial/rep100_n_obs200_seed", i, ".RDS")
    )
  }
}
