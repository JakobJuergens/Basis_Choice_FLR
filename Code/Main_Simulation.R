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
octane <- (gasoline$octane)
NIR <- as.matrix(gasoline$NIR)

##### Define Simulation Parameters #####
### set up "global" variables
n_obs <- 60
nharm <- 4
n_var <- 400
grid <- seq(0, 1, length = n_var + 1)

### set up coefficient "functions" / error terms
# smooth
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
# bumpy
f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) - exp(-0.5 * (grid - 0.8)^2 / 0.03^2)
# two different variances of error
sigma_eps_squared1_1 <- as.numeric((var(NIR %*% f1) / 0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 <- as.numeric((var(NIR %*% f1) / 0.6) - var(NIR %*% f1))
sigma_eps_squared2_1 <- as.numeric((var(NIR %*% f2) / 0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 <- as.numeric((var(NIR %*% f2) / 0.6) - var(NIR %*% f2))

##### import simulation functions #####
# this depends on the grid variable above
source("Simulation_Functions.R")
### load code from other source files
source("data_generator.R")
source("k_fold_CV_function.R")

##### Perform Simulations #####

# How to save files? --> Please safe them accordingly in their folder and name them with reps, obs, seed!
# conduct 5x500 reps for seed 100-seed 104

# Jakob 5x500
for (i in 100:104) {
  test_bspline_function <- bspline_function(
    rep = 500, my_data = NULL, n_obs = 200, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_bspline_function,
    file = paste0("Results/Partial/bspline_sim_partial/rep500_n_obs200_seed", i, ".RDS")
  )
}

# Jonghun 5x500
for (i in 100:104) {
  test_fourier_function <- fourier_function(
    rep = 500, my_data = NULL, n_obs = 200, seed = i, 
    even_basis = FALSE, debug = TRUE
  )
  saveRDS(
    object = test_fourier_function,
    file = paste0("Results/Partial/fourier_sim_partial/rep500_n_obs200_seed", i, ".RDS")
  )
}

# Jona 5x500
for (j in 2:4) {
  for (i in 100:104) {
    test_fpcr <- fpcr_function(
      rep = 500, my_data = NULL, n_obs = 200, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr,
      file = paste0("Results/Partial/pca_bspline_nharm", j, "_sim_partial/rep500_n_obs200_seed", i, ".RDS")
    )
  }
}

# Jakob 5x500
for (j in 2:4) {
  for (i in 100:104) {
    test_fpcr2 <- fpcr_fourier_function(
      rep = 500, my_data = NULL, nharm = j, n_obs = 200,
      seed = i, even_basis = FALSE, debug = TRUE
    )
    saveRDS(
      object = test_fpcr2,
      file = paste0("Results/Partial/pca_fourier_nharm", j, "_sim_partial/rep500_n_obs200_seed", i, ".RDS")
    )
  }
}
################################################################
##### The following simulations use the original NIR data ######
################################################################

# Jonghun 5x500
for (i in 100:104) {
  test_bspline_function_NIR <- bspline_function(
    rep = 500, my_data = NIR, n_obs = 60, seed = i, debug = TRUE
  )
  saveRDS(
    object = test_bspline_function_NIR,
    file = paste0("Results/Partial/bspline_NIR_partial/NIRrep500_n_obs60_seed", i, ".RDS")
  )
}

# Jona 5x500
for (i in 100:104) {
  test_fourier_function_NIR <- fourier_function(
    rep = 500, my_data = NIR, n_obs = 60, seed = i, debug = TRUE, even_basis = FALSE
  )
  saveRDS(
    object = test_fourier_function_NIR,
    file = paste0("Results/Partial/fourier_NIR_partial/NIR_fourier_rep500_n_obs60_seed", i, ".RDS")
  )
}

# Jona 5x500 (Jakob)
for (j in 4:4) {
  for (i in 104:104) {
    test_fpcr_NIR <- fpcr_function(
      rep = 500, my_data = NIR, n_obs = 60, nharm = j, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr_NIR,
      file = paste0("Results/Partial/pca_bspline_nharm", j, "_NIR_partial/NIRrep500_n_obs60_seed", i, ".RDS")
    )
  }
}

# Jonghun 5x500
for (j in 2:4) {
  for (i in 100:104) {
    test_fpcr2_NIR <- fpcr_fourier_function(
      rep = 500, my_data = NIR, n_obs = 60, nharm = j, 
      even_basis = FALSE, seed = i, debug = TRUE
    )
    saveRDS(
      object = test_fpcr2_NIR,
      file = paste0("Results/Partial/pca_fourier_nharm", j, "_NIR_partial/NIRrep500_n_obs60_seed", i, ".RDS")
    )
  }
}
