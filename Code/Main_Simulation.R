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

test_bspline_function <- bspline_function(rep = 500, my_data = NULL, n_obs = 200, seed = 10, debug = TRUE)
saveRDS(test_bspline_function,file = "Results/bspline_sim_partial/rep500_n_obs200_seed102.RDS")

test_fourier_function <- fourier_function(rep = 1000, my_data = NULL, n_obs = 200, seed = 100, debug = TRUE)
saveRDS(test_fourier_function,file = "Results/fourier_sim_partial/n_obs200_seed100.RDS")

test_fpcr <- fpcr_function(rep = 1000, my_data = NULL, n_obs = 200, seed = 100, debug = TRUE)
saveRDS(test_fpcr,file="Results/pca_bspline_sim_partial/n_obs200_seed100.RDS")

test_fpcr2 <- fpcr_fourier_function(rep = 1000, my_data = NULL, n_obs = 200, seed = 100, debug = TRUE)
saveRDS(test_fpcr2, file = "Results/pca_fourier_sim_partial/n_obs200_seed100.RDS")

test_bspline_function <- bspline_function(rep = 1000, my_data = NIR, n_obs = 60)
# write.table(test_bspline_function,file="Results/test_bspline_expansion_NIR.csv")

test_fourier_function <- fourier_function(rep = 1000, my_data = NIR, n_obs = 60)
# write.table(test_fourier_function,file="Results/test_fourier_expansion_NIR.csv")

test_fpcr <- fpcr_function(rep = 1000, my_data = NIR, n_obs = 60)
# write.table(test_fpcr,file="Results/test_fpcr_bsplines_NIR.csv")

test_fpcr2 <- fpcr_fourier_function(rep = 1000, my_data = NIR, n_obs = 60)
# write.table(test_fpcr2, file = "Results/test_fpcr_fourier_NIR.csv")
