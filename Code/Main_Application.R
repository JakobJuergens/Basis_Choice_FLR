##### Set up R session for application #####

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
source("Application_Functions.R")

### Run Application for basis expansion
bspline_appl <- bspline_appl_function(fold_size = 5, rep = 1, seed = 100, debug = TRUE)
saveRDS(object = bspline_appl, file = 'Results/Paper/Application/bspline_appl.RDS')
write_csv(x = bspline_appl, file = 'Results/Paper/Application/bspline_appl.csv')

fourier_appl <- fourier_appl_function(fold_size = 5, rep = 250, seed = 100, even_basis = FALSE, debug = TRUE)
saveRDS(object = fourier_appl, file = 'Results/Paper/Application/fourier_appl.RDS')
write_csv(x = fourier_appl, file = 'Results/Paper/Application/fourier_appl.csv')

monomial_appl <- monomial_appl_function(fold_size = 5, rep = 250, seed = 100, debug = TRUE)
saveRDS(object = monomial_appl, file = 'Results/Paper/Application/monomial_appl.RDS')
write_csv(x = monomial_appl, file = 'Results/Paper/Application/monomial_appl.csv')

### Run Application for fpcr
bspline_fpcr_appl <- bspline_fpcr_appl_function(fold_size = 5, rep = 250, seed = 100, debug = TRUE)
saveRDS(object = bspline_fpcr_appl, file = 'Results/Paper/Application/bspline_fpcr_appl.RDS')
write_csv(x = bspline_fpcr_appl, file = 'Results/Paper/Application/bspline_fpcr_appl.csv')

fourier_fpcr_appl <- fourier_fpcr_appl_function(fold_size = 5, rep = 250, seed = 100, even_basis = FALSE, debug = TRUE)
saveRDS(object = fourier_fpcr_appl, file = 'Results/Paper/Application/fourier_fpcr_appl.RDS')
write_csv(x = fourier_fpcr_appl, file = 'Results/Paper/Application/fourier_fpcr_appl.csv')

monomial_fpcr_appl <- monomial_fpcr_appl_function(fold_size = 5, rep = 250, seed = 100, debug = TRUE)
saveRDS(object = monomial_fpcr_appl, file = 'Results/Paper/Application/monomial_fpcr_appl.RDS')
write_csv(x = monomial_fpcr_appl, file = 'Results/Paper/Application/monomial_fpcr_appl.csv')

