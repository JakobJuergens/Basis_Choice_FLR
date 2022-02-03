##### Set up R session for application #####

### load libraries
suppressMessages(library(refund))
suppressMessages(library(MASS))
suppressMessages(library(tidyverse))
suppressMessages(library(fda))
suppressMessages(library(fpca))
suppressMessages(library(caret))

##### import simulation functions #####
source("Application_Functions.R")

### Run Application for basis expansion
bspline_appl <- bspline_appl_function(fold_size = 5, rep = 1000, seed = 100, debug = TRUE)
saveRDS(object = bspline_appl, file = "Results/Paper/Application/bspline_appl.RDS")
write_csv(x = bspline_appl, file = "Results/Paper/Application/bspline_appl.csv")

fourier_appl <- fourier_appl_function(fold_size = 5, rep = 1000, seed = 100, even_basis = FALSE, debug = TRUE)
saveRDS(object = fourier_appl, file = "Results/Paper/Application/fourier_appl.RDS")
write_csv(x = fourier_appl, file = "Results/Paper/Application/fourier_appl.csv")

monomial_appl <- monomial_appl_function(fold_size = 5, rep = 1000, seed = 100, debug = TRUE)
saveRDS(object = monomial_appl, file = "Results/Paper/Application/monomial_appl.RDS")
write_csv(x = monomial_appl, file = "Results/Paper/Application/monomial_appl.csv")

### Run Application for fpcr
for (j in 2:4) {
  bspline_fpcr_appl <- bspline_fpcr_appl_function(
    fold_size = 5, rep = 1000,
    seed = 100, nharm = j, debug = TRUE
  )
  saveRDS(object = bspline_fpcr_appl, file = paste0("Results/Paper/Application/bspline_fpcr_nharm", j, "_appl.RDS"))
  write_csv(x = bspline_fpcr_appl, file = paste0("Results/Paper/Application/bspline_fpcr_nharm", j, "_appl.csv"))
}

for (j in 2:4) {
  fourier_fpcr_appl <- fourier_fpcr_appl_function(
    fold_size = 5, rep = 1000,
    seed = 100, nharm = j,
    even_basis = FALSE, debug = TRUE
  )
  saveRDS(object = fourier_fpcr_appl, file = paste0("Results/Paper/Application/fourier_fpcr_nharm", j, "_appl.RDS"))
  write_csv(x = fourier_fpcr_appl, file = paste0("Results/Paper/Application/fourier_fpcr_nharm", j, "_appl.csv"))
}

for (j in 2:4) {
  monomial_fpcr_appl <- monomial_fpcr_appl_function(
    fold_size = 5, rep = 1000,
    seed = 100, nharm = j, debug = TRUE
  )
  saveRDS(object = monomial_fpcr_appl, file = paste0("Results/Paper/Application/monomial_fpcr_nharm", j, "_appl.RDS"))
  write_csv(x = monomial_fpcr_appl, file = paste0("Results/Paper/Application/monomial_fpcr_nharm", j, "_appl.csv"))
}
