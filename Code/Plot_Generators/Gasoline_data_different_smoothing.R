##### This script plots the true coefficient functions from Reiss and Ogden (2007)
# Clear Workspace
rm(list = ls())

# Load packages
library(tidyverse)
suppressMessages(library(refund))
suppressMessages(library(MASS))
suppressMessages(library(tidyverse))
suppressMessages(library(fda))
suppressMessages(library(fpca))
suppressMessages(library(caret))

##### import simulation functions #####
# this depends on the grid variable above
source("Simulation_Functions.R")
### load code from other source files
source("data_generator.R")
source("k_fold_CV_function.R")

grid = seq(0, 1, length = 401)
#smooth
f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
bumpy
f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -  4*exp(-0.5*(grid-0.45)^2/0.015^2) +  8*exp(-0.5*(grid-0.6)^2/0.02^2) -  exp(-0.5*(grid-0.8)^2/0.03^2)
#two different variances of error
plot(f1*1/401)

generated_curves <- NIR_curve_generator(n = 200, n_harmonics = 30)
my_data <- as.matrix(generated_curves[, -1])

data = t(my_data)

#print(dim(NIR))
smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(4), norder = 4)
smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
xfdlist = list(smooth_basis=smooth_basis_fd)

plot(smooth_basis_fd, xlab = "wavelength")

