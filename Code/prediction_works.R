rm(list = ls())
library(refund)
library(tidyverse)
library(fda)
set.seed(100)

data(gasoline)
data <- t(as.matrix(gasoline))

octane <- data[1, ]
waves <- data[-1, ]
wavelengths <- seq(from = 900, to = 1700, by = 2)
test <- seq(1, 30, 1)

waves_test <- waves[, test]
octane_test <- octane[test]
waves_train <- waves[, -test]
octane_train <- octane[-test]

# smallbasis <- create.bspline.basis(rangeval = c(min(wavelengths), max(wavelengths)), 
#                                    nbasis = 20, norder = 5)
# harmFdpar <- fdPar(smallbasis)
# 
# gasoline_fd_train <- smooth.basis(argvals = wavelengths,
#                                   y = waves_train, fdParobj = smallbasis)$fd
#                          
# gasoline_fd_test <- smooth.basis(argvals = wavelengths,
#                                  y = waves_test, fdParobj = smallbasis)$fd
# 
# gasoline_pcaObj <- pca.fd(fdobj = gasoline_fd_train, 
#                           nharm = 4, harmfdPar = harmFdpar, centerfns = FALSE)

smallbasis <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                   nbasis = 20, norder = 5)
harmFdpar <- fdPar(smallbasis)

gasoline_fd_train <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                                  y = waves_train, fdParobj = smallbasis)$fd

gasoline_fd_test <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                                 y = waves_test, fdParobj = smallbasis)$fd

gasoline_pcaObj <- pca.fd(fdobj = gasoline_fd_train,
                          nharm = 4, harmfdPar = harmFdpar, centerfns = FALSE)

# plot(gasoline_fd_train)

betabasis1 <- create.constant.basis(c(0, 30))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betafd2    <- create.bspline.basis(rangeval = c(0, 1700), nbasis = 20, 5)
betafdPar2  <- fdPar(betafd2)
betalist <- list(const=betafdPar1, gasoline_fd=betafdPar2)

f_regress <- fRegress(y = octane_train, 
                      xfdlist = list(#const = rep(1, 30), 
                        gasoline_fd = gasoline_fd_train),
                      betalist = list(gasoline_fd = betafdPar2)) #betalist)

prec <- predict(object = f_regress, 
                newdata = list(#const = fd(matrix(rep(1, 30), 1, 30), betabasis1), 
                               gasoline_fd = gasoline_fd_test))

#  plot the data and the fit
plot(f_regress2$betaestlist[[2]])
MSE <- mean((octane_test - prec)^2)
print(MSE)
plot(prec, octane_test, type = "p", pch = "o")
