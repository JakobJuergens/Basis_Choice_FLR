rm(list = ls())
library(refund)
library(tidyverse)
library(reshape2)
library(fda)
library(fdaACF)
library(fdapace)
library(fpca)

set.seed(100)

# set up "global" variables
n_obs <- 60
n_var <- 150 # number of vars
test <- seq(1, 20, 1)
grid <- seq(0, 1, length = n_var + 1)

# set up functions here!
# smooth
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
# bumpy
f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) -
  4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) +
  8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) -
  exp(-0.5 * (grid - 0.8)^2 / 0.03^2)

beta1 <- sin(grid * 2 * pi)
beta2 <- -dnorm(grid, mean = .2, sd = .03) + 3 * dnorm(grid, mean = .5, sd = .04) + dnorm(grid, mean = .75, sd = .05)
fun_1 <- 2 * sin(grid * 0.5 * pi) + 4 * sin(grid * 2.5 * pi) + 5 * sin(grid * 2.5 * pi)

nharm <- 4 # number of fpc

data_generation <- function(fun) {
  var1 <- 1
  var2 <- 2
  var3 <- 0.2

  # X <- matrix(0, nrow=n_obs, ncol=length(grid))
  # for(i2 in 1:n_obs){
  #     X[i2,]=X[i2,]+rnorm(length(grid), 0, var1)
  #     X[i2,]=X[i2,]+runif(1, 0, var2)
  #     X[i2,]=X[i2,]+rnorm(1, 0, var3)*grid

  #     for(j2 in 1:10){
  #         e =rnorm(2, 0, var1/j2^(2))
  #         X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
  #         X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
  #     }
  # }
  # return(X)
  # }
  # This is a "more clear" setup

  X <- matrix(0, nrow = n_obs, ncol = length(grid))
  for (i2 in 1:n_obs) {
    X[i2, ] <- X[i2, ] + rnorm(length(grid), 2, var1)
    X[i2, ] <- X[i2, ] + runif(1, 0, var2)
    X[i2, ] <- X[i2, ] + rnorm(1, 0, var3) * grid
    X[i2, ] <- X[i2, ] * fun

    for (j2 in 1:5) {
      e <- abs(rnorm(2, 0, var1 / j2^(2)))
      X[i2, ] <- X[i2, ] + e[1] * sin((2 * pi) * grid * j2)
      X[i2, ] <- X[i2, ] + e[2] * cos((2 * pi) * grid * j2)
    }
  }
  return(X)
}

X <- data_generation(beta2)
Y <- X %*% (beta1 / n_var) + rnorm(n_obs, 0, 0.125)
Y <- as.numeric(Y)
final_X <- t(X)
X_train <- final_X[, -test]
Y_train <- Y[-test]
X_test <- final_X[, test]
Y_test <- Y[test]


smallbasis <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = 12, 7)
# harmFdpar = fdPar(smallbasis)

smooth_basis_train <- smooth.basis(y = X_train, fdParobj = smallbasis)
smooth_basis_train_fd <- smooth_basis_train$fd
smooth_basis_test <- smooth.basis(y = X_test, fdParobj = smallbasis)
smooth_basis_test_fd <- smooth_basis_test$fd


# todo:use to choose appropriate level of explained variance

plot(smooth_basis_train_fd, ylab = "", xlab = "", col = "gray")


xfdlist <- list(smooth_basis = smooth_basis_train_fd)
betabasis1 <- create.constant.basis(c(0, (n_obs - length(test))))
# betabasis_test <- create.constant.basis(c(0, length(test)))
betafd1 <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betafd2 <- create.bspline.basis(rangeval = c(0, n_var + 1), nbasis = 12, 7)
betafdPar2 <- fdPar(betafd2)
betalist <- list(smooth_basis_fd = betafdPar2)

f_regress <- fRegress(y = Y_train, xfdlist, betalist, y2cMap = smooth_basis_train$y2cMap)
# also possible to use formula interface:
# f_regress2 <- fRegress(Y_train ~ smooth_basis_train_fd)
prec <- predict(
  object = f_regress,
  newdata = list(data = smooth_basis_test_fd)
)

# #  plot the data and the fit
# plot(f_regress$betaestlist[[2]])
# add true function plot
MSE_test <- mean((Y_test - prec)^2)
print(MSE_test)
plot(prec, Y_test, type = "p", pch = "o")
