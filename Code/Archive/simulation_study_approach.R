rm(list = ls())
library(refund)
library(fda)
library(tidyverse)

set.seed(100)

n_obs <- 60
n_var <- 150 # number of vars
test <- seq(1, 20, 1)
grid <- seq(0, 1, length = n_var + 1)
beta1 <- sin(grid * 2 * pi)
beta2 <- -dnorm(grid, mean = .2, sd = .03) + 3 * dnorm(grid,
  mean = .5,
  sd = .04
) + dnorm(grid, mean = .75, sd = .05)

nharm <- 4 # number of fpc

data_generation <- function() {
  var1 <- 1
  var2 <- 5
  var3 <- 0.2

  X <- matrix(0, nrow = n_obs, ncol = length(grid))
  for (i2 in 1:n_obs) {
    X[i2, ] <- X[i2, ] + rnorm(length(grid), 0, var1)
    X[i2, ] <- X[i2, ] + runif(1, 0, var2)
    X[i2, ] <- X[i2, ] + rnorm(1, 0, var3) * grid

    for (j2 in 1:10) {
      e <- rnorm(2, 0, var1 / j2^(2))
      X[i2, ] <- X[i2, ] + e[1] * sin((2 * pi) * grid * j2)
      X[i2, ] <- X[i2, ] + e[2] * cos((2 * pi) * grid * j2)
    }
  }
  return(X)

  # X <- matrix(0, nrow = n_obs, ncol = length(grid))
  # for (i2 in 1:n_obs) {
  #   X[i2, ] <- X[i2, ] + rnorm(length(grid), 0, var1)
  #   X[i2, ] <- X[i2, ] + runif(1, 0, var2)
  #   X[i2, ] <- X[i2, ] + rnorm(1, 0, var3) * grid
  # 
  #   e <- abs(rnorm(2, 0, 1**(2)))
  #   X[i2, ] <- X[i2, ] + e[1] * sin((2 * pi) * grid)
  #   X[i2, ] <- X[i2, ] + e[2] * cos((2 * pi) * grid)
  # }
  # return(X)
}

X <- data_generation()
Y <- X %*% beta1 * (1 / n_var) + rnorm(n_obs, 0, 0.125)

Y <- as.numeric(Y)

final_X <- t(X)
X_train <- final_X[, -test]
Y_train <- Y[-test]
X_test <- final_X[, test]
Y_test <- Y[test]

largebasis <- create.bspline.basis(rangeval = c(0, length(grid)),
                                   nbasis = 100, 7)

harmFdpar <- fdPar(largebasis)

smooth_basis_train <- smooth.basis(y = X_train,
                                   fdParobj = largebasis)

smooth_basis_train_fd <- smooth_basis_train$fd

simulated_pcaObj_train <- pca.fd(smooth_basis_train_fd, 
                                 nharm = nharm, 
                                 harmfdPar = harmFdpar, 
                                 centerfns = TRUE)

test_in_train_basis <- smooth.basis(y = X_test, 
                                    fdParobj = simulated_pcaObj_train$harmonics)

xfdlist <- list(smooth_basis = smooth_basis_train_fd)
betalist <- list(smooth_basis_fd = smallbasis)
cMap <- smooth_basis_train$y2cMap # a matrix that contains the linear transformation that takes the raw data values into the coefficients defining a smooth functional data object

f_regress <- fRegress(y = Y_train, xfdlist, betalist, y2cMap = cMap)

prec <- predict.fRegress(
  object = f_regress,
  newdata = list(data = smooth_basis_test_fd)
  )

# #  plot the data and the fit
plot(x = Y_test, y = prec)
MSE_test <- mean((Y_test - prec)^2)

## Trying to calculate confidence intervals, but our models miss betastderrlist, which s needed to calculate fregress.stderr

## Tryto replicate from https://github.com/cran/fda/blob/master/demo/weatherANOVA.R line 62 ff"

yhatfd <- fd(f_regress$yhatfdobj)
argval <- (0:n_var) + 0.5
yhatmat <- eval.fd(argval, yhatfd)
ymat <- eval.fd(seq(0, n_var, 1), smooth_basis_train_fd)
tempresmat <- ymat[, 1:(n_obs - length(test))] - yhatmat[, 1:(n_obs - length(test))]
SigmaE <- var(t(tempresmat))

# #  plot covariance surface for errors

# par(mfrow=c(1,1))
# contour(SigmaE, xlab=\"Day\", ylab=\"Day\")
# lines(dayrange,dayrange,lty=4)

# #  plot standard deviation of errors

# par(mfrow=c(1,1), mar=c(5,5,3,2), pty=\"m\", ask=F)
# stddevE <- sqrt(diag(SigmaE))
# plot(daytime, stddevE, type=\"l\",
#      xlab=\"Day\", ylab=\"Standard error (deg C)\")

# #  Repeat regression, this time outputting results for
# #  confidence intervals

# stderrList <- fRegressStderr(fRegressList, tempy2cMap, SigmaE)

# betastderrlist <- stderrList$betastderrlist"

f_std <- fRegress.stderr(y = f_regress, 
                         y2cMap = cMap, 
                         SigmaE = f_regress$Cmatinv)

xfdlist <- list(smooth_basis = smooth_basis_train_fd)

betafdPar2 <- fdPar(simulated_pcaObj_train$harmonics)
betalist_pca <- list(smooth_basis_fd = betafdPar2)

data <- simulated_pcaObj_train$scores

f_regress_pca <- lm(Y_train ~ data)
(pca_coeffs <- summary(f_regress_pca)$coef)

betafd <- pca_coeffs[2, 1] * simulated_pcaObj_train$harmonics[1] +
  pca_coeffs[3, 1] * simulated_pcaObj_train$harmonics[2] +
  pca_coeffs[4, 1] * simulated_pcaObj_train$harmonics[3] +
  pca_coeffs[5, 1] * simulated_pcaObj_train$harmonics[4]

coefvar <- pca_coeffs[, 2]**2

betavar <- (pca_coeffs[2] * simulated_pcaObj_train$harmonics[1])**2 +
  (pca_coeffs[3] * simulated_pcaObj_train$harmonics[2])**2 +
  (pca_coeffs[4] * simulated_pcaObj_train$harmonics[3])**2 +
  (pca_coeffs[5] * simulated_pcaObj_train$harmonics[4])**2

#   This all has to become flexible! S.t we can choose the desired components that explain a
#   certain amount of the variance and then this will become dynamic

plot(betafd)
lines(betafd + 2 * sqrt(betavar), lty = 2, lwd = 1)
lines(betafd - 2 * sqrt(betavar), lty = 2, lwd = 1)

data <- simulated_pcaObj_test$scores
prec_pca <- predict.lm(object = f_regress_pca, data.frame(data = data))

prec_pca2 <- predict(object = f_regress_pca)

MSE_test <- mean((Y_test - prec_pca)^2)
print(MSE_test)
