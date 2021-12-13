# Clear Workspace
rm(list = ls())

# load libraries
library(fda)

# Define Data generator
data_generation <- function(fun) {
  # Variables for the data generation
  var1 <- 1
  var2 <- 2
  var3 <- 0.2

  # Container for observations
  X <- matrix(data = 0, nrow = n_obs, ncol = length(grid))

  for (i2 in 1:n_obs) {
    X[i2, ] <- {
      rnorm(n = length(grid), mean = 2, sd = var1)
      +runif(n = 1, min = 0, max = var2)
      +rnorm(n = 1, mean = 0, sd = var3) * grid
    } * fun

    for (j2 in 1:5) {
      e <- abs(rnorm(n = 2, mean = 0, sd = var1 / j2^(2)))
      X[i2, ] <- X[i2, ] + {
        e[1] * sin((2 * pi) * grid * j2)
        +e[2] * cos((2 * pi) * grid * j2)
      }
    }
  }
  return(X)
}

data(gasoline)
octane <- (gasoline$octane)
NIR <- as.matrix(gasoline$NIR)

###################################################################################################################
# set up "global" variables
###################################################################################################################

set.seed(100)

n_obs <- 60
nharm <- 4
n_var <- 400
grid <- seq(0, 1, length = n_var + 1) # why do we use exactly n_var = 1? (Jakob)

# smooth
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
# bumpy
f2 <- {
  1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2)
  -4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2)
  +8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2)
  -1 * exp(-0.5 * (grid - 0.8)^2 / 0.03^2)
}

# two different variances of error
sigma_eps_squared1_1 <- as.numeric((var(NIR %*% f1) / 0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 <- as.numeric((var(NIR %*% f1) / 0.6) - var(NIR %*% f1))
sigma_eps_squared2_1 <- as.numeric((var(NIR %*% f2) / 0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 <- as.numeric((var(NIR %*% f2) / 0.6) - var(NIR %*% f2))

# Function for simulation using b-splines
bspline_function <- function(rep) {
  CV_container_spline <- c()

  for (order in seq(8, 9, 1)) {
    for (j in seq(11, 12, 1)) {
      CV_container <- matrix(NaN, nrow = rep, ncol = 4)

      for (i in 1:rep) {
        # Generating simulation data from gasoline data set by adding different
        # error terms - matrix multiplication  yields a scalar response variable
        # each true beta and variance
        Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
        Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
        Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
        Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)

        data <- t(NIR)

        # print(dim(NIR))
        
        # Create b-spline basis of chosen specification
        smallbasis <- create.bspline.basis(
          rangeval = c(0, length(grid)), nbasis = as.numeric(j),
          norder = as.numeric(order)
        )
        
        # Convert data to basis representation and extract fd-object
        smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
        
        # Specify regressors as functional parameter objects
        xfdlist <- list(smooth_basis = smooth_basis_fd)
        
        ### As far as I can see (Jakob) we don't use these objects, do we?
        # betabasis1 <- create.constant.basis(c(0, 60))
        # betafd1 <- fd(0, betabasis1)
        # betafdPar1 <- fdPar(betafd1)
        # betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
        
        # Extract parameters in the form used by fRegress
        betafdPar2 <- fdPar(smallbasis)
        betalist <- list(smooth_basis_fd_data = betafdPar2)
  
        # Perform functional regressions including the calculation
        # of CV-scores using 
        f_regress1_1 <- fRegress.CV(y = Y1_1, xfdlist, betalist)
        f_regress1_2 <- fRegress.CV(y = Y1_2, xfdlist, betalist)
        f_regress2_1 <- fRegress.CV(y = Y2_1, xfdlist, betalist)
        f_regress2_2 <- fRegress.CV(y = Y2_2, xfdlist, betalist)

        CV_container[i, 1] <- f_regress1_1$SSE.CV
        CV_container[i, 2] <- f_regress1_2$SSE.CV
        CV_container[i, 3] <- f_regress2_1$SSE.CV
        CV_container[i, 4] <- f_regress2_2$SSE.CV
      }
      scaled_MSE <- colMeans(CV_container)
      scaled_MSE[5] <- j
      scaled_MSE[6] <- order
      CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    }
  }
  colnames(CV_container_spline) <- c("f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline", "n_basis", "n_order")
  return(CV_container_spline)
}

fpcr_function <- function(rep) {
  CV_container_spline <- c()
  train.control <- trainControl(method = "cv", number = 10)
  for (order in seq(8, 9, 1)) {
    for (j in seq(11, 12, 1)) {
      CV_container <- matrix(NaN, nrow = rep, ncol = 5)

      for (i in 1:rep) {
        # each true beta and variance
        Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
        Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
        Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
        Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)

        data <- t(NIR)

        # print(dim(NIR))
        smallbasis <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order))
        smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd

        simulated_pcaObj <- pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)


        dataframe1_1 <- data.frame(Y1_1, simulated_pcaObj$scores)
        CV_container[i, 1] <- train(Y1_1 ~ ., data = dataframe1_1, method = "lm", trControl = train.control)$results[[2]]

        dataframe1_2 <- data.frame(Y1_2, simulated_pcaObj$scores)
        CV_container[i, 2] <- train(Y1_2 ~ ., data = dataframe1_2, method = "lm", trControl = train.control)$results[[2]]

        dataframe2_1 <- data.frame(Y2_1, simulated_pcaObj$scores)
        CV_container[i, 3] <- train(Y2_1 ~ ., data = dataframe2_1, method = "lm", trControl = train.control)$results[[2]]

        dataframe2_2 <- data.frame(Y2_2, simulated_pcaObj$scores)
        CV_container[i, 4] <- train(Y2_2 ~ ., data = dataframe2_2, method = "lm", trControl = train.control)$results[[2]]

        CV_container[i, 5] <- sum(simulated_pcaObj$varprop)
      }
      scaled_MSE <- colMeans(CV_container)
      scaled_MSE[6] <- j
      scaled_MSE[7] <- order
      CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    }
  }
  colnames(CV_container_spline) <- c("f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr", "varprop", "n_basis", "n_order")
  return(CV_container_spline)
}
