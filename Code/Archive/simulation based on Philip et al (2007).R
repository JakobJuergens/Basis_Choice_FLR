rm(list=ls())

data(gasoline)
octane <- (gasoline$octane)
NIR    <- as.matrix(gasoline$NIR)

# set up "global" variables
set.seed(100)
n_obs = 60
n_var = 400
grid = seq(0, 1, length = n_var+1)

#smooth
f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
#bumpy
f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -
  4*exp(-0.5*(grid-0.45)^2/0.015^2) +
  8*exp(-0.5*(grid-0.6)^2/0.02^2) -
  exp(-0.5*(grid-0.8)^2/0.03^2)

#two different variances of error
sigma_eps_squared1 = as.numeric(2/3*var(NIR %*% f1))
sigma_eps_squared2 = as.numeric(1/9*var(NIR %*% f1))

rep      <- 300
MSE      <- matrix(0, nrow = rep, ncol = 4)
MSE_avrg <- matrix(0, nrow = 5, ncol = 4)

for(j in c(5, 12, 20 ,35, 50)){
  for(i in 1 : rep){
    #each true beta and variance
    Y1_1 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared1)
    Y1_2 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared2)
    Y2_1 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared1)
    Y2_2 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared2)
    
    Y1_1 <- as.numeric(Y1_1)
    Y1_2 <- as.numeric(Y1_2)
    Y2_1 <- as.numeric(Y2_1)
    Y2_2 <- as.numeric(Y2_2)
    
    smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, 5)
    smooth_basis_fd <- smooth.basis(y = t(NIR), fdParobj=smallbasis)$fd
    
    f_regress1_1 <- fRegress(Y1_1 ~ smooth_basis_fd)
    f_regress1_2 <- fRegress(Y1_2 ~ smooth_basis_fd)
    f_regress2_1 <- fRegress(Y2_1 ~ smooth_basis_fd)
    f_regress2_2 <- fRegress(Y2_2 ~ smooth_basis_fd)
    
    prediction1_1 <- as.numeric(f_regress1_1$yhatfdobj)
    prediction1_2 <- as.numeric(f_regress1_2$yhatfdobj)
    prediction2_1 <- as.numeric(f_regress2_1$yhatfdobj)
    prediction2_2 <- as.numeric(f_regress2_2$yhatfdobj)
    
    MSE[i,1] <- mean((prediction1_1 - as.numeric(octane))^2)
    MSE[i,2] <- mean((prediction1_2 - as.numeric(octane))^2)
    MSE[i,3] <- mean((prediction2_1 - as.numeric(octane))^2)
    MSE[i,4] <- mean((prediction2_2 - as.numeric(octane))^2)
    
    scaled_MSE <- colMeans(MSE/as.numeric(var(NIR %*% f1)))
  }
  if(j == 5){
    MSE_avrg <- scaled_MSE
  }
  else{
    MSE_avrg <- rbind(MSE_avrg, scaled_MSE)
  }
}

MSE_avrg


#######
# Additions by Jona, try it out!
#####

# rm(list=ls())

# data(gasoline)
# octane <- (gasoline$octane)
# NIR    <- as.matrix(gasoline$NIR)
# test = seq(1,20,1)


# # set up "global" variables
# set.seed(100)
# n_obs = 60
# n_var = 400
# grid = seq(0, 1, length = n_var+1)

# #smooth
# f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
# #bumpy
# f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -
#   4*exp(-0.5*(grid-0.45)^2/0.015^2) +
#   8*exp(-0.5*(grid-0.6)^2/0.02^2) -
#   exp(-0.5*(grid-0.8)^2/0.03^2)

# #two different variances of error
# sigma_eps_squared1 = as.numeric(2/3*var(NIR %*% f1))
# sigma_eps_squared2 = as.numeric(1/9*var(NIR %*% f1))

# rep      <- 300
# MSE      <- matrix(NaN, nrow = rep, ncol = 4)
# MSE_avrg = c()

# for(j in c(12, 20, 25, 30 ,35)){
#   for(i in 1 : rep){
#     #each true beta and variance
#     Y1_1 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared1)
#     Y1_2 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared2)
#     Y2_1 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared1)
#     Y2_2 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared2)
    
#     Y1_1 <- as.numeric(Y1_1)
#     Y1_2 <- as.numeric(Y1_2)
#     Y2_1 <- as.numeric(Y2_1)
#     Y2_2 <- as.numeric(Y2_2)

#     Y1_1train <- Y1_1[-test]
#     Y1_2train <- Y1_2[-test]
#     Y2_1train <- Y2_1[-test]
#     Y2_2train <- Y2_2[-test]

#     Y1_1test <- Y1_1[test]
#     Y1_2test <- Y1_2[test]
#     Y2_1test <- Y2_1[test]
#     Y2_2test <- Y2_2[test]
#     data = t(NIR)
#     X_train = data[,-test]
#     X_test = data[,test]

#     #print(dim(NIR))
#     smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, 5)
#     smooth_basis_fd_train <- smooth.basis(y = X_train, fdParobj=smallbasis)$fd
#     smooth_basis_fd_test <- smooth.basis(y = X_test, fdParobj=smallbasis)$fd

#     xfdlist = list(smooth_basis=smooth_basis_fd_train)
#     betabasis1 <- create.constant.basis(c(0, (60-length(test))))
#     betafd1    <- fd(0, betabasis1)
#     betafdPar1 <- fdPar(betafd1)
#     betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, 5)
#     betafdPar2  <- fdPar(betafd2)
#     betalist <- list(smooth_basis_fd=betafdPar2)

#     f_regress1_1 <- fRegress(y = Y1_1train, xfdlist, betalist)
#     f_regress1_2 <- fRegress(y = Y1_2train, xfdlist, betalist)
#     f_regress2_1 <- fRegress(y = Y2_1train, xfdlist, betalist)
#     f_regress2_2 <- fRegress(y = Y2_2train, xfdlist, betalist)
    
#     prediction1_1 <- predict.fRegress(object = f_regress1_1, newdata = list(data = smooth_basis_fd_test))
#     prediction1_2 <- predict.fRegress(object = f_regress1_2, newdata = list(data = smooth_basis_fd_test))
#     prediction2_1 <- predict.fRegress(object = f_regress2_1, newdata = list(data = smooth_basis_fd_test))
#     prediction2_2 <- predict.fRegress(object = f_regress2_2, newdata = list(data = smooth_basis_fd_test))
    
#     MSE[i,1] <- mean((prediction1_1 - as.numeric(Y1_1test))^2)
#     MSE[i,2] <- mean((prediction1_2 - as.numeric(Y1_2test))^2)
#     MSE[i,3] <- mean((prediction2_1 - as.numeric(Y2_1test))^2)
#     MSE[i,4] <- mean((prediction2_2 - as.numeric(Y2_2test))^2)
    
#     scaled_MSE <- colMeans(MSE/as.numeric(var(NIR %*% f1)))
#   }#end of reps

  
#     MSE_avrg <- rbind(MSE_avrg, scaled_MSE)
  
# }

# MSE_avrg
