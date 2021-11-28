rm(list=ls())

data_generation <-function(fun){
    var1 = 1
    var2 = 2
    var3 = 0.2

    X <- matrix(0, nrow=n_obs, ncol=length(grid))
    for(i2 in 1:n_obs){
        X[i2,]=X[i2,]+rnorm(length(grid), 2, var1)
        X[i2,]=X[i2,]+runif(1, 0, var2)
        X[i2,]=X[i2,]+rnorm(1, 0, var3)*grid
        X[i2,]=X[i2,]*fun

        for(j2 in 1:5){
           e =abs(rnorm(2, 0, var1/j2^(2)))
           X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
           X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
        }
        
    }
    return(X)
}

data(gasoline)
octane <- (gasoline$octane)
NIR    <- as.matrix(gasoline$NIR)
test_size = 20

# set up "global" variables
set.seed(100)
n_obs = 60
nharm = 4
n_var = 400
grid = seq(0, 1, length = n_var+1)
#smooth
f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
#bumpy
f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -  4*exp(-0.5*(grid-0.45)^2/0.015^2) +  8*exp(-0.5*(grid-0.6)^2/0.02^2) -  exp(-0.5*(grid-0.8)^2/0.03^2)

#two different variances of error
#sigma_eps_squared1 = as.numeric(2/3*var(NIR %*% f1))
#sigma_eps_squared2 = as.numeric(1/9*var(NIR %*% f1))
sigma_eps_squared1_1 = as.numeric((var(NIR %*% f1)/0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 = as.numeric((var(NIR %*% f1)/0.6) - var(NIR %*% f1) )
sigma_eps_squared2_1 = as.numeric((var(NIR %*% f2)/0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 = as.numeric((var(NIR %*% f2)/0.6) - var(NIR %*% f2) )

rep      <- 500
MSE      <- matrix(NaN, nrow = rep, ncol = 9)
MSE_avrg = c()


for(j in seq(10,20,1)){
  for(i in 1 : rep){

    test = sample(1:n_obs, size = test_size, replace = FALSE)

    #each true beta and variance
    Y1_1 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared1_1)
    Y1_2 <- NIR %*% f1 + rnorm(n_obs, 0, sigma_eps_squared1_2)
    Y2_1 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared2_1)
    Y2_2 <- NIR %*% f2 + rnorm(n_obs, 0, sigma_eps_squared2_2)
    
    Y1_1 <- as.numeric(Y1_1)
    Y1_2 <- as.numeric(Y1_2)
    Y2_1 <- as.numeric(Y2_1)
    Y2_2 <- as.numeric(Y2_2)

    Y1_1train <- Y1_1[-test]
    Y1_2train <- Y1_2[-test]
    Y2_1train <- Y2_1[-test]
    Y2_2train <- Y2_2[-test]

    Y1_1test <- Y1_1[test]
    Y1_2test <- Y1_2[test]
    Y2_1test <- Y2_1[test]
    Y2_2test <- Y2_2[test]
    data = t(NIR)
    X_train = data[,-test]
    X_test = data[,test]

    #print(dim(NIR))
    smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, 8)
    smooth_basis_fd_train <- smooth.basis(y = X_train, fdParobj=smallbasis)$fd
    smooth_basis_fd_test <- smooth.basis(y = X_test, fdParobj=smallbasis)$fd

    xfdlist = list(smooth_basis=smooth_basis_fd_train)
    betabasis1 <- create.constant.basis(c(0, (60-length(test))))
    betafd1    <- fd(0, betabasis1)
    betafdPar1 <- fdPar(betafd1)
    betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, 8)
    betafdPar2  <- fdPar(betafd2)
    betalist <- list(smooth_basis_fd=betafdPar2)

    f_regress1_1 <- fRegress(y = Y1_1train, xfdlist, betalist)
    f_regress1_2 <- fRegress(y = Y1_2train, xfdlist, betalist)
    f_regress2_1 <- fRegress(y = Y2_1train, xfdlist, betalist)
    f_regress2_2 <- fRegress(y = Y2_2train, xfdlist, betalist)
    
    prediction1_1 <- predict.fRegress(object = f_regress1_1, newdata = list(data = smooth_basis_fd_test))
    prediction1_2 <- predict.fRegress(object = f_regress1_2, newdata = list(data = smooth_basis_fd_test))
    prediction2_1 <- predict.fRegress(object = f_regress2_1, newdata = list(data = smooth_basis_fd_test))
    prediction2_2 <- predict.fRegress(object = f_regress2_2, newdata = list(data = smooth_basis_fd_test))
    
    MSE[i,1] <- mean((prediction1_1 - Y1_1test)^2)
    MSE[i,2] <- mean((prediction1_2 - Y1_2test)^2)
    MSE[i,3] <- mean((prediction2_1 - Y2_1test)^2)
    MSE[i,4] <- mean((prediction2_2 - Y2_2test)^2)
    
    

################################################################################

    simulated_pcaObj_train = pca.fd(smooth_basis_fd_train, nharm = nharm,  centerfns = TRUE)
    simulated_pcaObj_test = pca.fd(smooth_basis_fd_test, nharm = nharm, centerfns = TRUE)
    data = simulated_pcaObj_train$scores
    f_regress_pca <- lm(Y1_1train ~ data)
    data = simulated_pcaObj_test$scores
    prediction1_1 <- predict.lm(object = f_regress_pca, data.frame(data = data))

    data = simulated_pcaObj_train$scores
    f_regress_pca <- lm(Y1_2train ~ data)
    data = simulated_pcaObj_test$scores
    prediction1_2 <- predict.lm(object = f_regress_pca, data.frame(data = data))

    data = simulated_pcaObj_train$scores
    f_regress_pca <- lm(Y2_1train ~ data)
    data = simulated_pcaObj_test$scores
    prediction2_1 <- predict.lm(object = f_regress_pca, data.frame(data = data))

    data = simulated_pcaObj_train$scores
    f_regress_pca <- lm(Y2_2train ~ data)
    data = simulated_pcaObj_test$scores
    prediction2_2 <- predict.lm(object = f_regress_pca, data.frame(data = data))
    
    MSE[i,5] <- mean((prediction1_1 - Y1_1test)^2)
    MSE[i,6] <- mean((prediction1_2 - Y1_2test)^2)
    MSE[i,7] <- mean((prediction2_1 - Y2_1test)^2)
    MSE[i,8] <- mean((prediction2_2 - Y2_2test)^2)
    
    
  }#end of reps

    scaled_MSE <- colMeans(MSE)
    # scalling by var(X%*% tur_fun)
    for (i in seq(1,7,2)){
        scaled_MSE[i] = scaled_MSE[i] / var(NIR %*% f1)
    }
    for (i in seq(2,8,2)){
        scaled_MSE[i] = scaled_MSE[i] / var(NIR %*% f2)
    }
    scaled_MSE[9] = j
    
    
  
    MSE_avrg <- rbind(MSE_avrg, scaled_MSE)
  
}
colnames(MSE_avrg) = c("f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline", "f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr", "n_basis")
MSE_avrg