suppressMessages(library(refund))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2)) 
suppressMessages(library(fda))
#suppressMessages(library(fdaACF))
suppressMessages(library(fpca))
suppressMessages(library(caret))


source("C:/Users/Jonathan/Desktop/RM_Stats/Code/data_generator.R") 
source("C:/Users/Jonathan/Desktop/RM_Stats/Code/auxiliary_functions.R") 
set.seed(100)

###################################################

data(gasoline)
octane <- (gasoline$octane)
NIR    <- as.matrix(gasoline$NIR)

###################################################################################################################
# set up "global" variables
###################################################################################################################
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
sigma_eps_squared1_1 = as.numeric((var(NIR %*% f1)/0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 = as.numeric((var(NIR %*% f1)/0.6) - var(NIR %*% f1) )
sigma_eps_squared2_1 = as.numeric((var(NIR %*% f2)/0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 = as.numeric((var(NIR %*% f2)/0.6) - var(NIR %*% f2) )



bspline_function <- function(rep, NIR, n_obs){
  CV_container_spline  <- c()
  for(j in seq(5,25,1)){
    
    CV_container  <- matrix(NaN, nrow = rep, ncol = 4)
    success_count = 0
    for(i in 1 : rep){
      print(success_count)
      tryCatch({
        Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
        Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
        Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
        Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
        
        data = t(NIR)
        
        #print(dim(NIR))
        smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = 4)
        smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
        xfdlist = list(smooth_basis=smooth_basis_fd)
        betabasis1 <- create.constant.basis(c(0, 60))
        betafd1    <- fd(0, betabasis1)
        betafdPar1 <- fdPar(betafd1)
        #betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
        betafdPar2  <- fdPar(smallbasis)
        betalist <- list(smooth_basis_fd_data=betafdPar2)
        
        f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist, betalist)
        f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist, betalist)
        f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist, betalist)
        f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist, betalist)
        
        CV_container[i,1] <- f_regress1_1$SSE.CV
        CV_container[i,2] <- f_regress1_2$SSE.CV
        CV_container[i,3] <- f_regress2_1$SSE.CV
        CV_container[i,4] <- f_regress2_2$SSE.CV
        success_count = success_count +1
      }, error = function(e){print("not succesfull!")})
      
    }
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[5] = j
    scaled_MSE[6] = success_count
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    
  }
  colnames(CV_container_spline) = c("f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline", "n_basis", "success_count")
  return(CV_container_spline)
}


fourier_function <- function(rep, NIR, n_obs){
  CV_container_spline  <- c()
  
  for(j in seq(1,25,1)){
    success_count = 0
    CV_container  <- matrix(NaN, nrow = rep, ncol = 4)
    
    for(i in 1 : rep){
      print(success_count)
      
      tryCatch({
      #each true beta and variance
      Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
      Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
      Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
      Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
      
      data = t(NIR)
      
      #print(dim(NIR))
      smallbasis      <- create.fourier.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j))
      smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
      xfdlist = list(smooth_basis=smooth_basis_fd)
      betabasis1 <- create.constant.basis(c(0, 60))
      betafd1    <- fd(0, betabasis1)
      betafdPar1 <- fdPar(betafd1)
      #betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
      betafdPar2  <- fdPar(smallbasis)
      betalist <- list(smooth_basis_fd_data=betafdPar2)
      
      f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist, betalist)
      f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist, betalist)
      f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist, betalist)
      f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist, betalist)
      
      CV_container[i,1] <- f_regress1_1$SSE.CV
      CV_container[i,2] <- f_regress1_2$SSE.CV
      CV_container[i,3] <- f_regress2_1$SSE.CV
      CV_container[i,4] <- f_regress2_2$SSE.CV
      success_count = success_count +1
      }, error = function(e){print("not succesfull!")})
    }
    
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[5] = j
    scaled_MSE[6] = success_count
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    
  }
  colnames(CV_container_spline) = c("f1_e1_fourier", "f1_e2_fourier", "f2_e1_fourier", "f2_e2_fourier", "n_basis", "success_count")
  return(CV_container_spline)
}

fpcr_function <- function(rep, NIR, n_obs){
  CV_container_spline  <- c()
  train.control <- caret::trainControl(method = "cv", number = 10)
  for(j in seq(5,25,1)){
    success_count = 0
    CV_container  <- matrix(NaN, nrow = rep, ncol = 5)
    
    for(i in 1 : rep){
      tryCatch({
        print(success_count)
        #each true beta and variance
        Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
        Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
        Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
        Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
        
        data = t(NIR)
        
        #print(dim(NIR))
        smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = 4)
        smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
        
        simulated_pcaObj = pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)
        
        
        dataframe1_1 = data.frame(Y1_1, simulated_pcaObj$scores)
        CV_container[i,1] <- train(Y1_1 ~., data = dataframe1_1, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe1_2 = data.frame(Y1_2, simulated_pcaObj$scores)
        CV_container[i,2] <- train(Y1_2 ~., data = dataframe1_2, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe2_1 = data.frame(Y2_1, simulated_pcaObj$scores)
        CV_container[i,3] <- train(Y2_1 ~., data = dataframe2_1, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe2_2 = data.frame(Y2_2, simulated_pcaObj$scores)
        CV_container[i,4] <- train(Y2_2 ~., data = dataframe2_2, method = "lm", trControl = train.control)$results[[2]]
        
        CV_container[i,5] <- sum(simulated_pcaObj$varprop)
        
        success_count = success_count +1
      }, error = function(e){print("not succesfull!")})
    }
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[6] = j
    scaled_MSE[7] = success_count
    
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    
  }
  colnames(CV_container_spline) = c("f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr","varprop" ,"n_basis", "success_count")
  return(CV_container_spline)
}


fpcr_fourier_function <- function(rep, NIR, n_obs){
  CV_container_spline  <- c()
  train.control <- caret::trainControl(method = "cv", number = 10)
  for(j in seq(1,25,1)){
    success_count = 0
    CV_container  <- matrix(NaN, nrow = rep, ncol = 5)
    
    for(i in 1 : rep){
      print(success_count)
      
      tryCatch({
      
        #each true beta and variance
        Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
        Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
        Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
        Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
        
        data = t(NIR)
        
        #print(dim(NIR))
        smallbasis      <- create.fourier.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j))
        smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
        
        simulated_pcaObj = pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)
        
        
        dataframe1_1 = data.frame(Y1_1, simulated_pcaObj$scores)
        CV_container[i,1] <- train(Y1_1 ~., data = dataframe1_1, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe1_2 = data.frame(Y1_2, simulated_pcaObj$scores)
        CV_container[i,2] <- train(Y1_2 ~., data = dataframe1_2, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe2_1 = data.frame(Y2_1, simulated_pcaObj$scores)
        CV_container[i,3] <- train(Y2_1 ~., data = dataframe2_1, method = "lm", trControl = train.control)$results[[2]]
        
        dataframe2_2 = data.frame(Y2_2, simulated_pcaObj$scores)
        CV_container[i,4] <- train(Y2_2 ~., data = dataframe2_2, method = "lm", trControl = train.control)$results[[2]]
        
        CV_container[i,5] <- sum(simulated_pcaObj$varprop)
        success_count = success_count +1
      }, error = function(e){print("not succesfull!")})
      
    }
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[6] = j
    scaled_MSE[7] = success_count
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
    
  }
  colnames(CV_container_spline) = c("f1_e1_fpcr_fourier", "f1_e2_fpcr_fourier", "f2_e1_fpcr_fourier", "f2_e2_fpcr_fourier","varprop" ,"n_basis", "success_count")
  return(CV_container_spline)
}


#############
generated_curves = NIR_curve_generator(n=1000, n_harmonics = 30)
input_data = as.matrix(generated_curves[,-1])
#############

#############################################################################
###             Run function using the NIR data / generated data          ###
#############################################################################


test_bspline_function = bspline_function(1000, input_data, 1000)
#write.table(test_bspline_function,file="Results/test_bspline_expansion_1000.csv")
test_bspline_function

test_fourier_function = fourier_function(1000, input_data, 1000)
#write.table(test_fourier_function,file="Results/test_fourier_expansion_1000.csv")
test_fourier_function

test_fpcr = fpcr_function(1000, input_data, 1000)
#write.table(test_fpcr,file="Results/test_fpcr_bsplines_1000.csv")
test_fpcr

test_fpcr2 = fpcr_fourier_function(1000, input_data, 1000)
write.table(test_fpcr2,file="Results/test_fpcr_fourier_1000.csv")
test_fpcr2

test_bspline_function = bspline_function(1000, NIR, 60)
#write.table(test_bspline_function,file="Results/test_bspline_expansion_NIR.csv")
test_bspline_function

test_fourier_function = fourier_function(1000, NIR, 60)
#write.table(test_fourier_function,file="Results/test_fourier_expansion_NIR.csv")
test_fourier_function

test_fpcr = fpcr_function(1000, NIR, 60)
#write.table(test_fpcr,file="Results/test_fpcr_bsplines_NIR.csv")
test_fpcr

test_fpcr2 = fpcr_fourier_function(1000, NIR, 60)
write.table(test_fpcr2,file="Results/test_fpcr_fourier_NIR.csv")
test_fpcr2

#######################################
###             Test-Area          ###
#######################################
data = t(NIR)
Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
smallbasis      <- create.fourier.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(4))
smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
xfdlist = list(smooth_basis=smooth_basis_fd)
betabasis1 <- create.constant.basis(c(0, 60))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
#betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
betafdPar2  <- fdPar(smallbasis)
betalist <- list(smooth_basis_fd_data=betafdPar2)
f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist, betalist)



f_regress1_1original <- fRegress.CVoriginal(y = Y1_1, xfdlist, betalist)
f_regress1_1original$SSE.CV



f_regressn<- fRegress(y = Y1_1, xfdlist, betalist)

