rm(list=ls())

suppressMessages(library(refund))
suppressMessages(library(MASS))
suppressMessages(library(tidyverse))
suppressMessages(library(fda))
suppressMessages(library(fpca))
suppressMessages(library(caret))

### load data set
data(gasoline)
octane <- (gasoline$octane)
NIR <- as.matrix(gasoline$NIR)

##### Define Simulation Parameters #####
### set up "global" variables
n_obs <- 60
nharm <- 4
n_var <- 400
grid <- seq(0, 1, length = n_var + 1)

### set up coefficient "functions" / error terms
# smooth
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
# bumpy
f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) - exp(-0.5 * (grid - 0.8)^2 / 0.03^2)
# two different variances of error
sigma_eps_squared1_1 <- as.numeric((var(NIR %*% f1) / 0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 <- as.numeric((var(NIR %*% f1) / 0.6) - var(NIR %*% f1))
sigma_eps_squared2_1 <- as.numeric((var(NIR %*% f2) / 0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 <- as.numeric((var(NIR %*% f2) / 0.6) - var(NIR %*% f2))


### load code from other source files
source("data_generator.R")
source("k_fold_CV_function.R")


### FPCA simulation function - bspline basis
fpcr_function <- function(rep, my_data = NULL, n_obs, nharm, seed, debug = FALSE) {
  set.seed(seed)
  ### set up coefficient "functions" / error terms
  # smooth
  f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
  # bumpy
  f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) - exp(-0.5 * (grid - 0.8)^2 / 0.03^2)
  # two different variances of error
  sigma_eps_squared1_1 <- as.numeric((var(NIR %*% f1) / 0.9) - var(NIR %*% f1))
  sigma_eps_squared1_2 <- as.numeric((var(NIR %*% f1) / 0.6) - var(NIR %*% f1))
  sigma_eps_squared2_1 <- as.numeric((var(NIR %*% f2) / 0.9) - var(NIR %*% f2))
  sigma_eps_squared2_2 <- as.numeric((var(NIR %*% f2) / 0.6) - var(NIR %*% f2))
  
  # check if data (my_data) is provided
  if (is.null(my_data)) {
    data_cond <- FALSE
  } else {
    data_cond <- TRUE
  }
  
  # specfiy number of basis functions that should be considered
  n_basis <- seq(from = 4, to = 25, by = 1)
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 7)
  )
  
  colnames(CV_container) <- c(
    "f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr",
    "varprop", "n_basis", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))
  
  # create basis functions
  basis_functions <- map(
    .x = n_basis,
    .f = function(j) create.bspline.basis(rangeval = c(0, 1), nbasis = j, norder = 4)
  )
  
  
  # loop over repetitions
  for (i in 1:rep) {
    if (debug) {
      print(paste0("Repetition: ", i, " out of ", rep))
    }
    
    # if no data (my_data) is provided generate new curves
    if (data_cond == FALSE) {
      generated_curves <- NIR_curve_generator(n = n_obs, n_harmonics = 30)
      my_data <- as.matrix(generated_curves[, -1])
    }
    
    # transpose data for use with fda package
    data <- t(my_data)
    
    # calculate responses
    Y1_1 <- as.numeric(my_data %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
    Y1_2 <- as.numeric(my_data %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
    Y2_1 <- as.numeric(my_data %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
    Y2_2 <- as.numeric(my_data %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
    
    # df1_1 <- data.frame(Y1_1, my_data)
    # df1_2 <- data.frame(Y1_2, my_data)
    # df2_1 <- data.frame(Y2_1, my_data)
    # df2_2 <- data.frame(Y2_2, my_data)
    
    # loop over repetitions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      smallbasis <- basis_functions[[j]]
      
      tryCatch(
        {
          MSPE_mat <- as.data.frame(matrix(0, nrow = 10, ncol = 4))
          colnames(MSPE_mat) <- c('MSPE1_1', 'MSPE1_2', 'MSPE2_1', 'MSPE2_2')
          
          for(m in 1 : 10){
            sampling <- sort(sample(1:n_obs, n_obs/10, replace = FALSE))
            
            # test_data <- data[, sampling]
            # train_data <- data[, -sampling]
            
            # df1_1_train <- df1_1[-sampling, ]
            # df1_2_train <- df1_2[-sampling, ]
            # df2_1_train <- df2_1[-sampling, ]
            # df2_2_train <- df2_2[-sampling, ]
            # 
            # df1_1_test <- df1_1[sampling, ]
            # df1_2_test <- df1_1[sampling, ]
            # df2_1_test <- df1_1[sampling, ]
            # df2_2_test <- df1_1[sampling, ]
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd
            smooth_basis_fd_train <- smooth_basis_fd[-sampling]
            smooth_basis_fd_test  <- smooth_basis_fd[sampling]
            
            # perform fPCA
            train_fd <- pca.fd(fdobj = smooth_basis_fd_train, 
                               nharm = nharm, centerfns = TRUE)
            
            # extract scores for training set
            train_scores <- train_fd$scores
            
            scores_vec  <- c()
            scores_mat  <- matrix(0, nrow = n_obs/10, ncol = nharm)
            
            tmp_grid <- seq(0, 1, length.out = 1000)
            basis_eval <- eval.basis(evalarg = tmp_grid, basisobj = smallbasis)
            integral_matrix <- matrix(data = 0, nrow = n_basis[j], ncol = n_basis[j])
            
            for(l in 1:n_basis[j]){
              for(k in 1:n_basis[j]){
                integral_matrix[l,k] <- integrate(function(t){
                    approx(x = tmp_grid,
                           y = basis_eval[, l] * t(basis_eval[, k]), xout = t)$y}, lower = 0, upper = 1
                  )$value
              }
            }
            
            sc <-  function(x){
              scores_vec <- c()
              
              for(k in 1 : nharm){
                newfun <- (smooth_basis_fd_test[x]$coefs - train_fd$meanfd$coefs) %*% t(train_fd$harmonics[k]$coefs)
                
                weighted_basis_integrals_products <- newfun * integral_matrix
                
                estimated_score <- sum(weighted_basis_integrals_products)
                  
                scores_vec    <- c(scores_vec, estimated_score)
              }
              return(scores_vec)
            }
            
            for (p in 1:(n_obs/10)){
              scores_mat[p,] <- sc(p)
            }
            
            scores_tib <- as_tibble(scores_mat)
            names(scores_tib) <- paste0('harm_', 1:nharm)
            
            # # combine objects into data frames for linear regression
            dataframe1_1 <- as_tibble(cbind(df1_1_train[,1], train_fd$scores))
            dataframe1_2 <- as_tibble(cbind(df1_2_train[,1], train_fd$scores))
            dataframe2_1 <- as_tibble(cbind(df2_1_train[,1], train_fd$scores))
            dataframe2_2 <- as_tibble(cbind(df2_2_train[,1], train_fd$scores))
            
            # give correct names
            names(dataframe1_1) <- c('Y1_1', paste0('harm_', 1:nharm))
            names(dataframe1_2) <- c('Y1_2', paste0('harm_', 1:nharm))
            names(dataframe2_1) <- c('Y2_1', paste0('harm_', 1:nharm))
            names(dataframe2_2) <- c('Y2_2', paste0('harm_', 1:nharm))
            
            model1_1 <- lm(Y1_1 ~., dataframe1_1) 
            model1_2 <- lm(Y1_2 ~., dataframe1_2)
            model2_1 <- lm(Y2_1 ~., dataframe2_1)
            model2_2 <- lm(Y2_2 ~., dataframe2_2)
            
            test_prediction1_1 <- unname(predict(object = model1_1, newdata = scores_tib))
            test_prediction1_2 <- unname(predict(object = model1_2, newdata = scores_tib))
            test_prediction2_1 <- unname(predict(object = model2_1, newdata = scores_tib))
            test_prediction2_2 <- unname(predict(object = model2_2, newdata = scores_tib))
            
            rmse_test1_1 <- sqrt(1/10 * sum((test_prediction1_1 - Y1_1)^2))
            rmse_test1_2 <- sqrt(1/10 * sum((test_prediction1_2 - Y1_2)^2))
            rmse_test2_1 <- sqrt(1/10 * sum((test_prediction2_1 - Y2_1)^2))
            rmse_test2_2 <- sqrt(1/10 * sum((test_prediction2_2 - Y2_2)^2))
            
            # MSPE1_1 <- sum((model1_1$residuals[(n_obs * 9 / 10 + 1):n_obs])^2)/(n_obs/10)
            # MSPE1_2 <- sum((model1_2$residuals[(n_obs * 9 / 10 + 1):n_obs])^2)/(n_obs/10)
            # MSPE2_1 <- sum((model2_1$residuals[(n_obs * 9 / 10 + 1):n_obs])^2)/(n_obs/10)
            # MSPE2_2 <- sum((model2_2$residuals[(n_obs * 9 / 10 + 1):n_obs])^2)/(n_obs/10)
            
            MSPE_mat[m, 1] <- rmse_test1_1
            MSPE_mat[m, 2] <- rmse_test1_2
            MSPE_mat[m, 3] <- rmse_test2_1
            MSPE_mat[m, 4] <- rmse_test2_2
          }
          
          MSPE_mat <- colMeans(MSPE_mat)
          
          # generate tmp variable for current number of successful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_mat[1] * 1 / (tmp_sr + 1)
          
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_mat[2] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_mat[3] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_mat[4] * 1 / (tmp_sr + 1)
          
          # update varprop
          CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(train_fd$varprop) * 1 / (tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", j, " bspline basis functions for FPCA Run.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}


test_fpcr_NIR2 <- fpcr_function(
  rep = 100, my_data = NIR, n_obs = 60, nharm = 2, seed = 100, debug = TRUE
)
test_fpcr_NIR3 <- fpcr_function(
  rep = 100, my_data = NIR, n_obs = 60, nharm = 3, seed = 100, debug = TRUE
)
test_fpcr_NIR4 <- fpcr_function(
  rep = 100, my_data = NIR, n_obs = 60, nharm = 4, seed = 100, debug = TRUE
)