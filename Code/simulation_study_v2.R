##### Set up R session for simulation #####
### Clear Workspace
rm(list = ls())

### load libraries
suppressMessages(library(refund))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(fda))
# suppressMessages(library(fdaACF))
suppressMessages(library(fpca))
suppressMessages(library(caret))
suppressMessages(library(parallel))

### load code from other source files
source("data_generator.R")
source("auxiliary_functions.R")

### set seed for reproducibility
set.seed(100)

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

##### Define Simulation Functions #####
### b_spline simulation function
bspline_function <- function(rep, NIR = NULL, n_obs, debug = FALSE) {
  # check if data (NIR) is provided
  if (is.null(NIR)) {
    data_cond <- FALSE
  } else {
    data_cond <- TRUE
  }
  
  # specfiy number of basis functions that should be considered
  n_basis <- seq(from = 5, to = 25, by = 1)
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 6)
  )
  
  colnames(CV_container) <- c(
    "f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline",
    "n_basis", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))
  
  # create basis functions
  basis_functions <- map(
    .x = n_basis,
    .f = function(j) create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, norder = 4)
  )
  
  # prepare objects for functional linear regression
  betafdPar2_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) fdPar(basis_functions[[i]])
  )
  
  betalist_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) list(smooth_basis_fd_data = betafdPar2_list[[i]])
  )
  
  # loop over repetitions
  for (i in 1:rep) {    
    if(debug){print(paste0('Repetition: ', i, ' out of ', rep))}
    
    # if no data (NIR) is provided generate new curves
    if (data_cond == FALSE) {
      generated_curves <- NIR_curve_generator(n = n_obs, n_harmonics = 30)
      NIR <- as.matrix(generated_curves[, -1])
    }
    
    # transpose data for use with fda package
    data <- t(NIR)
    
    # loop over repetitions
    for (j in 1:length(n_basis)) { 
      
      if(debug){print(paste0('Basis specification: ', j, ' out of ', length(n_basis)))}
      
      tryCatch(
        {
          # calculate responses
          Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
          Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
          Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
          Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
          
          # construct basis for functional representation
          smallbasis <- basis_functions[[j]]
          
          # express NIR data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
          
          # prepare objects for functional regression
          xfdlist <- list(smooth_basis = smooth_basis_fd)
          
          # perform functional regression with cross validation for all 4 scenarios
          f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # save Cross validation criterion in provided container
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress1_1$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress1_2$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress2_1$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress2_2$SSE.CV * 1/(tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", n_basis[j], " b-spline basis functions.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}

### Fourier simulation function
fourier_function <- function(rep, NIR = NULL, n_obs, debug = FALSE) {
  # check if data (NIR) is provided
  if (is.null(NIR)) {
    data_cond <- FALSE
  } else {
    data_cond <- TRUE
  }
  
  # specfiy number of basis functions that should be considered
  n_basis <- seq(from = 5, to = 25, by = 1)
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 6)
  )
  
  colnames(CV_container) <- c(
    "f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline",
    "n_basis", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))
  
  # create basis functions
  basis_functions <- map(
    .x = n_basis,
    .f = function(j) create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j, norder = 4)
  )
  
  # prepare objects for functional linear regression
  betafdPar2_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) fdPar(basis_functions[[i]])
  )
  
  betalist_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) list(smooth_basis_fd_data = betafdPar2_list[[i]])
  )
  
  # loop over repetitions
  for (i in 1:rep) {    
    if(debug){print(paste0('Repetition: ', i, ' out of ', rep))}
    
    # if no data (NIR) is provided generate new curves
    if (data_cond == FALSE) {
      generated_curves <- NIR_curve_generator(n = n_obs, n_harmonics = 30)
      NIR <- as.matrix(generated_curves[, -1])
    }
    
    # transpose data for use with fda package
    data <- t(NIR)
    
    # loop over repetitions
    for (j in 1:length(n_basis)) { 
      
      if(debug){print(paste0('Basis specification: ', j, ' out of ', length(n_basis)))}
      
      tryCatch(
        {
          # calculate responses
          Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
          Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
          Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
          Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
          
          # construct basis for functional representation
          smallbasis <- basis_functions[[j]]
          
          # express NIR data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
          
          # prepare objects for functional regression
          xfdlist <- list(smooth_basis = smooth_basis_fd)
          
          # perform functional regression with cross validation for all 4 scenarios
          f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # save Cross validation criterion in provided container
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress1_1$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress1_2$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress2_1$SSE.CV * 1/(tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr/(tmp_sr + 1) + f_regress2_2$SSE.CV * 1/(tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", n_basis[j], " b-spline basis functions.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}

### FPCA simulation function - bspline basis
fpcr_function <- function(rep, NIR = NULL, n_obs) {
  # check if data (NIR) is provided
  if (is.null(NIR)) {
    data_cond <- FALSE
  } else {
    data_cond <- TRUE
  }

  # set up container for cross validation scores
  CV_container_spline <- c()

  # determine cv-method for later runs
  train.control <- caret::trainControl(method = "cv", number = 10)

  # loop over number of basis functions
  for (j in seq(5, 25, 1)) {
    # create container for cross validation scores for each iteration
    CV_container <- matrix(NaN, nrow = rep, ncol = 5)
    success_count <- 0

    # loop over repetitions
    for (i in 1:rep) {
      print(paste0("Current number of basis functions: ", j, " Success_count: ", success_count))

      # if no data (NIR) is provided generate new curves
      if (data_cond == FALSE) {
        generated_curves <- NIR_curve_generator(n = n_obs, n_harmonics = 30)
        NIR <- as.matrix(generated_curves[, -1])
      }

      tryCatch(
        {
          # calculate responses
          Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
          Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
          Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
          Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)

          # transpose data for use in fda package
          data <- t(NIR)

          # construct basis for functional representation
          smallbasis <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = 4)

          # express NIR data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd

          # perform fPCA
          simulated_pcaObj <- pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)

          # Use results from FPCA to perform linear regression and obtain CV scores
          dataframe1_1 <- data.frame(Y1_1, simulated_pcaObj$scores)
          CV_container[i, 1] <- train(
            Y1_1 ~ .,
            data = dataframe1_1, method = "lm",
            trControl = train.control
          )$results[[2]]

          dataframe1_2 <- data.frame(Y1_2, simulated_pcaObj$scores)
          CV_container[i, 2] <- train(
            Y1_2 ~ .,
            data = dataframe1_2, method = "lm", trControl = train.control
          )$results[[2]]

          dataframe2_1 <- data.frame(Y2_1, simulated_pcaObj$scores)
          CV_container[i, 3] <- train(
            Y2_1 ~ .,
            data = dataframe2_1, method = "lm", trControl = train.control
          )$results[[2]]

          dataframe2_2 <- data.frame(Y2_2, simulated_pcaObj$scores)
          CV_container[i, 4] <- train(
            Y2_2 ~ .,
            data = dataframe2_2, method = "lm", trControl = train.control
          )$results[[2]]

          CV_container[i, 5] <- sum(simulated_pcaObj$varprop)

          success_count <- success_count + 1
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
    # calculate average cv-scores and combine in container with other information
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[6] <- j
    scaled_MSE[7] <- success_count
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
  }
  colnames(CV_container_spline) <- c("f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr", "varprop", "n_basis", "success_count")

  # return object containing averaged cv_scores
  return(CV_container_spline)
}

### FPCA simulation function - Fourier basis
fpcr_fourier_function <- function(rep, NIR = NULL, n_obs) {
  # check if data (NIR) is provided
  if (is.null(NIR)) {
    data_cond <- FALSE
  } else {
    data_cond <- TRUE
  }

  # set up container for cross validation scores
  CV_container_spline <- c()

  # determine cv-method for later runs
  train.control <- caret::trainControl(method = "cv", number = 10)

  # loop over number of basis functions
  for (j in seq(1, 25, 1)) {
    # create container for cross validation scores for each iteration
    CV_container <- matrix(NaN, nrow = rep, ncol = 5)
    success_count <- 0

    # loop over repetitions
    for (i in 1:rep) {
      print(paste0("Current number of basis functions: ", j, " Success_count: ", success_count))

      # if no data (NIR) is provided generate new curves
      if (data_cond == FALSE) {
        generated_curves <- NIR_curve_generator(n = n_obs, n_harmonics = 30)
        NIR <- as.matrix(generated_curves[, -1])
      }

      tryCatch(
        {
          # calculate responses
          Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
          Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
          Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
          Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)

          # transpose data for use in fda package
          data <- t(NIR)

          # construct basis for functional representation
          smallbasis <- create.fourier.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j))

          # express NIR data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd

          # perform fPCA
          simulated_pcaObj <- pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)

          # Use results from FPCA to perform linear regression and obtain CV scores
          dataframe1_1 <- data.frame(Y1_1, simulated_pcaObj$scores)
          CV_container[i, 1] <- train(
            Y1_1 ~ .,
            data = dataframe1_1, method = "lm", trControl = train.control
          )$results[[2]]

          dataframe1_2 <- data.frame(Y1_2, simulated_pcaObj$scores)
          CV_container[i, 2] <- train(
            Y1_2 ~ .,
            data = dataframe1_2, method = "lm", trControl = train.control
          )$results[[2]]

          dataframe2_1 <- data.frame(Y2_1, simulated_pcaObj$scores)
          CV_container[i, 3] <- train(
            Y2_1 ~ .,
            data = dataframe2_1, method = "lm", trControl = train.control
          )$results[[2]]

          dataframe2_2 <- data.frame(Y2_2, simulated_pcaObj$scores)
          CV_container[i, 4] <- train(
            Y2_2 ~ .,
            data = dataframe2_2, method = "lm", trControl = train.control
          )$results[[2]]

          CV_container[i, 5] <- sum(simulated_pcaObj$varprop)
          success_count <- success_count + 1
        },
        error = function(e) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", j, " bspline basis functions for FPCA Run.",
            " Error message: ", cond
          ))
        }
      )
    }
    # calculate average cv-scores and combine in container with other information
    scaled_MSE <- colMeans(CV_container)
    scaled_MSE[6] <- j
    scaled_MSE[7] <- success_count
    CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
  }
  colnames(CV_container_spline) <- c("f1_e1_fpcr_fourier", "f1_e2_fpcr_fourier", "f2_e1_fpcr_fourier", "f2_e2_fpcr_fourier", "varprop", "n_basis", "success_count")

  # return object containing averaged cv_scores
  return(CV_container_spline)
}

##### Define Parallelized versions #####

### ATTENTION! You have to make sure that the cluster workers are appropriately randomized !!!
### Otherwise this whole thing is completely useless

### parallelized b_spline simulation function
bspline_function_par <- function(cl, rep, NIR = NULL, n_obs) {
  # determine available number of cores
  n_cores <- length(cl)

  # determine how the reps are split between cores
  core_reps <- rep(x = rep %/% n_cores, times = n_cores) + c(rep(x = 1, times = rep %% n_cores), rep(x = 0, times = n_cores - (rep %% n_cores)))

  # determine weights for the averaging process
  core_weights <- core_reps / rep

  # perform simulation in parallel
  par_results <- clusterApply(
    cl = cl, x = core_reps,
    fun = function(tmp_rep) {
      bspline_function(rep = tmp_rep, NIR = NIR, n_obs = n_obs, debug = FALSE)
    }
  )

  # Generate object for the results
  n_row <- dim(par_results[[1]])[1]
  n_col <- dim(par_results[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(par_results[[1]])
  colnames(results) <- colnames(par_results[[1]])
  results$n_basis <- par_results[[1]]$n_basis

  # Combine elements of the list par_results to be one result.
  # for cv_scores do weighted averaging
  for (i in 1:4) {
    for (j in 1:n_row) {
      tmp <- unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]][j, i]
        )
      )
      results[j, i] <- sum(tmp * core_weights)
    }
  }
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]]$success_count[j]
        )
      )
    )
  }
  
  # return averaged object
  return(results)
}

### parallelized Fourier simulation function
fourier_function_par <- function(clrep, NIR, n_obs) {
  # determine available number of cores
  n_cores <- length(cl)
  
  # determine how the reps are split between cores
  core_reps <- rep(x = rep %/% n_cores, times = n_cores) + c(rep(x = 1, times = rep %% n_cores), rep(x = 0, times = n_cores - (rep %% n_cores)))
  
  # determine weights for the averaging process
  core_weights <- core_reps / rep
  
  # perform simulation in parallel
  par_results <- clusterApply(
    cl = cl, x = core_reps,
    fun = function(tmp_rep) {
      fourier_function(rep = tmp_rep, NIR = NIR, n_obs = n_obs, debug = FALSE)
    }
  )
  
  # Generate object for the results
  n_row <- dim(par_results[[1]])[1]
  n_col <- dim(par_results[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(par_results[[1]])
  colnames(results) <- colnames(par_results[[1]])
  results$n_basis <- par_results[[1]]$n_basis
  
  # Combine elements of the list par_results to be one result.
  # for cv_scores do weighted averaging
  for (i in 1:4) {
    for (j in 1:n_row) {
      tmp <- unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]][j, i]
        )
      )
      results[j, i] <- sum(tmp * core_weights)
    }
  }
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]]$success_count[j]
        )
      )
    )
  }
  
  # return averaged object
  return(results)
}

### parallelized FPCA simulation function - bspline basis
fpcr_function_par <- function(cl, rep, NIR, n_obs) {
  # determine available number of cores
  n_cores <- length(cl)
  
  # determine how the reps are split between cores
  core_reps <- rep(x = rep %/% n_cores, times = n_cores) + c(rep(x = 1, times = rep %% n_cores), rep(x = 0, times = n_cores - (rep %% n_cores)))
  
  # determine weights for the averaging process
  core_weights <- core_reps / rep
  
  # perform simulation in parallel
  par_results <- clusterApply(
    cl = cl, x = core_reps,
    fun = function(tmp_rep) {
      fpcr_function(rep = tmp_rep, NIR = NIR, n_obs = n_obs)
    }
  )
  
  # Generate object for the results
  n_row <- dim(par_results[[1]])[1]
  n_col <- dim(par_results[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(par_results[[1]])
  colnames(results) <- colnames(par_results[[1]])
  results$n_basis <- par_results[[1]]$n_basis
  
  # Combine elements of the list par_results to be one result.
  # for cv_scores do weighted averaging
  for (i in 1:4) {
    for (j in 1:n_row) {
      tmp <- unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]][j, i]
        )
      )
      results[j, i] <- sum(tmp * core_weights)
    }
  }
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]]$success_count[j]
        )
      )
    )
  }
  
  # return averaged object
  return(results)
}

### parallelized FPCA simulation function - Fourier basis
fpcr_fourier_function_par <- function(cl, rep, NIR, n_obs) {
  # determine available number of cores
  n_cores <- length(cl)
  
  # determine how the reps are split between cores
  core_reps <- rep(x = rep %/% n_cores, times = n_cores) + c(rep(x = 1, times = rep %% n_cores), rep(x = 0, times = n_cores - (rep %% n_cores)))
  
  # determine weights for the averaging process
  core_weights <- core_reps / rep
  
  # perform simulation in parallel
  par_results <- clusterApply(
    cl = cl, x = core_reps,
    fun = function(tmp_rep) {
      fpcr_fourier_function(rep = tmp_rep, NIR = NIR, n_obs = n_obs)
    }
  )
  
  # Generate object for the results
  n_row <- dim(par_results[[1]])[1]
  n_col <- dim(par_results[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(par_results[[1]])
  colnames(results) <- colnames(par_results[[1]])
  results$n_basis <- par_results[[1]]$n_basis
  
  # Combine elements of the list par_results to be one result.
  # for cv_scores do weighted averaging
  for (i in 1:4) {
    for (j in 1:n_row) {
      tmp <- unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]][j, i]
        )
      )
      results[j, i] <- sum(tmp * core_weights)
    }
  }
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:n_cores,
          .f = function(index) par_results[[index]]$success_count[j]
        )
      )
    )
  }
  
  # return averaged object
  return(results)
}

##### Generate Data for Simulations #####
generated_curves <- NIR_curve_generator(n = 1000, n_harmonics = 30)
input_data <- as.matrix(generated_curves[, -1])

##### Perform Simulations #####

test_bspline_function <- bspline_function(rep = 1000, NIR = input_data, n_obs = 1000)
# write.table(test_bspline_function,file="Results/test_bspline_expansion_1000.csv")
test_bspline_function

test_fourier_function <- fourier_function(rep = 1000, NIR = input_data, n_obs = 1000)
# write.table(test_fourier_function,file="Results/test_fourier_expansion_1000.csv")
test_fourier_function

test_fpcr <- fpcr_function(rep = 1000, NIR = input_data, n_obs = 1000)
# write.table(test_fpcr,file="Results/test_fpcr_bsplines_1000.csv")
test_fpcr

test_fpcr2 <- fpcr_fourier_function(rep = 1000, NIR = input_data, n_obs = 1000)
write.table(test_fpcr2, file = "Results/test_fpcr_fourier_1000.csv")
test_fpcr2

test_bspline_function <- bspline_function(rep = 1000, NIR = NIR, n_obs = 60)
# write.table(test_bspline_function,file="Results/test_bspline_expansion_NIR.csv")
test_bspline_function

test_fourier_function <- fourier_function(rep = 1000, NIR = NIR, n_obs = 60)
# write.table(test_fourier_function,file="Results/test_fourier_expansion_NIR.csv")
test_fourier_function

test_fpcr <- fpcr_function(rep = 1000, NIR = NIR, n_obs = 60)
# write.table(test_fpcr,file="Results/test_fpcr_bsplines_NIR.csv")
test_fpcr

test_fpcr2 <- fpcr_fourier_function(rep = 1000, NIR = NIR, n_obs = 60)
# write.table(test_fpcr2, file = "Results/test_fpcr_fourier_NIR.csv")
test_fpcr2

#######################################
###             Test-Area          ###
#######################################
data <- t(NIR)
Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
smallbasis <- create.fourier.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(4))
smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
xfdlist <- list(smooth_basis = smooth_basis_fd)
betabasis1 <- create.constant.basis(c(0, 60))
betafd1 <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
# betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
betafdPar2 <- fdPar(smallbasis)
betalist <- list(smooth_basis_fd_data = betafdPar2)
f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist, betalist)



f_regress1_1original <- fRegress.CVoriginal(y = Y1_1, xfdlist, betalist)
f_regress1_1original$SSE.CV



f_regressn <- fRegress(y = Y1_1, xfdlist, betalist)
