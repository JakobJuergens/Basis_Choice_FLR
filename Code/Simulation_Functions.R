##### Define Simulation Functions #####
### b_spline simulation function
bspline_function <- function(rep, my_data = NULL, n_obs, seed, debug = FALSE) {
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
  data_cond <- TRUE
  if (is.null(my_data)) {
    data_cond <- FALSE
  }
  
  # specify number of basis functions that should be considered
  n_basis <- seq(from = 4, to = 25, by = 1)
  
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
    
    # loop over diferent number of basis functions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      tryCatch(
      {
      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # express my_data data in functional basis
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
      CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress1_1$SSE.CV * 1 / (tmp_sr + 1)
      CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress1_2$SSE.CV * 1 / (tmp_sr + 1)
      CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress2_1$SSE.CV * 1 / (tmp_sr + 1)
      CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress2_2$SSE.CV * 1 / (tmp_sr + 1)
      
      # count successful runs
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
fourier_function <- function(rep, my_data = NULL, n_obs, seed, even_basis = FALSE, debug = FALSE) {
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
  if(even_basis == TRUE){
    n_basis <- seq(from = 1, to = 25, by = 1)
  } else{
    n_basis <- seq(from = 1, to = 25, by = 2)
  }
  
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
  if(even_basis == TRUE){
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) {
        if(j %% 2 == 1){
          return(create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j))
        } else{
          return(create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j, dropind = j))
        }
      }
    )
  } else{
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j)
    )
  }

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
    
    # loop over repetitions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      tryCatch(
        {
          # construct basis for functional representation
          smallbasis <- basis_functions[[j]]
          
          # express my_data data in functional basis
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
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress1_1$SSE.CV * 1 / (tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress1_2$SSE.CV * 1 / (tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress2_1$SSE.CV * 1 / (tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr / (tmp_sr + 1) + f_regress2_2$SSE.CV * 1 / (tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", n_basis[j], " fourier basis functions.",
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
    .f = function(j) create.bspline.basis(rangeval = c(0, length(grid)), nbasis = j, norder = 4)
  )
  
  # determine cv-method for later runs
  train.control <- caret::trainControl(method = "cv", number = 10)
  
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
    
    # loop over repetitions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      tryCatch(
        {
          # construct basis for functional representation
          smallbasis <- basis_functions[[j]]
          
          # express my_data data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
          
          # perform fPCA
          simulated_pcaObj <- pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)
          
          # combine objects into data frames for linear regression
          dataframe1_1 <- data.frame(Y1_1, simulated_pcaObj$scores)
          dataframe1_2 <- data.frame(Y1_2, simulated_pcaObj$scores)
          dataframe2_1 <- data.frame(Y2_1, simulated_pcaObj$scores)
          dataframe2_2 <- data.frame(Y2_2, simulated_pcaObj$scores)
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y1_1 ~ .,
            data = dataframe1_1, method = "lm",
            trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y1_2 ~ .,
            data = dataframe1_2, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y2_1 ~ .,
            data = dataframe2_1, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y2_2 ~ .,
            data = dataframe2_2, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          # update varprop
          CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(simulated_pcaObj$varprop) * 1 / (tmp_sr + 1)
          
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

### FPCA simulation function - Fourier basis
fpcr_fourier_function <- function(rep, my_data = NULL, n_obs, nharm, seed, even_basis = FALSE, debug = FALSE) {
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
  if(even_basis == TRUE){
    n_basis <- seq(from = 1, to = 25, by = 1)
  } else{
    n_basis <- seq(from = 1, to = 25, by = 2)
  }
  
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
  if(even_basis == TRUE){
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) {
        if(j %% 2 == 1){
          return(create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j))
        } else{
          return(create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j, dropind = j))
        }
      }
    )
  } else{
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) create.fourier.basis(rangeval = c(0, length(grid)), nbasis = j)
    )
  }
  
  # determine cv-method for later runs
  train.control <- caret::trainControl(method = "cv", number = 10)
  
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
    
    # loop over repetitions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      tryCatch(
        {
          # construct basis for functional representation
          smallbasis <- basis_functions[[j]]
          
          # express my_data data in functional basis
          smooth_basis_fd <- smooth.basis(y = data, fdParobj = smallbasis)$fd
          
          # perform fPCA
          simulated_pcaObj <- pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)
          
          # combine objects into data frames for linear regression
          dataframe1_1 <- data.frame(Y1_1, simulated_pcaObj$scores)
          dataframe1_2 <- data.frame(Y1_2, simulated_pcaObj$scores)
          dataframe2_1 <- data.frame(Y2_1, simulated_pcaObj$scores)
          dataframe2_2 <- data.frame(Y2_2, simulated_pcaObj$scores)
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y1_1 ~ .,
            data = dataframe1_1, method = "lm",
            trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y1_2 ~ .,
            data = dataframe1_2, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y2_1 ~ .,
            data = dataframe2_1, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + train(
            Y2_2 ~ .,
            data = dataframe2_2, method = "lm", trControl = train.control
          )$results[[2]] * 1 / (tmp_sr + 1)
          
          # update varprop
          CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(simulated_pcaObj$varprop) * 1 / (tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", j, " fourier basis functions for FPCA Run.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}
