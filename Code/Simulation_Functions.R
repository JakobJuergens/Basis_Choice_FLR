##### Define Simulation Functions #####
### b_spline simulation function
bspline_function <- function(rep, my_data = NULL, n_obs, seed, debug = FALSE) {
  set.seed(seed)

  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  
  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
  
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
    .f = function(j) create.bspline.basis(rangeval = c(0, 1), nbasis = j, norder = 4)
  )

  # prepare objects for functional linear regression
  betafdPar1 <- fdPar(fd(0, create.constant.basis(c(0, 1))))
  
  betafdPar2_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) fdPar(basis_functions[[i]])
  )

  betalist_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) list(const = betafdPar1, 
                          smooth_basis = betafdPar2_list[[i]]
                          )
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
          smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd

          # prepare objects for functional regression
          xfdlist <- list(const = rep(x = 1, times = n_obs), 
                          smooth_basis = smooth_basis_fd)

          # perform functional regression with cross validation for all 4 scenarios
          f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist = xfdlist, betalist = betalist_list[[j]])

          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]

          # save Cross validation criterion in provided container
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_2$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_2$SSE.CV) * 1 / (tmp_sr + 1)

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

  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)

  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
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
  if (even_basis == TRUE) {
    n_basis <- seq(from = 1, to = 25, by = 1)
  } else {
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
  if (even_basis == TRUE) {
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) {
        if (j %% 2 == 1) {
          return(create.fourier.basis(rangeval = c(0, 1), nbasis = j))
        } else {
          return(create.fourier.basis(rangeval = c(0, 1), nbasis = j, dropind = j))
        }
      }
    )
  } else {
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) create.fourier.basis(rangeval = c(0, 1), nbasis = j)
    )
  }
  
  # prepare objects for functional linear regression
  betafdPar1 <- fdPar(fd(0, create.constant.basis(c(0, 1))))
  
  betafdPar2_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) fdPar(basis_functions[[i]])
  )
  
  betalist_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) list(const = betafdPar1, 
                          smooth_basis = betafdPar2_list[[i]]
    )
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
          smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd

          # prepare objects for functional regression
          xfdlist <- list(const = rep(x = 1, times = n_obs), 
                          smooth_basis = smooth_basis_fd)

          # perform functional regression with cross validation for all 4 scenarios
          f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]

          # save Cross validation criterion in provided container
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_2$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_2$SSE.CV) * 1 / (tmp_sr + 1)

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

### monomial simulation function
monomial_function <- function(rep, my_data = NULL, n_obs, seed, debug = FALSE) {
  set.seed(seed)
  
  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  
  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
  
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
  n_basis <- seq(from = 1, to = 6, by = 1)
  
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
    .f = function(j) create.monomial.basis(rangeval = c(0,1), nbasis = j)
  )
  
  # prepare objects for functional linear regression
  betafdPar1 <- fdPar(fd(0, create.constant.basis(c(0, 1))))
  
  betafdPar2_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) fdPar(basis_functions[[i]])
  )
  
  betalist_list <- map(
    .x = 1:length(n_basis),
    .f = function(i) list(const = betafdPar1, 
                          smooth_basis = betafdPar2_list[[i]]
    )
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
          smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd
          
          # prepare objects for functional regression
          xfdlist <- list(const = rep(x = 1, times = n_obs), 
                          smooth_basis = smooth_basis_fd)
          
          # perform functional regression with cross validation for all 4 scenarios
          f_regress1_1 <- fRegress.CVk(y = Y1_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress1_2 <- fRegress.CVk(y = Y1_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_1 <- fRegress.CVk(y = Y2_1, xfdlist = xfdlist, betalist = betalist_list[[j]])
          f_regress2_2 <- fRegress.CVk(y = Y2_2, xfdlist = xfdlist, betalist = betalist_list[[j]])
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # save Cross validation criterion in provided container
          CV_container$f1_e1_spline[j] <- CV_container$f1_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f1_e2_spline[j] <- CV_container$f1_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress1_2$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e1_spline[j] <- CV_container$f2_e1_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_1$SSE.CV) * 1 / (tmp_sr + 1)
          CV_container$f2_e2_spline[j] <- CV_container$f2_e2_spline[j] * tmp_sr / (tmp_sr + 1) + sqrt(f_regress2_2$SSE.CV) * 1 / (tmp_sr + 1)
          
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

### FPCA simulation function - bspline basis
fpcr_function <- function(rep, my_data = NULL, n_obs, nharm, seed, debug = FALSE) {
  set.seed(seed)

  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  
  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
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
    matrix(data = 0, nrow = length(n_basis), ncol = 6)
  )

  colnames(CV_container) <- c(
    "f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr",
    # "varprop", 
    "n_basis", "success_count"
  )

  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))

  # create basis functions
  basis_functions <- map(
    .x = n_basis,
    .f = function(j) create.bspline.basis(rangeval = c(0, 1), nbasis = j, norder = 4)
  )

  # find number of elements in fold
  n_elem_fold <- floor(n_obs / 10)
  
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

    # loop over different basis specifications
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }

      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # create container for values
      MSPE_mat <- as.data.frame(matrix(0, nrow = 10, ncol = 4))
      colnames(MSPE_mat) <- c('MSPE1_1', 'MSPE1_2', 'MSPE2_1', 'MSPE2_2')
      
      tryCatch(
        {
          
          #Randomly reorder samples
          shuffled <- 1:n_obs # sample(x = 1:n_obs, size = n_obs, replace = FALSE)
          
          for(m in 1:10){
            #Choose samples for test data with ordering
            sampling <- sort(x = shuffled[((n_elem_fold)*(m-1)+1):((n_elem_fold)*m)])
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd
            smooth_basis_fd_train <- smooth_basis_fd[-sampling]
            smooth_basis_fd_test  <- smooth_basis_fd[sampling]
            
            # perform fPCA on training data
            train_fd <- pca.fd(fdobj = smooth_basis_fd_train, 
                               nharm = nharm, centerfns = TRUE)
            
            # extract scores for training set
            train_scores <- train_fd$scores
            
            # create container for estimated scores in test data
            scores_mat  <- matrix(0, nrow = n_elem_fold, ncol = nharm)
            
            # grid for approximation of intergrals on [0,1]
            tmp_grid <- seq(0, 1, length.out = 1000)
            # get values of basis functions at approximation grid
            basis_eval <- eval.basis(evalarg = tmp_grid, basisobj = smallbasis)
            # create a container for the integrals over the pairwise products 
            # of the basis functions
            integral_matrix <- matrix(data = 0, nrow = n_basis[j], ncol = n_basis[j])
            
            # calculate integrals over pairwise products of basis functions
            # and save in provided matrix
            for (l in 1:n_basis[j]) {
              for (k in 1:n_basis[j]) {
                integral_matrix[l, k] <- integrate(function(t) {
                  approx(
                    x = tmp_grid,
                    y = basis_eval[, l] * t(basis_eval[, k]), xout = t
                  )$y
                }, lower = 0, upper = 1)$value
              }
            }
            
            # iterate over curves in the test set
            for (p in 1:n_elem_fold){
              # iterate over fpcs
              for(k in 1 : nharm){
                # calculate weight matrix for the integrated pairwise products of 
                # basis functions (integral_matrix)
                newfun <- (smooth_basis_fd_test[p]$coefs - train_fd$meanfd$coefs) %*% t(train_fd$harmonics[k]$coefs)
                # do element wise product between matrices
                weighted_basis_integrals_products <- newfun * integral_matrix
                # sum over all elements of the matrix to obtain the score estimate
                scores_mat[p, k] <- sum(weighted_basis_integrals_products)
              }
            }
            
            # transform estimated scores to tibble for prediction
            scores_tib <- as_tibble(scores_mat)
            names(scores_tib) <- paste0('harm_', 1:nharm)
            
            # combine objects into data frames for linear regression
            dataframe1_1 <- as_tibble(cbind(Y1_1[-sampling], train_fd$scores))
            dataframe1_2 <- as_tibble(cbind(Y1_2[-sampling], train_fd$scores))
            dataframe2_1 <- as_tibble(cbind(Y2_1[-sampling], train_fd$scores))
            dataframe2_2 <- as_tibble(cbind(Y2_2[-sampling], train_fd$scores))
            
            # give correct names
            names(dataframe1_1) <- c('Y1_1', paste0('harm_', 1:nharm))
            names(dataframe1_2) <- c('Y1_2', paste0('harm_', 1:nharm))
            names(dataframe2_1) <- c('Y2_1', paste0('harm_', 1:nharm))
            names(dataframe2_2) <- c('Y2_2', paste0('harm_', 1:nharm))
            
            # estimate linear models
            model1_1 <- lm(Y1_1 ~., dataframe1_1) 
            model1_2 <- lm(Y1_2 ~., dataframe1_2)
            model2_1 <- lm(Y2_1 ~., dataframe2_1)
            model2_2 <- lm(Y2_2 ~., dataframe2_2)
            
            # do prediction using estimated scores
            test_prediction1_1 <- unname(predict(object = model1_1, newdata = scores_tib))
            test_prediction1_2 <- unname(predict(object = model1_2, newdata = scores_tib))
            test_prediction2_1 <- unname(predict(object = model2_1, newdata = scores_tib))
            test_prediction2_2 <- unname(predict(object = model2_2, newdata = scores_tib))
            
            # calculate RMSE from predicted values
            MSPE_mat[m, 1] <- sqrt(1/10 * sum((test_prediction1_1 - Y1_1)^2))
            MSPE_mat[m, 2] <- sqrt(1/10 * sum((test_prediction1_2 - Y1_2)^2))
            MSPE_mat[m, 3] <- sqrt(1/10 * sum((test_prediction2_1 - Y2_1)^2))
            MSPE_mat[m, 4] <- sqrt(1/10 * sum((test_prediction2_2 - Y2_2)^2))
          }
          
          MSPE_average <- colMeans(MSPE_mat)

          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]

          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[1] * 1 / (tmp_sr + 1)
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[2] * 1 / (tmp_sr + 1)
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[3] * 1 / (tmp_sr + 1)
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[4] * 1 / (tmp_sr + 1)
          
          # update varprop
          # CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(train_fd$varprop) * 1 / (tmp_sr + 1)
          
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

  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  
  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
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
  if (even_basis == TRUE) {
    n_basis <- seq(from = 1, to = 25, by = 1)
  } else {
    n_basis <- seq(from = 1, to = 25, by = 2)
  }

  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 6)
  )

  colnames(CV_container) <- c(
    "f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr",
    # "varprop", 
    "n_basis", "success_count"
  )

  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))

  # create basis functions
  if (even_basis == TRUE) {
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) {
        if (j %% 2 == 1) {
          return(create.fourier.basis(rangeval = c(0, 1), nbasis = j))
        } else {
          return(create.fourier.basis(rangeval = c(0, 1), nbasis = j, dropind = j))
        }
      }
    )
  } else {
    basis_functions <- map(
      .x = n_basis,
      .f = function(j) create.fourier.basis(rangeval = c(0, 1), nbasis = j)
    )
  }

  # find number of elements in fold
  n_elem_fold <- floor(n_obs / 10)
  
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
    
    # loop over different basis specifications
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # create container for values
      MSPE_mat <- as.data.frame(matrix(0, nrow = 10, ncol = 4))
      colnames(MSPE_mat) <- c('MSPE1_1', 'MSPE1_2', 'MSPE2_1', 'MSPE2_2')
      
      tryCatch(
        {
          
          #Randomly reorder samples
          shuffled <- 1:n_obs # sample(x = 1:n_obs, size = n_obs, replace = FALSE)
          
          for(m in 1:10){
            #Choose samples for test data with ordering
            sampling <- sort(x = shuffled[((n_elem_fold)*(m-1)+1):((n_elem_fold)*m)])
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd
            smooth_basis_fd_train <- smooth_basis_fd[-sampling]
            smooth_basis_fd_test  <- smooth_basis_fd[sampling]
            
            # perform fPCA on training data
            train_fd <- pca.fd(fdobj = smooth_basis_fd_train, 
                               nharm = nharm, centerfns = TRUE)
            
            # extract scores for training set
            train_scores <- train_fd$scores
            
            # create container for estimated scores in test data
            scores_mat  <- matrix(0, nrow = n_elem_fold, ncol = nharm)
            
            # grid for approximation of intergrals on [0,1]
            tmp_grid <- seq(0, 1, length.out = 1000)
            # get values of basis functions at approximation grid
            basis_eval <- eval.basis(evalarg = tmp_grid, basisobj = smallbasis)
            # create a container for the integrals over the pairwise products 
            # of the basis functions
            integral_matrix <- matrix(data = 0, nrow = n_basis[j], ncol = n_basis[j])
            
            # calculate integrals over pairwise products of basis functions
            # and save in provided matrix
            for (l in 1:n_basis[j]) {
              for (k in 1:n_basis[j]) {
                integral_matrix[l, k] <- integrate(function(t) {
                  approx(
                    x = tmp_grid,
                    y = basis_eval[, l] * t(basis_eval[, k]), xout = t
                  )$y
                }, lower = 0, upper = 1)$value
              }
            }
            
            # iterate over curves in the test set
            for (p in 1:n_elem_fold){
              # iterate over fpcs
              for(k in 1 : nharm){
                # calculate weight matrix for the integrated pairwise products of 
                # basis functions (integral_matrix)
                newfun <- (smooth_basis_fd_test[p]$coefs - train_fd$meanfd$coefs) %*% t(train_fd$harmonics[k]$coefs)
                # do element wise product between matrices
                weighted_basis_integrals_products <- newfun * integral_matrix
                # sum over all elements of the amtrix to obtain the score estimate
                scores_mat[p, k] <- sum(weighted_basis_integrals_products)
              }
            }
            
            # transform estimated scores to tibble for prediction
            scores_tib <- as_tibble(scores_mat)
            names(scores_tib) <- paste0('harm_', 1:nharm)
            
            # combine objects into data frames for linear regression
            dataframe1_1 <- as_tibble(cbind(Y1_1[-sampling], train_fd$scores))
            dataframe1_2 <- as_tibble(cbind(Y1_2[-sampling], train_fd$scores))
            dataframe2_1 <- as_tibble(cbind(Y2_1[-sampling], train_fd$scores))
            dataframe2_2 <- as_tibble(cbind(Y2_2[-sampling], train_fd$scores))
            
            # give correct names
            names(dataframe1_1) <- c('Y1_1', paste0('harm_', 1:nharm))
            names(dataframe1_2) <- c('Y1_2', paste0('harm_', 1:nharm))
            names(dataframe2_1) <- c('Y2_1', paste0('harm_', 1:nharm))
            names(dataframe2_2) <- c('Y2_2', paste0('harm_', 1:nharm))
            
            # estimate linear models
            model1_1 <- lm(Y1_1 ~., dataframe1_1) 
            model1_2 <- lm(Y1_2 ~., dataframe1_2)
            model2_1 <- lm(Y2_1 ~., dataframe2_1)
            model2_2 <- lm(Y2_2 ~., dataframe2_2)
            
            # do prediction using estimated scores
            test_prediction1_1 <- unname(predict(object = model1_1, newdata = scores_tib))
            test_prediction1_2 <- unname(predict(object = model1_2, newdata = scores_tib))
            test_prediction2_1 <- unname(predict(object = model2_1, newdata = scores_tib))
            test_prediction2_2 <- unname(predict(object = model2_2, newdata = scores_tib))
            
            # calculate RMSE from predicted values
            MSPE_mat[m, 1] <- sqrt(1/10 * sum((test_prediction1_1 - Y1_1)^2))
            MSPE_mat[m, 2] <- sqrt(1/10 * sum((test_prediction1_2 - Y1_2)^2))
            MSPE_mat[m, 3] <- sqrt(1/10 * sum((test_prediction2_1 - Y2_1)^2))
            MSPE_mat[m, 4] <- sqrt(1/10 * sum((test_prediction2_2 - Y2_2)^2))
          }
          
          MSPE_average <- colMeans(MSPE_mat)
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[1] * 1 / (tmp_sr + 1)
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[2] * 1 / (tmp_sr + 1)
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[3] * 1 / (tmp_sr + 1)
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[4] * 1 / (tmp_sr + 1)
          
          # update varprop
          # CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(train_fd$varprop) * 1 / (tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", n_basis[j], " fourier basis functions for FPCA Run.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}

### FPCA simulation function - monomial basis
fpcr_monomial_function <- function(rep, my_data = NULL, n_obs, nharm, seed, debug = FALSE) {
  set.seed(seed)
  
  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  
  ### set up coefficient "functions" / error terms
  grid <- seq(0, 1, length.out = 401)
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
  n_basis <- seq(from = 4, to = 12, by = 1)
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 6)
  )
  
  colnames(CV_container) <- c(
    "f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr",
    # "varprop", 
    "n_basis", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$success_count <- rep(x = 0, times = length(n_basis))
  
  # create basis functions
  basis_functions <- map(
    .x = n_basis,
    .f = function(j) create.monomial.basis(rangeval = c(0, 1), nbasis = j)
  )
  
  # find number of elements in fold
  n_elem_fold <- floor(n_obs / 10)
  
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
    
    # loop over different basis specifications
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # create container for values
      MSPE_mat <- as.data.frame(matrix(0, nrow = 10, ncol = 4))
      colnames(MSPE_mat) <- c('MSPE1_1', 'MSPE1_2', 'MSPE2_1', 'MSPE2_2')
      
      tryCatch(
        {
          
          #Randomly reorder samples
          shuffled <- 1:n_obs # sample(x = 1:n_obs, size = n_obs, replace = FALSE)
          
          for(m in 1:10){
            #Choose samples for test data with ordering
            sampling <- sort(x = shuffled[((n_elem_fold)*(m-1)+1):((n_elem_fold)*m)])
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)$fd
            smooth_basis_fd_train <- smooth_basis_fd[-sampling]
            smooth_basis_fd_test  <- smooth_basis_fd[sampling]
            
            # perform fPCA on training data
            train_fd <- pca.fd(fdobj = smooth_basis_fd_train, 
                               nharm = nharm, centerfns = TRUE)
            
            # extract scores for training set
            train_scores <- train_fd$scores
            
            # create container for estimated scores in test data
            scores_mat  <- matrix(0, nrow = n_elem_fold, ncol = nharm)
            
            # grid for approximation of intergrals on [0,1]
            tmp_grid <- seq(0, 1, length.out = 1000)
            # get values of basis functions at approximation grid
            basis_eval <- eval.basis(evalarg = tmp_grid, basisobj = smallbasis)
            # create a container for the integrals over the pairwise products 
            # of the basis functions
            integral_matrix <- matrix(data = 0, nrow = n_basis[j], ncol = n_basis[j])
            
            # calculate integrals over pairwise products of basis functions
            # and save in provided matrix
            for (l in 1:length(n_basis[j])) {
              for (k in 1:length(n_basis[j])) {
                integral_matrix[l, k] <- integrate(function(t) {
                  approx(
                    x = tmp_grid,
                    y = basis_eval[, l] * t(basis_eval[, k]), xout = t
                  )$y
                }, lower = 0, upper = 1)$value
              }
            }
            
            # iterate over curves in the test set
            for (p in 1:n_elem_fold){
              # iterate over fpcs
              for(k in 1 : nharm){
                # calculate weight matrix for the integrated pairwise products of 
                # basis functions (integral_matrix)
                newfun <- (smooth_basis_fd_test[p]$coefs - train_fd$meanfd$coefs) %*% t(train_fd$harmonics[k]$coefs)
                # do element wise product between matrices
                weighted_basis_integrals_products <- newfun * integral_matrix
                # sum over all elements of the matrix to obtain the score estimate
                scores_mat[p, k] <- sum(weighted_basis_integrals_products)
              }
            }
            
            # transform estimated scores to tibble for prediction
            scores_tib <- as_tibble(scores_mat)
            names(scores_tib) <- paste0('harm_', 1:nharm)
            
            # combine objects into data frames for linear regression
            dataframe1_1 <- as_tibble(cbind(Y1_1[-sampling], train_fd$scores))
            dataframe1_2 <- as_tibble(cbind(Y1_2[-sampling], train_fd$scores))
            dataframe2_1 <- as_tibble(cbind(Y2_1[-sampling], train_fd$scores))
            dataframe2_2 <- as_tibble(cbind(Y2_2[-sampling], train_fd$scores))
            
            # give correct names
            names(dataframe1_1) <- c('Y1_1', paste0('harm_', 1:nharm))
            names(dataframe1_2) <- c('Y1_2', paste0('harm_', 1:nharm))
            names(dataframe2_1) <- c('Y2_1', paste0('harm_', 1:nharm))
            names(dataframe2_2) <- c('Y2_2', paste0('harm_', 1:nharm))
            
            # estimate linear models
            model1_1 <- lm(Y1_1 ~., dataframe1_1) 
            model1_2 <- lm(Y1_2 ~., dataframe1_2)
            model2_1 <- lm(Y2_1 ~., dataframe2_1)
            model2_2 <- lm(Y2_2 ~., dataframe2_2)
            
            # do prediction using estimated scores
            test_prediction1_1 <- unname(predict(object = model1_1, newdata = scores_tib))
            test_prediction1_2 <- unname(predict(object = model1_2, newdata = scores_tib))
            test_prediction2_1 <- unname(predict(object = model2_1, newdata = scores_tib))
            test_prediction2_2 <- unname(predict(object = model2_2, newdata = scores_tib))
            
            # calculate RMSE from predicted values
            MSPE_mat[m, 1] <- sqrt(1/10 * sum((test_prediction1_1 - Y1_1)^2))
            MSPE_mat[m, 2] <- sqrt(1/10 * sum((test_prediction1_2 - Y1_2)^2))
            MSPE_mat[m, 3] <- sqrt(1/10 * sum((test_prediction2_1 - Y2_1)^2))
            MSPE_mat[m, 4] <- sqrt(1/10 * sum((test_prediction2_2 - Y2_2)^2))
          }
          
          MSPE_average <- colMeans(MSPE_mat)
          
          # generate tmp variable for current number of succesful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$f1_e1_fpcr[j] <- CV_container$f1_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[1] * 1 / (tmp_sr + 1)
          CV_container$f1_e2_fpcr[j] <- CV_container$f1_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[2] * 1 / (tmp_sr + 1)
          CV_container$f2_e1_fpcr[j] <- CV_container$f2_e1_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[3] * 1 / (tmp_sr + 1)
          CV_container$f2_e2_fpcr[j] <- CV_container$f2_e2_fpcr[j] * tmp_sr / (tmp_sr + 1) + MSPE_average[4] * 1 / (tmp_sr + 1)
          
          # update varprop
          # CV_container$varprop[j] <- CV_container$varprop[j] * tmp_sr / (tmp_sr + 1) + sum(train_fd$varprop) * 1 / (tmp_sr + 1)
          
          # count succesfull runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", j, " monomial basis functions for FPCA Run.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}