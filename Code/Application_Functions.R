##### Define Application Functions #####
### b_spline Application function
bspline_appl_function <- function(rep, fold_size, seed, debug = FALSE) {
  set.seed(seed)

  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  octane <- gasoline$octane

  # calculate number of folds
  n_folds <- 60 %/% fold_size

  # set up grid
  grid <- seq(0, 1, length.out = 401)

  # specify number of basis functions that should be considered
  n_basis <- seq(from = 4, to = 15, by = 1)

  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 5)
  )

  colnames(CV_container) <- c(
    "CV_score", "n_basis", "fold_size", "n_folds", "success_count"
  )

  CV_container$n_basis <- n_basis
  CV_container$fold_size <- fold_size
  CV_container$n_folds <- n_folds
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
    .f = function(i) {
      list(
        const = betafdPar1,
        smooth_basis = betafdPar2_list[[i]]
      )
    }
  )

  # loop over repetitions
  for (i in 1:rep) {
    if (debug) {
      print(paste0("Repetition: ", i, " out of ", rep))
    }

    # shuffle data set and transpose for use with fda package
    shuffling <- sample(x = 1:60, size = 60, replace = FALSE)
    data <- t(NIR[shuffling, ])
    responses <- octane[shuffling]

    # loop over diferent number of basis functions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }

      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]

      # create container for values
      MSPE_vec <- vector(mode = "numeric", length = n_folds)

      tryCatch(
        {
          for (m in 1:n_folds) {
            # Choose samples for test data with ordering
            sampling <- shuffling[((fold_size) * (m - 1) + 1):((fold_size) * m)]

            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)
            smooth_basis_fd_train <- smooth.basis(argvals = grid, y = data[, -sampling], fdParobj = smallbasis)$fd
            smooth_basis_fd_test <- smooth.basis(argvals = grid, y = data[, sampling], fdParobj = smallbasis)$fd

            # prepare objects for functional regression
            xfdlist_train <- list(
              const = rep(x = 1, times = 60 - fold_size),
              smooth_basis = smooth_basis_fd_train
            )

            # do fRegress for training data set
            fReg <- fRegress(y = responses[-sampling], xfdlist = xfdlist_train, betalist = betalist_list[[j]])

            # extract constants
            const <- as.numeric(fReg$betaestlist$const$fd$coefs)

            # get estimated coefficients for basis representation of coefficient function
            est_coef <- fReg$betaestlist$smooth_basis$fd$coefs

            # get coefficients for basis representation of observations
            obs_test_coefs <- smooth_basis_fd_test$coefs

            # multiply
            coefficient_matrix <- map(
              .x = 1:fold_size,
              .f = function(i) est_coef %*% obs_test_coefs[, i]
            )

            # integral matrix
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

            # multiply matrices
            regressor_inf <- unlist(
              map(
                .x = 1:fold_size,
                .f = function(i) sum(coefficient_matrix[[i]] * integral_matrix)
              )
            )

            # predict values
            pred <- const + regressor_inf
            # calculate errors
            error <- responses[sampling] - pred
            # calculate MSE
            MSPE_vec[m] <- mean(error^2)
          }

          MSPE_average <- mean(MSPE_vec)

          # generate tmp variable for current number of successful runs
          tmp_sr <- CV_container$success_count[j]

          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$CV_score[j] <- CV_container$CV_score[j] * tmp_sr / (tmp_sr + 1) + MSPE_average / (tmp_sr + 1)

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

### fourier Application function
fourier_appl_function <- function(rep, fold_size, seed, even_basis = FALSE, debug = FALSE) {
  set.seed(seed)
  
  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  octane <- gasoline$octane
  
  # calculate number of folds
  n_folds <- 60 %/% fold_size
  
  # set up grid
  grid <- seq(0, 1, length.out = 401)
  
  # specfiy number of basis functions that should be considered
  if (even_basis == TRUE) {
    n_basis <- seq(from = 1, to = 19, by = 1)
  } else {
    n_basis <- seq(from = 1, to = 19, by = 2)
  }
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 5)
  )
  
  colnames(CV_container) <- c(
    "CV_score", "n_basis", "fold_size", "n_folds", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$fold_size <- fold_size
  CV_container$n_folds <- n_folds
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
    .f = function(i) {
      list(
        const = betafdPar1,
        smooth_basis = betafdPar2_list[[i]]
      )
    }
  )
  
  # loop over repetitions
  for (i in 1:rep) {
    if (debug) {
      print(paste0("Repetition: ", i, " out of ", rep))
    }
    
    # shuffle data set and transpose for use with fda package
    shuffling <- sample(x = 1:60, size = 60, replace = FALSE)
    data <- t(NIR[shuffling, ])
    responses <- octane[shuffling]
    
    # loop over diferent number of basis functions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # create container for values
      MSPE_vec <- vector(mode = "numeric", length = n_folds)
      
      tryCatch(
        {
          for (m in 1:n_folds) {
            # Choose samples for test data with ordering
            sampling <- shuffling[((fold_size) * (m - 1) + 1):((fold_size) * m)]
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)
            smooth_basis_fd_train <- smooth.basis(argvals = grid, y = data[, -sampling], fdParobj = smallbasis)$fd
            smooth_basis_fd_test <- smooth.basis(argvals = grid, y = data[, sampling], fdParobj = smallbasis)$fd
            
            # prepare objects for functional regression
            xfdlist_train <- list(
              const = rep(x = 1, times = 60 - fold_size),
              smooth_basis = smooth_basis_fd_train
            )
            
            # do fRegress for training data set
            fReg <- fRegress(y = responses[-sampling], xfdlist = xfdlist_train, betalist = betalist_list[[j]])
            
            # extract constants
            const <- as.numeric(fReg$betaestlist$const$fd$coefs)
            
            # get estimated coefficients for basis representation of coefficient function
            est_coef <- fReg$betaestlist$smooth_basis$fd$coefs
            
            # get coefficients for basis representation of observations
            obs_test_coefs <- smooth_basis_fd_test$coefs
            
            # multiply
            coefficient_matrix <- map(
              .x = 1:fold_size,
              .f = function(i) est_coef %*% obs_test_coefs[, i]
            )
            
            # integral matrix
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
            
            # multiply matrices
            regressor_inf <- unlist(
              map(
                .x = 1:fold_size,
                .f = function(i) sum(coefficient_matrix[[i]] * integral_matrix)
              )
            )
            
            # predict values
            pred <- const + regressor_inf
            # calculate errors
            error <- responses[sampling] - pred
            # calculate MSE
            MSPE_vec[m] <- mean(error^2)
          }
          
          MSPE_average <- mean(MSPE_vec)
          
          # generate tmp variable for current number of successful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$CV_score[j] <- CV_container$CV_score[j] * tmp_sr / (tmp_sr + 1) + MSPE_average / (tmp_sr + 1)
          
          # count successful runs
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

### monomial Application function
monomial_appl_function <- function(rep, fold_size, seed, debug = FALSE) {
  set.seed(seed)
  
  # load NIR data
  data(gasoline)
  NIR <- as.matrix(gasoline$NIR)
  octane <- gasoline$octane
  
  # calculate number of folds
  n_folds <- 60 %/% fold_size
  
  # set up grid
  grid <- seq(0, 1, length.out = 401)
  
  # specify number of basis functions that should be considered
  n_basis <- seq(from = 1, to = 6, by = 1)
  
  # set up container for averaged cross validation scores
  CV_container <- as.data.frame(
    matrix(data = 0, nrow = length(n_basis), ncol = 5)
  )
  
  colnames(CV_container) <- c(
    "CV_score", "n_basis", "fold_size", "n_folds", "success_count"
  )
  
  CV_container$n_basis <- n_basis
  CV_container$fold_size <- fold_size
  CV_container$n_folds <- n_folds
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
    .f = function(i) {
      list(
        const = betafdPar1,
        smooth_basis = betafdPar2_list[[i]]
      )
    }
  )
  
  # loop over repetitions
  for (i in 1:rep) {
    if (debug) {
      print(paste0("Repetition: ", i, " out of ", rep))
    }
    
    # shuffle data set and transpose for use with fda package
    shuffling <- sample(x = 1:60, size = 60, replace = FALSE)
    data <- t(NIR[shuffling, ])
    responses <- octane[shuffling]
    
    # loop over diferent number of basis functions
    for (j in 1:length(n_basis)) {
      if (debug) {
        print(paste0("Basis specification: ", j, " out of ", length(n_basis)))
      }
      
      # construct basis for functional representation
      smallbasis <- basis_functions[[j]]
      
      # create container for values
      MSPE_vec <- vector(mode = "numeric", length = n_folds)
      
      tryCatch(
        {
          for (m in 1:n_folds) {
            # Choose samples for test data with ordering
            sampling <- shuffling[((fold_size) * (m - 1) + 1):((fold_size) * m)]
            
            # express my_data data in functional basis
            smooth_basis_fd <- smooth.basis(argvals = grid, y = data, fdParobj = smallbasis)
            smooth_basis_fd_train <- smooth.basis(argvals = grid, y = data[, -sampling], fdParobj = smallbasis)$fd
            smooth_basis_fd_test <- smooth.basis(argvals = grid, y = data[, sampling], fdParobj = smallbasis)$fd
            
            # prepare objects for functional regression
            xfdlist_train <- list(
              const = rep(x = 1, times = 60 - fold_size),
              smooth_basis = smooth_basis_fd_train
            )
            
            # do fRegress for training data set
            fReg <- fRegress(y = responses[-sampling], xfdlist = xfdlist_train, betalist = betalist_list[[j]])
            
            # extract constants
            const <- as.numeric(fReg$betaestlist$const$fd$coefs)
            
            # get estimated coefficients for basis representation of coefficient function
            est_coef <- fReg$betaestlist$smooth_basis$fd$coefs
            
            # get coefficients for basis representation of observations
            obs_test_coefs <- smooth_basis_fd_test$coefs
            
            # multiply
            coefficient_matrix <- map(
              .x = 1:fold_size,
              .f = function(i) est_coef %*% obs_test_coefs[, i]
            )
            
            # integral matrix
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
            
            # multiply matrices
            regressor_inf <- unlist(
              map(
                .x = 1:fold_size,
                .f = function(i) sum(coefficient_matrix[[i]] * integral_matrix)
              )
            )
            
            # predict values
            pred <- const + regressor_inf
            # calculate errors
            error <- responses[sampling] - pred
            # calculate MSE
            MSPE_vec[m] <- mean(error^2)
          }
          
          MSPE_average <- mean(MSPE_vec)
          
          # generate tmp variable for current number of successful runs
          tmp_sr <- CV_container$success_count[j]
          
          # Use results from FPCA to perform linear regression and obtain CV scores
          CV_container$CV_score[j] <- CV_container$CV_score[j] * tmp_sr / (tmp_sr + 1) + MSPE_average / (tmp_sr + 1)
          
          # count successful runs
          CV_container$success_count[j] <- CV_container$success_count[j] + 1
        },
        error = function(cond) {
          # in Case an error occurs issue a warning with the corresponding message
          warning(paste0(
            "Problem occured in run ", i, " for ", n_basis[j], " monomial functions.",
            " Error message: ", cond
          ))
        }
      )
    }
  }
  
  # return object containing averaged cv_scores
  return(CV_container)
}

### b_spline fpcr Application function
bspline_fpcr_appl_function <- function(rep, fold_size, nharm, seed, debug = FALSE){
  
}
  
### fourier fpcr Application function
fourier_fpcr_appl_function <- function(rep, fold_size, nharm, seed, even_basis = FALSE, debug = FALSE){
  
}
    
### monomial fpcr Application function
monomial_fpcr_appl_function <- function(rep, fold_size, nharm, seed, debug = FALSE){
  
}