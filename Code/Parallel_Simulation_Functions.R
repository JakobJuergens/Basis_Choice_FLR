##### Define Parallelized versions #####

### parallelized b_spline simulation function
bspline_function_stitch <- function(path) {

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
fourier_function_stitch <- function(rep, my_data, n_obs) {

  
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
fpcr_function_stitch <- function(cl, rep, my_data, n_obs) {
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
      fpcr_function(rep = tmp_rep, my_data = my_data, n_obs = n_obs, debug = FALSE)
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
  for (i in 1:5) {
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
fpcr_fourier_function_stitch <- function(cl, rep, my_data, n_obs) {
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
      fpcr_fourier_function(rep = tmp_rep, my_data = my_data, n_obs = n_obs, debug = FALSE)
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
  for (i in 1:5) {
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