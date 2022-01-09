##### Functions for stitching together the results of different simulation runs #####

### b_spline simulation function
simulation_stitch <- function(path) {

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

### FPCA simulation function - bspline basis
fpcr_simulation_stitch <- function(path) {

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