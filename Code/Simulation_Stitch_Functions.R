##### Functions for stitching together the results of different simulation runs #####

### b_spline simulation function
simulation_stitch <- function(path) {

  # get files in path
  path_files <- list.files(path = path)

  # read into list
  output_dfs <- map(
    .x = path_files,
    .f = function(obj) readRDS(paste0(path, "/", obj))
  )

  # Generate object for the results
  n_row <- dim(output_dfs[[1]])[1]
  n_col <- dim(output_dfs[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(output_dfs[[1]])
  colnames(results) <- colnames(output_dfs[[1]])
  results$n_basis <- output_dfs[[1]]$n_basis

  # Combine elements of the list par_results to be one result.
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:length(output_dfs),
          .f = function(index) output_dfs[[index]]$success_count[j]
        )
      )
    )
  }

  # for cv_scores do weighted averaging
  for (i in 1:4) {
    for (j in 1:n_row) {
      results[j, i] <- sum(
        unlist(
          map(
            .x = 1:length(output_dfs),
            .f = function(index) {
              output_dfs[[index]][j, i] * ifelse(results$success_count[j] != 0, (output_dfs[[index]]$success_count[j] / results$success_count[j]), 0)
            }
          )
        ),
        na.rm = TRUE
      )
    }
  }

  # return averaged object
  return(results)
}

### FPCA simulation function - bspline basis
fpcr_simulation_stitch <- function(path) {

  # get files in path
  path_files <- list.files(path = path)

  # read into list
  output_dfs <- map(
    .x = path_files,
    .f = function(obj) readRDS(paste0(path, "/", obj))
  )

  # Generate object for the results
  n_row <- dim(output_dfs[[1]])[1]
  n_col <- dim(output_dfs[[1]])[2]
  results <- as.data.frame(
    matrix(data = NA, nrow = n_row, ncol = n_col)
  )
  rownames(results) <- rownames(output_dfs[[1]])
  colnames(results) <- colnames(output_dfs[[1]])
  results$n_basis <- output_dfs[[1]]$n_basis

  # Combine elements of the list par_results to be one result.
  # for success_count do summation
  for (j in 1:n_row) {
    results$success_count[j] <- sum(
      unlist(
        map(
          .x = 1:length(output_dfs),
          .f = function(index) output_dfs[[index]]$success_count[j]
        )
      )
    )
  }

  # for cv_scores do weighted averaging
  for (i in 1:5) {
    for (j in 1:n_row) {
      results[j, i] <- sum(
        unlist(
          map(
            .x = 1:length(output_dfs),
            .f = function(index) {
              output_dfs[[index]][j, i] * ifelse(results$success_count[j] != 0, (output_dfs[[index]]$success_count[j] / results$success_count[j]), 0)
            }
          )
        ),
        na.rm = TRUE
      )
    }
  }

  # return averaged object
  return(results)
}
