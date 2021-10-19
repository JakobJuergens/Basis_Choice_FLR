library(tidyverse)
library(parallel)

# set seed for overall script (this fixes the seeds for the underlying parallel sessions)
set.seed(42)

##### example for approximating pi #####

# This helper function is what is later called in parallel
parallel_helper <- function(seed = 42, n_iter = 1){
  # due to handling of random states in cluster, this has to happen in the function itself
  set.seed(seed)
  
  # generate random coordinates on [-1, 1]x[-1, 1] (bivariate uniform distribution)
  points <- tibble(x = runif(n = n_iter, min = -1, max = 1),
                   y = runif(n = n_iter, min = -1, max = 1)) %>% 
  # check if points are in cricle around (0,0) with radius 1
    mutate(in_circle = x^2 + y^2 <= 1)
  
  # number of points in circle
  n_in_circle <- sum(points$in_circle)
  
  # return
  return(n_in_circle = n_in_circle)
}

# this function actually performs the simulation by calling the helper function
# in parallel and using its outputs
simulation_function <- function(cl, n_iter){
  
  # get number of cores in cluster
  n_cores <- length(cl)
  
  # seeds for different cores
  seeds <- sample(x = 1:10e6, size = n_cores)
  
  # fix distribution of tasks for each core
  iters <- rep(n_iter %/% n_cores, times = n_cores) + 
    c(rep(1, times = n_iter %% n_cores), rep(0, times = n_cores - n_iter %% n_cores))
  
  # perform parallel runs of function parallel_helper
  par_results <- unlist(
    clusterApply(cl = cl,
                 x = 1:n_cores,
                 fun = function(i) 
                   parallel_helper(seed = seeds[i], n_iter = iters[i])
                 )
    )
  
  # Aggregate results as desired
  pi_approx <- 4 * sum(par_results) / n_iter
  
  # return approximation of pi
  return(pi_approx)
}

# detects number of cores we can use for parallelization
n_cores <- detectCores()

# make cluster (this one should also work on Windows)
cl <- makePSOCKcluster(n_cores)

# export function parallel_helper to cluster sessions
clusterExport(cl = cl, varlist = list('parallel_helper'), envir = .GlobalEnv)

# load tidyverse in cluster sessions
clusterEvalQ(cl = cl, library(tidyverse))

# use function to perform MC-Simulation  
pi_approx <- simulation_function(cl = cl, n_iter = 100000)

# stop cluster
stopCluster(cl)

