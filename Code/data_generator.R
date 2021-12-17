# load libraries
library(fda)
library(refund)
library(MASS)
library(tidyverse)

# gasoline load data set
data(gasoline)
NIR <- as.matrix(gasoline$NIR)

# get dimensions
n_wl <- dim(NIR)[2]

# define curve-generator function
NIR_curve_generator <- function(n = 1, n_harmonics = 4, n_order = 4,
                                n_knots = 50, plot = FALSE) {

  # define container for output
  curve_matrix <- matrix(data = NA, ncol = n_wl, nrow = n)

  ### create basis object for fpca
  # generate vector of knots
  breaks <- seq(from = 0, to = n_wl, length.out = n_knots)
  # generate bspline basis accordingly
  bspline_basis <- create.bspline.basis(
    rangeval = c(0, n_wl), norder = n_order, breaks = breaks
  )

  # express NIR in terms of the bspline basis
  NIR_bspline_fd <- smooth.basis(y = t(NIR), fdParobj = bspline_basis)$fd

  # perform fpca
  NIR_pcaObj <- pca.fd(
    fdobj = NIR_bspline_fd, nharm = n_harmonics,
    centerfns = TRUE
  )

  # extract eigenvalues
  eigenvalues <- NIR_pcaObj$values[1:n_harmonics]
  var_mat <- diag(eigenvalues)

  # draw realizations of scores
  scores_realizations <- mvrnorm(
    n = n, mu = rep(x = 0, times = n_harmonics), Sigma = var_mat
  )

  # extract coefficients of PCs in bspline-basis representation
  pc_basis_bspline_coef <- NIR_pcaObj$harmonics$coefs

  # get values of bspline basis functions at grid points
  bspline_basis_vals <- eval.basis(evalarg = 1:n_wl, basisobj = bspline_basis)

  # get values of PCs at grid points
  pc_grid_vals <- bspline_basis_vals %*% pc_basis_bspline_coef

  # get mean_function
  mean_function_coefs <- NIR_pcaObj$meanfd$coefs
  mean_function <- as.vector(bspline_basis_vals %*% mean_function_coefs)

  # get new realizations with realizations of scores
  new_realizations <- cbind(
    tibble(observation = paste0("Observation_", 1:n)),
    as_tibble(
      (scores_realizations %*% t(pc_grid_vals)) +
        matrix(data = rep(mean_function, times = n), nrow = n, byrow = TRUE)
    )
  )

  names(new_realizations) <- c("observation", colnames(NIR))

  # if plot == TRUE plot the new observations
  if (plot) {
    new_realizations_plot <- new_realizations %>%
      pivot_longer(cols = !observation, names_to = "wavelength", values_to = "value") %>%
      mutate(wavelength = as.numeric(str_split(wavelength, " ") %>% map_chr(., 1)))

    print(ggplot(data = new_realizations_plot) +
      geom_line(aes(x = wavelength, y = value, col = observation)) +
      theme_light() +
      theme(legend.position = "none"))
  }

  return(new_realizations)
}

##### plot original data #####
plot_original <- FALSE
if (plot_original) {
  NIR_tibble <- cbind(
    tibble(observation = paste0("Observation_", 1:(dim(NIR)[1]))),
    as_tibble(NIR)
  )

  names(NIR_tibble) <- c("observation", colnames(NIR))

  NIR_plot <- NIR_tibble %>%
    pivot_longer(cols = !observation, names_to = "wavelength", values_to = "value") %>%
    mutate(wavelength = as.numeric(str_split(wavelength, " ") %>% map_chr(., 1)))

  ggplot(data = NIR_plot) +
    geom_line(aes(x = wavelength, y = value, col = observation)) +
    theme_light() +
    theme(legend.position = "none")
}
