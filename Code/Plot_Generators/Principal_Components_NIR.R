##### This script generates a plot showing of a chosen number of principal components #
# of the NIR-spectrum data from the gasoline data set and representations of the
# original curves approximated by different numbers of these curves #####

# Clear Workspace
rm(list = ls())

# load libraries
library(fda)
library(refund)
library(tidyverse)

# gasoline load data set
data(gasoline)
NIR <- as.matrix(gasoline$NIR)

# set up parameters
n_basis_fun <- 20
n_order <- 4
n_harm <- 5
stepsize <- 0.1

# generate fitting functional basis (b-spline basis)
smallbasis <- create.bspline.basis(
  rangeval = c(0, dim(NIR)[2]), nbasis = n_basis_fun,
  norder = n_order
)

# evaluate basis at grid points
basis_values <- eval.basis(
  evalarg = seq(from = 0, to = dim(NIR)[2], by = stepsize),
  basisobj = smallbasis
)

# translate observations in matrix into basis representation
NIR_fd <- smooth.basis(y = t(NIR), fdParobj = smallbasis)$fd

# get principal components
NIR_pcaObj <- pca.fd(fdobj = NIR_fd, nharm = n_harm, centerfns = TRUE)

# get coefficients for principal components
NIR_pc_coefs <- NIR_pcaObj$harmonics$coefs

# calculate values of principal components
NIR_pc_values <- basis_values %*% NIR_pc_coefs

# extract real wavelengths
wavelengths <- unlist(
  map(
    .x = colnames(NIR),
    .f = function(str) strtoi(strsplit(x = str, split = " nm"))
  )
)

nm <- (seq(from = 0, to = dim(NIR)[2], by = stepsize) / dim(NIR)[2]) * (max(wavelengths) - min(wavelengths)) + min(wavelengths)

### this is for plotting the principal components themselves ###
# form into tibble
NIR_pc_tibble <- cbind(
  nm, as_tibble(NIR_pc_values, .name_repair = "unique")
)
names(NIR_pc_tibble) <- c("nm", paste0("PC_", 1:n_harm))

# bring into format for plotting with ggplot
NIR_plotting_tibble <- NIR_pc_tibble %>%
  pivot_longer(cols = !nm, names_to = "PC", values_to = "Value")

# generate plot
principal_component_plot <- ggplot(data = NIR_plotting_tibble) +
  geom_line(aes(x = nm, y = Value, col = PC)) +
  # ggtitle("B-spline Basis") +
  theme_light() +
  # scale_colour_brewer(palette="Set1") +
  theme(
    legend.position = "bottom",
    # plot.title = element_text(size = 30),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20)
  )

# save plot in appropriate folder
ggsave(
  filename = "../Graphics/principal_components.pdf", plot = principal_component_plot,
  width = 20, height = 8, units = "in", dpi = 600
)

### from here on out it's about the approximated curves ###
# extract scores - weights for the approximation
NIR_pc_scores <- NIR_pcaObj$scores

# define helper function to get matrix format as wanted for further processing
approx_helper <- function(obs_index, score_mat, pc_value_mat) {
  # get scores for current observation
  obs_scores <- score_mat[obs_index, ]
  # create matrix container for approximations
  approx_mat <- matrix(data = NA, nrow = dim(pc_value_mat)[1], ncol = dim(pc_value_mat)[2])
  # fill in approximations
  for (i in 1:dim(pc_value_mat)[2]) {
    approx_mat[, i] <- colSums(obs_scores[1:i] %*% t(pc_value_mat[, 1:i]))
  }
  # return matrix
  return(approx_mat)
}

# get different approximations with increasing number of pc's included
NIR_approx <- map(
  .x = 1:dim(NIR)[1],
  .f = function(i) {
    approx_helper(
      obs_index = i, score_mat = NIR_pc_scores, pc_value_mat = NIR_pc_values
    )
  }
)

# extract single approximation matrix
approx_obs_tibble <- cbind(
  nm, as_tibble(NIR_approx[[1]], .name_repair = "unique")
)

names(approx_obs_tibble) <- c("nm", 1:n_harm)

# bring into format for plotting with ggplot
plotting_approx_obs_tibble <- approx_obs_tibble %>%
  pivot_longer(cols = !nm, names_to = "Num_PCs", values_to = "Value")

# generate plot
pc_approx_plot <- ggplot(data = plotting_approx_obs_tibble) +
  geom_line(aes(x = nm, y = Value, col = Num_PCs)) +
  # ggtitle("B-spline Basis") +
  theme_light() +
  # scale_colour_brewer(palette="Set1") +
  theme(
    legend.position = "bottom",
    # plot.title = element_text(size = 30),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20)
  )

# save plot in appropriate folder
ggsave(
  filename = "../Graphics/pc_approx.pdf", plot = pc_approx_plot,
  width = 20, height = 8, units = "in", dpi = 600
)