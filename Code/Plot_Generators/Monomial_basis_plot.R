##### This script generates a plot showing of a chosen number of fourier basis functions #####

# Load packages
library(tidyverse)
library(fda)

# setup for plot generation
n_points <- 1000

# Create points for calculation of basis functions
points <- seq(from = 0, to = 1, length.out = n_points)

# function factory that generates the basis functions
monomial_basis <- create.monomial.basis(rangeval = c(0,1), nbasis = 8)

# Calculate basis functions at points
res_matrix <- eval.basis(evalarg = points, basisobj = monomial_basis)

# convert to tibble
res_tibble <- cbind(x = points, as_tibble(res_matrix))
names(res_tibble) <- c('x', paste0('S_', 1:dim(res_matrix)[2]))

# pivot_longer to allow for easy plotting
plot_tibble <- res_tibble %>% 
  pivot_longer(cols = !x, names_to = 'Function', values_to = 'Value')

# generate plot
monomial_plot <- ggplot(data = plot_tibble) +
  geom_line(aes(x = x, y = Value, col = Function)) +
  #ggtitle("B-spline Basis") +
  theme_light() +
  #scale_colour_brewer(palette="Set1") + 
  theme(legend.position = "none",
        #plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

# save plot in appropriate folder
ggsave(filename = '../Graphics/Monomial_Basis.pdf', plot = monomial_plot,
       width = 20, height = 8, units = "in", dpi = 600)
