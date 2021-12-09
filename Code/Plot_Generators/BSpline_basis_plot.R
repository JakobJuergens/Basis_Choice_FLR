##### This script generates a plot showing of a chosen number of fourier basis functions #####

# Load packages
library(tidyverse)
library(fda)

# setup for plot generation
n_points <- 1000
knots <- seq(from = 0, to = 1, length.out = 8)
order <- 3

# Create points for calculation of basis functions
points <- seq(from = 0, to = 1, length.out = n_points)

# function factory that generates the basis functions
bspline_basis <- create.bspline.basis(rangeval = c(0,1), norder = order,
                                        breaks = knots)

# Calculate basis functions at points
res_matrix <- eval.basis(evalarg = points, basisobj = bspline_basis)

# convert to tibble
res_tibble <- cbind(x = points, as_tibble(res_matrix))
names(res_tibble) <- c('x', paste0('S_', 1:dim(res_matrix)[2]))

# pivot_longer to allow for easy plotting
plot_tibble <- res_tibble %>% 
  pivot_longer(cols = !x, names_to = 'Function', values_to = 'Value')

# generate plot
bspline_plot <- ggplot(data = plot_tibble) +
  geom_line(aes(x = x, y = Value, col = Function)) +
  ggtitle("B-spline Basis") +
  theme_light() +
  #scale_colour_brewer(palette="Set1") + 
  theme(legend.position="bottom",
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

# save plot in appropriate folder
ggsave(filename = '../Graphics/Bspline_Basis.pdf', plot = bspline_plot,
       width = 20, height = 8, units = "in", dpi = 600)
