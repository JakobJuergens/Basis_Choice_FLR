##### This script generates an iid data set from wiener processes and plots it

# Load packages
library(tidyverse)

# setup for plot generation
n_points <- 10000
grid <- seq(from = 0, to = 1, length.out = n_points)
realizations <- 25
#smooth
f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
#bumpy
f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -  4*exp(-0.5*(grid-0.45)^2/0.015^2) +  8*exp(-0.5*(grid-0.6)^2/0.02^2) -  exp(-0.5*(grid-0.8)^2/0.03^2)
#two different variances of error

data_mat <- matrix(data = f2,
                   nrow = n_points, ncol = realizations, byrow = FALSE)

data_tibble <- cbind(tibble(x = points), as_tibble(data_mat, .name_repair = 'unique'))
names(data_tibble) <- c('t', paste0('Realization_', 1:realizations))

plot_tibble <- data_tibble %>% 
  pivot_longer(cols = !t, names_to = 'Realization', values_to = 'Value')

# generate plot
Wiener_plot <- ggplot(data = plot_tibble) +
  geom_line(aes(x = t, y = Value, col = Realization)) +
  #ggtitle("B-spline Basis") +
  theme_light() +
  #scale_colour_brewer(palette="Set1") + 
  theme(legend.position = "none",
        #plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

Wiener_plot

# save plot in appropriate folder
ggsave(filename = '../Graphics/Wiener_plot.pdf', plot = Wiener_plot,
       width = 20, height = 8, units = "in", dpi = 600)
