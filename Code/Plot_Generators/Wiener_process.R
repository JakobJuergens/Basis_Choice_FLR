##### This script generates an iid data set from wiener processes and plots it

# Load packages
library(tidyverse)

# setup for plot generation
n_points <- 10000
points <- seq(from = 0, to = 1, length.out = n_points)
realizations <- 25

Wiener_process <- function(grid){
  values <- rep(0, times = length(grid))
  for(i in 2:length(grid)){
    values[i] <- values[i-1] + rnorm(n = 1, mean = 0, sd = sqrt(grid[i] - grid[i-1]))
  }
  return(values)
}

data_mat <- matrix(data = unlist(map(.x = 1:realizations,
                                     .f = function(i) Wiener_process(points))),
                   nrow = n_points, ncol = realizations, byrow = FALSE)

data_tibble <- cbind(tibble(x = points), as_tibble(data_mat, .name_repair = 'unique'))
names(data_tibble) <- c('t', paste0('Realization_', 1:realizations))

plot_tibble <- data_tibble %>% 
  pivot_longer(cols = !x, names_to = 'Realization', values_to = 'Value')

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

# save plot in appropriate folder
ggsave(filename = '../Graphics/Wiener_plot.pdf', plot = Wiener_plot,
       width = 20, height = 8, units = "in", dpi = 600)
