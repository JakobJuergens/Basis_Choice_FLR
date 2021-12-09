##### This script generates a plot showing of a chosen number of fourier basis functions #####

# Load packages
library(tidyverse)

# setup for plot generation
n_funcs <- 6
n_points <- 1000

# function factory that generates the basis functions
basis_func <- function(i){
  if(i == 1){
    ret_func <- function(x){return(1)}
  } else
  if(i %% 2 == 0){
      ret_func <- function(x){return(sqrt(2)*cos(pi * i * x))}
  } else
  if(i > 1 && i %% 2 == 1){
    ret_func <- function(x){return(sqrt(2)*sin(pi * (i-1) * x))}
  } else{
    stop('Malformed Index chosen.')
  }
}

# Create points for calculation of basis functions
points <- seq(from = 0, to = 1, length.out = n_points)

# Calculate basis functions at points
res_list <- list()
for(i in 1:n_funcs){
  tmp_func <- basis_func(i)
  res_list[[i]] <- unlist(
    map(.x = points, .f = tmp_func)
    )
}

# convert to tibble
res_tibble <- cbind(tibble(x = points),
                    as_tibble(res_list, .name_repair = 'unique'))
names(res_tibble) <- c('x', paste0('S_', 1:n_funcs))

# pivot_longer to allow for easy plotting
plot_tibble <- res_tibble %>% 
  pivot_longer(cols = !x, names_to = 'Function', values_to = 'Value')

# generate plot
fourier_plot <- ggplot(data = plot_tibble) +
  geom_line(aes(x = x, y = Value, col = Function)) +
  #ggtitle("Fourier Basis") +
  theme_light() +
  #scale_colour_brewer(palette="Set1") + 
  theme(legend.position = "none",
        #plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

# save plot in appropriate folder
ggsave(filename = '../Graphics/Fourier_Basis.pdf', plot = fourier_plot,
       width = 20, height = 8, units = "in", dpi = 600)
