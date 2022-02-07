# Load packages
library(tidyverse)

# define polynomial functions
Legendre <- function(degree, x){
  if(degree == 0){return(1)}
  if(degree == 1){return(x)}
  if(degree == 2){return(1/2 * (3*x^2 - 1))}
  if(degree == 3){return(1/2 * (5*x^3 - 3*x))}
  if(degree == 4){return(1/8 * (35*x^4 - 30*x^2 + 3))}
  if(degree == 5){return(1/8 * (63*x^5 - 70*x^3 + 15*x))}
  if(degree == 6){return(1/16 * (231*x^6 - 315*x^4 + 105*x^2 - 5))}
  if(degree == 7){return(1/16 * (429*x^7 - 693*x^5 + 315*x^3 - 35*x))}
  else{stop('Higher degrees than 7 not implemented.')}
}

# create plotting tibble
plotting <- tibble(x = seq(0,1, length.out = 1000),
                   P_0 = Legendre(degree = 0, x = x),
                   P_1 = Legendre(degree = 1, x = x),
                   P_2 = Legendre(degree = 2, x = x),
                   P_3 = Legendre(degree = 3, x = x),
                   P_4 = Legendre(degree = 4, x = x),
                   P_5 = Legendre(degree = 5, x = x),
                   P_6 = Legendre(degree = 6, x = x),
                   P_7 = Legendre(degree = 7, x = x)) %>% 
  pivot_longer(cols = !x, names_to = 'Polynomial', values_to = 'y')


Legendre_plot <- ggplot(data = plotting) +
  geom_line(aes(x = x, y = y, col = Polynomial)) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20))

# save plot in appropriate folder
ggsave(filename = '../Graphics/Legendre_plot.pdf', plot = Legendre_plot,
       width = 20, height = 8, units = "in", dpi = 600)