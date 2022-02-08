# Clear Workspace
rm(list = ls())

# load libraries
library(fda)
library(refund)
library(tidyverse)

data(gasoline)
NIR <- gasoline$NIR

my_basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 100, norder = 4)

NIR_func <- smooth.basis(argvals = seq(0, 1, length.out = 401), 
                         y = t(NIR), fdParobj = my_basis)

plot_grid <- seq(0, 1, length.out = 1000)
NIR_plot_tibble <- cbind(
  tibble(x = plot_grid),
  as_tibble(eval.fd(evalarg = plot_grid, fdobj = NIR_func$fd)))

colnames(NIR_plot_tibble) <- c('x', paste0('Observation_', 1:60))

NIR_plot_tibble <- NIR_plot_tibble %>% 
  pivot_longer(cols = !x, names_to = 'Observation', values_to = 'y')

NIR_plot <- ggplot(data = NIR_plot_tibble) +
  geom_line(aes(x = x, y = y, col = Observation), alpha = 0.5) +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
  )

# save plot in appropriate folder
ggsave(
  filename = "../Graphics/NIR_data.pdf", plot = NIR_plot,
  width = 20, height = 12, units = "in", dpi = 600
)