library(tidyverse)

f1 <- function(t){
  return(401*(2*sin(0.5*pi*t) + 4*sin(1.5*pi*t) + 5*sin(2.5*pi*t)))
}

f2 <- function(t){
  return(401*(1.5 * exp(-0.5 * (t - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (t - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (t - 0.6)^2 / 0.02^2) - exp(-0.5 * (t - 0.8)^2 / 0.03^2)))
}

grid <- seq(0, 1, length.out = 1000)

plot_tibble_1 <- tibble(x = grid,
                        y = f1(grid))

plot_tibble_2 <- tibble(x = grid,
                        y = f2(grid))

f1_plot <- ggplot(data = plot_tibble_1) +
  geom_line(aes(x = x, y = y), col = 'red') +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40),
        legend.title = element_text(size = 48),
        legend.text = element_text(size = 40))

# save plot in appropriate folder
ggsave(filename = '../Graphics/f1_plot.pdf', plot = f1_plot,
       width = 20, height = 12, units = "in", dpi = 600)

f2_plot <- ggplot(data = plot_tibble_2) +
  geom_line(aes(x = x, y = y), col = 'red') +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40),
        legend.title = element_text(size = 48),
        legend.text = element_text(size = 40))

# save plot in appropriate folder
ggsave(filename = '../Graphics/f2_plot.pdf', plot = f2_plot,
       width = 20, height = 12, units = "in", dpi = 600)
