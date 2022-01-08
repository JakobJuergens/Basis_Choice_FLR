##### This script plots the true coefficient functions from Reiss and Ogden (2007)
# Clear Workspace
rm(list = ls())

# Load packages
library(tidyverse)

# setup for plot generation
n_points <- 10000
grid <- seq(from = 0, to = 1, length.out = n_points)
realizations <- 25
# smooth
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
# bumpy
f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) - exp(-0.5 * (grid - 0.8)^2 / 0.03^2)

# Create tibble containing grid points and function values
coef_func_tibble <- tibble(t = grid, `f1(t)` = f1, `f2(t)` = f2)

# Pivot tibble for plotting
plot_tibble <- coef_func_tibble %>%
  pivot_longer(cols = !t, names_to = "Coefficient-Function", values_to = "Value")

# generate plot
Coef_Functions_plot <- ggplot(data = plot_tibble) +
  geom_line(aes(x = t, y = Value, col = `Coefficient-Function`)) +
  # ggtitle("B-spline Basis") +
  theme_light() +
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
  filename = "../Graphics/Coef_Functions_plot.pdf", plot = Coef_Functions_plot,
  width = 20, height = 8, units = "in", dpi = 600
)
