##### Load packages #####
library(tidyverse)
library(fda)
library(refund)
library(rlist)

###### load auxiliary files #####
source("data_generator.R")

###### original functions #####
f_1 <- function(t) {
  return(401*(2 * sin(0.5 * pi * t) + 4 * sin(1.5 * pi * t) + 5 * sin(2.5 * pi * t)))
}

f_2 <- function(t) {
  return(401*(1.5 * exp(-0.5 * (t - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (t - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (t - 0.6)^2 / 0.02^2) - exp(-0.5 * (t - 0.8)^2 / 0.03^2)))
}

# grid for plotting
plot_grid <- seq(0, 1, length.out = 1000)

###### generate response variables #####
grid <- seq(0, 1, length.out = 401)
f1 <- 2 * sin(0.5 * pi * grid) + 4 * sin(1.5 * pi * grid) + 5 * sin(2.5 * pi * grid)
f2 <- 1.5 * exp(-0.5 * (grid - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (grid - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (grid - 0.6)^2 / 0.02^2) - exp(-0.5 * (grid - 0.8)^2 / 0.03^2)

data(gasoline)
NIR <- gasoline$NIR

sigma_eps_squared1_1 <- as.numeric((var(NIR %*% f1) / 0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 <- as.numeric((var(NIR %*% f1) / 0.6) - var(NIR %*% f1))
sigma_eps_squared2_1 <- as.numeric((var(NIR %*% f2) / 0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 <- as.numeric((var(NIR %*% f2) / 0.6) - var(NIR %*% f2))

generated_curves <- NIR_curve_generator(n = 200, n_harmonics = 30)
my_data <- as.matrix(generated_curves[, -1])

Y1_1 <- as.numeric(my_data %*% f1 + rnorm(200, 0, 1) * sigma_eps_squared1_1)
Y1_2 <- as.numeric(my_data %*% f1 + rnorm(200, 0, 1) * sigma_eps_squared1_2)
Y2_1 <- as.numeric(my_data %*% f2 + rnorm(200, 0, 1) * sigma_eps_squared2_1)
Y2_2 <- as.numeric(my_data %*% f2 + rnorm(200, 0, 1) * sigma_eps_squared2_2)

##### create basis objects from optimal specifications #####

# constant for functional linear regression
betafdPar1 <- fdPar(fd(0, create.constant.basis(c(0, 1))))

# basis functions
basis_functions <- list(
  bspline_nharm2_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm3_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm4_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm5_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm6_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm7_50 = create.bspline.basis(rangeval = c(0, 1), nbasis = 50, norder = 4),
  bspline_nharm2_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4),
  bspline_nharm3_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4),
  bspline_nharm4_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4),
  bspline_nharm5_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4),
  bspline_nharm6_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4),
  bspline_nharm7_70 = create.bspline.basis(rangeval = c(0, 1), nbasis = 70, norder = 4)
)

##### evaluate bases at plot_grid #####

eval_bases <- map(
  .x = basis_functions,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

##### Construct fpcs with chosen specifications #####
pca <- map(
  .x = 1:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(my_data), fdParobj = basis_functions[[i]])$fd,
      nharm = (i - 1) %% 6 + 2, centerfns = TRUE
    )
  }
)

##### Extract Coefficients for fpcs #####
pca_basis_coef <- map(
  .x = pca,
  .f = function(fpca_obj) fpca_obj$harmonics$coefs
)

#### Multiply coefficients with basis_evaluations #####
pc_values <- map(
  .x = 1:12,
  .f = function(i) eval_bases[[i]] %*% pca_basis_coef[[i]]
)

##### Extract Scores from fpcs #####
pca_scores_coef <- map(
  .x = pca,
  .f = function(fpca_obj) fpca_obj$scores
)

##### Bind scores and responses to data frame for lm #####
pca_df_1_1 <- map(
  .x = pca_scores_coef,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y1_1, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_1_2 <- map(
  .x = pca_scores_coef,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y1_2, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_2_1 <- map(
  .x = pca_scores_coef,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y2_1, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_2_2 <- map(
  .x = pca_scores_coef,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y2_2, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

##### Perform lm with scores #####
pca_reg_1_1 <- map(
  .x = pca_df_1_1,
  .f = function(df) lm(formula = Response ~ ., data = df)
)

pca_reg_1_2 <- map(
  .x = pca_df_1_2,
  .f = function(df) lm(formula = Response ~ ., data = df)
)

pca_reg_2_1 <- map(
  .x = pca_df_2_1,
  .f = function(df) lm(formula = Response ~ ., data = df)
)

pca_reg_2_2 <- map(
  .x = pca_df_2_2,
  .f = function(df) lm(formula = Response ~ ., data = df)
)

##### Extract Estimates without intercept from lm objects ####
pca_reg_coef_1_1 <- map(
  .x = pca_reg_1_1,
  .f = function(lm_obj) unname(lm_obj$coefficients)[-1]
)

pca_reg_coef_1_2 <- map(
  .x = pca_reg_1_2,
  .f = function(lm_obj) unname(lm_obj$coefficients)[-1]
)

pca_reg_coef_2_1 <- map(
  .x = pca_reg_2_1,
  .f = function(lm_obj) unname(lm_obj$coefficients)[-1]
)

pca_reg_coef_2_2 <- map(
  .x = pca_reg_2_2,
  .f = function(lm_obj) unname(lm_obj$coefficients)[-1]
)

##### Multiply Estimates with principal components, scale
# and put into tibbles for plotting #####
pca_beta_est_1_1 <- as_tibble(list.cbind(map(
  .x = 1:12,
  .f = function(i) t(pca_reg_coef_1_1[[i]] %*% t(pc_values[[i]])) # / 401
)))

colnames(pca_beta_est_1_1) <- c(
  paste0("bspline_nharm", 2:7, "_50"),
  paste0("bspline_nharm", 2:7, "_70")
)

pca_beta_est_1_2 <- as_tibble(list.cbind(map(
  .x = 1:12,
  .f = function(i) t(pca_reg_coef_1_2[[i]] %*% t(pc_values[[i]])) # / 401
)))

colnames(pca_beta_est_1_2) <- c(
  paste0("bspline_nharm", 2:7, "_50"),
  paste0("bspline_nharm", 2:7, "_70")
)

pca_beta_est_2_1 <- as_tibble(list.cbind(map(
  .x = 1:12,
  .f = function(i) t(pca_reg_coef_2_1[[i]] %*% t(pc_values[[i]])) # / 401
)))

colnames(pca_beta_est_2_1) <- c(
  paste0("bspline_nharm", 2:7, "_50"),
  paste0("bspline_nharm", 2:7, "_70")
)

pca_beta_est_2_2 <- as_tibble(list.cbind(map(
  .x = 1:12,
  .f = function(i) t(pca_reg_coef_2_2[[i]] %*% t(pc_values[[i]])) # / 401
)))

colnames(pca_beta_est_2_2) <- c(
  paste0("bspline_nharm", 2:7, "_50"),
  paste0("bspline_nharm", 2:7, "_70")
)

##### Extract Intercepts from lm models #####
pca_reg_interc_1_1 <- map(
  .x = pca_reg_1_1,
  .f = function(lm_obj) unname(lm_obj$coefficients)[1]
)

pca_reg_interc_1_2 <- map(
  .x = pca_reg_1_2,
  .f = function(lm_obj) unname(lm_obj$coefficients)[1]
)

pca_reg_interc_2_1 <- map(
  .x = pca_reg_2_1,
  .f = function(lm_obj) unname(lm_obj$coefficients)[1]
)

pca_reg_interc_2_2 <- map(
  .x = pca_reg_2_2,
  .f = function(lm_obj) unname(lm_obj$coefficients)[1]
)

##### bind tibbles for plotting #####
plot_tibble_1_1 <- cbind(
  tibble(
    x = plot_grid,
    f1 = f_1(plot_grid)
  ),
  pca_beta_est_1_1
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f1", paste0("bspline_nharm", 2:7, "_50"), paste0("bspline_nharm", 2:7, "_70")
  )))

plot_tibble_1_2 <- cbind(
  tibble(
    x = plot_grid,
    f1 = f_1(plot_grid)
  ),
  pca_beta_est_1_2
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f1", paste0("bspline_nharm", 2:7, "_50"), paste0("bspline_nharm", 2:7, "_70")
  )))

plot_tibble_2_1 <- cbind(
  tibble(
    x = plot_grid,
    f2 = f_2(plot_grid)
  ),
  pca_beta_est_2_1
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f2", paste0("bspline_nharm", 2:7, "_50"), paste0("bspline_nharm", 2:7, "_70")
  )))

plot_tibble_2_2 <- cbind(
  tibble(
    x = plot_grid,
    f2 = f_2(plot_grid)
  ),
  pca_beta_est_2_2
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f2", paste0("bspline_nharm", 2:7, "_50"), paste0("bspline_nharm", 2:7, "_70")
  )))

##### generate plots for fpcr nharm = 2 #####
fpcr_nharm2_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm2_50", "bspline_nharm2_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm2_1_1.pdf", plot = fpcr_nharm2_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm2_50", "bspline_nharm2_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm2_1_2.pdf", plot = fpcr_nharm2_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm2_50", "bspline_nharm2_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm2_2_1.pdf", plot = fpcr_nharm2_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm2_50", "bspline_nharm2_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm2_2_2.pdf", plot = fpcr_nharm2_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 3 #####
fpcr_nharm3_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm3_50", "bspline_nharm3_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm3_1_1.pdf", plot = fpcr_nharm3_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm3_50", "bspline_nharm3_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm3_1_2.pdf", plot = fpcr_nharm3_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm3_50", "bspline_nharm3_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm3_2_1.pdf", plot = fpcr_nharm3_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm3_50", "bspline_nharm3_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm3_2_2.pdf", plot = fpcr_nharm3_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 4 #####
fpcr_nharm4_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm4_50", "bspline_nharm4_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm4_1_1.pdf", plot = fpcr_nharm4_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm4_50", "bspline_nharm4_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm4_1_2.pdf", plot = fpcr_nharm4_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm4_50", "bspline_nharm4_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm4_2_1.pdf", plot = fpcr_nharm4_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm4_50", "bspline_nharm4_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm4_2_2.pdf", plot = fpcr_nharm4_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 5 #####
fpcr_nharm5_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm5_50", "bspline_nharm5_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm5_1_1.pdf", plot = fpcr_nharm5_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm5_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm5_50", "bspline_nharm5_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm5_1_2.pdf", plot = fpcr_nharm5_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm5_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm5_50", "bspline_nharm5_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm5_2_1.pdf", plot = fpcr_nharm5_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm5_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm5_50", "bspline_nharm5_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm5_2_2.pdf", plot = fpcr_nharm5_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 6 #####
fpcr_nharm6_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm6_50", "bspline_nharm6_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm6_1_1.pdf", plot = fpcr_nharm6_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm6_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm6_50", "bspline_nharm6_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm6_1_2.pdf", plot = fpcr_nharm6_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm6_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm6_50", "bspline_nharm6_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm6_2_1.pdf", plot = fpcr_nharm6_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm6_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm6_50", "bspline_nharm6_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm6_2_2.pdf", plot = fpcr_nharm6_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 7 #####
fpcr_nharm7_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm7_50", "bspline_nharm7_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm7_1_1.pdf", plot = fpcr_nharm7_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm7_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm7_50", "bspline_nharm7_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm7_1_2.pdf", plot = fpcr_nharm7_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm7_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm7_50", "bspline_nharm7_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm7_2_1.pdf", plot = fpcr_nharm7_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm7_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm7_50", "bspline_nharm7_70"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 10)))

ggsave(
  filename = "../Graphics/Curve_Estimates/Large/fpcr_nharm7_2_2.pdf", plot = fpcr_nharm7_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)
