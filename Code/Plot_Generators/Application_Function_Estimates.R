##### Load packages #####
library(tidyverse)
library(fda)
library(refund)
library(rlist)

# load data
data(gasoline)
NIR <- gasoline$NIR
octane <- gasoline$octane
grid <- seq(0, 1, length.out = 401)

# grid for plotting
plot_grid <- seq(0, 1, length.out = 1000)

# basis functions
basis_functions <- list(
  bspline_exp = create.bspline.basis(rangeval = c(0, 1), nbasis = 10, norder = 4),
  monomial_exp = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_exp = create.fourier.basis(rangeval = c(0, 1), nbasis = 9),
  bspline_nharm2 = create.bspline.basis(rangeval = c(0, 1), nbasis = 19, norder = 4),
  bspline_nharm3 = create.bspline.basis(rangeval = c(0, 1), nbasis = 23, norder = 4),
  bspline_nharm4 = create.bspline.basis(rangeval = c(0, 1), nbasis = 25, norder = 4),
  monomial_nharm2 = create.monomial.basis(rangeval = c(0, 1), nbasis = 3),
  monomial_nharm3 = create.monomial.basis(rangeval = c(0, 1), nbasis = 6),
  monomial_nharm4 = create.monomial.basis(rangeval = c(0, 1), nbasis = 11),
  fourier_nharm2 = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_nharm3 = create.fourier.basis(rangeval = c(0, 1), nbasis = 9),
  fourier_nharm4 = create.fourier.basis(rangeval = c(0, 1), nbasis = 7)
)

##### evaluate bases at plot_grid #####
eval_bases <- map(
  .x = basis_functions,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

##### create basis objects from optimal specifications #####

# constant for functional linear regression
betafdPar1 <- fdPar(fd(0, create.constant.basis(c(0, 1))))

##### create betalist_list for basis expansion regression #####
betalist_list <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = betafdPar1,
      smooth_basis = fdPar(basis_functions[[i]])
    )
  }
)

##### create xfdlist_list for basis expansion regression #####
xfdlist_list <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = rep(x = 1, times = 60),
      smooth_basis = smooth.basis(
        argvals = grid, y = t(NIR), fdParobj = basis_functions[[i]]
      )$fd
    )
  }
)

##### estimate basis expansion models #####
fReg <- map(
  .x = 1:3,
  .f = function(i) {
    fRegress(
      y = octane, xfdlist = xfdlist_list[[i]], betalist = betalist_list[[i]]
    )
  }
)

##### extract intercepts parameters from basis expansion models #####
fReg_interc <- map(
  .x = fReg,
  .f = function(model) model$betaestlist$const$fd$coefs[1]
)

##### extract coefficients from basis expansion models #####
fReg_coef <- map(
  .x = fReg,
  .f = function(model) model$betaestlist$smooth_basis$fd$coefs
)

###### calculate estimated coefficient functions for basis expansion models
##### and bind to data.frame #####
beta_hat <- as_tibble(list.cbind(map(
  .x = 1:3,
  .f = function(i) (eval_bases[[i]] %*% fReg_coef[[i]]) # / 401
)))

names(beta_hat) <- c("B-Spline", "Monomial", "Fourier")

##### Construct fpcs with chosen specifications #####
pca <- map(
  .x = 4:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(NIR), fdParobj = basis_functions[[i]])$fd,
      nharm = (i - 4) %% 3 + 2, centerfns = TRUE
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
  .x = 1:9,
  .f = function(i) eval_bases[[i + 3]] %*% pca_basis_coef[[i]]
)

##### Extract Scores from fpcs #####
pca_scores_coef <- map(
  .x = pca,
  .f = function(fpca_obj) fpca_obj$scores
)

##### Bind scores and responses to data frame for lm #####
pca_df <- map(
  .x = pca_scores_coef,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(octane, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

##### Perform lm with scores #####
pca_reg <- map(
  .x = pca_df,
  .f = function(df) lm(formula = Response ~ ., data = df)
)

##### Extract Estimates without intercept from lm objects ####
pca_reg_coef <- map(
  .x = pca_reg,
  .f = function(lm_obj) unname(lm_obj$coefficients)[-1]
)

##### Multiply Estimates with principal components, scale
# and put into tibbles for plotting #####
pca_beta_est <- as_tibble(list.cbind(map(
  .x = 1:9,
  .f = function(i) t(pca_reg_coef[[i]] %*% t(pc_values[[i]])) # / 401
)))

colnames(pca_beta_est) <- c(
  paste0("bspline_nharm", 2:4),
  paste0("monomial_nharm", 2:4),
  paste0("fourier_nharm", 2:4)
)

##### Extract Intercepts from lm models #####
pca_reg_interc <- map(
  .x = pca_reg,
  .f = function(lm_obj) unname(lm_obj$coefficients)[1]
)

##### bind tibbles for plotting #####
plot_tibble <- cbind(
  tibble(
    x = plot_grid,
  ),
  beta_hat,
  pca_beta_est
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "B-Spline", "Monomial", "Fourier",
    "bspline_nharm2", "bspline_nharm3", "bspline_nharm4",
    "monomial_nharm2", "monomial_nharm3", "monomial_nharm4",
    "fourier_nharm2", "fourier_nharm3", "fourier_nharm4"
  )))

##### generate plots for basis expansion #####
basis_expansion_plot <- ggplot(data = plot_tibble %>%
                                     filter(Curve %in% c("B-Spline", "Monomial", "Fourier"))) +
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
  filename = "../Graphics/Appl_Curve_Estimates/appl_basis_expansion.pdf", plot = basis_expansion_plot,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 2 #####
fpcr_nharm2_plot <- ggplot(data = plot_tibble %>%
                                 filter(Curve %in% c("bspline_nharm2", "monomial_nharm2", "fourier_nharm2"))) +
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
  filename = "../Graphics/Appl_Curve_Estimates/fpcr_nharm2.pdf", plot = fpcr_nharm2_plot,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 3 #####
fpcr_nharm3_plot <- ggplot(data = plot_tibble %>%
                                 filter(Curve %in% c("bspline_nharm3", "monomial_nharm3", "fourier_nharm3"))) +
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
  filename = "../Graphics/Appl_Curve_Estimates/fpcr_nharm3.pdf", plot = fpcr_nharm3_plot,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 4 #####
fpcr_nharm4_plot <- ggplot(data = plot_tibble %>%
                                 filter(Curve %in% c("bspline_nharm4", "monomial_nharm4", "fourier_nharm4"))) +
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
  filename = "../Graphics/Appl_Curve_Estimates/fpcr_nharm4.pdf", plot = fpcr_nharm4_plot,
  width = 20, height = 12, units = "in", dpi = 600
)