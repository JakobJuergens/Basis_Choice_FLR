##### Load packages #####
library(tidyverse)
library(fda)
library(refund)
library(rlist)

###### load auxiliary files #####
source("data_generator.R")

###### original functions #####
f_1 <- function(t) {
  return(2 * sin(0.5 * pi * t) + 4 * sin(1.5 * pi * t) + 5 * sin(2.5 * pi * t))
}

f_2 <- function(t) {
  return(1.5 * exp(-0.5 * (t - 0.3)^2 / 0.02^2) - 4 * exp(-0.5 * (t - 0.45)^2 / 0.015^2) + 8 * exp(-0.5 * (t - 0.6)^2 / 0.02^2) - exp(-0.5 * (t - 0.8)^2 / 0.03^2))
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
basis_functions_1_1 <- list(
  bspline_exp = create.bspline.basis(rangeval = c(0, 1), nbasis = 5, norder = 4),
  monomial_exp = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_exp = create.fourier.basis(rangeval = c(0, 1), nbasis = 5),
  bspline_nharm2 = create.bspline.basis(rangeval = c(0, 1), nbasis = 5, norder = 4),
  bspline_nharm3 = create.bspline.basis(rangeval = c(0, 1), nbasis = 5, norder = 4),
  bspline_nharm4 = create.bspline.basis(rangeval = c(0, 1), nbasis = 6, norder = 4),
  monomial_nharm2 = create.monomial.basis(rangeval = c(0, 1), nbasis = 4),
  monomial_nharm3 = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  monomial_nharm4 = create.monomial.basis(rangeval = c(0, 1), nbasis = 6),
  fourier_nharm2 = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_nharm3 = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_nharm4 = create.fourier.basis(rangeval = c(0, 1), nbasis = 5)
)

basis_functions_1_2 <- list(
  bspline_exp = create.bspline.basis(rangeval = c(0, 1), nbasis = 4, norder = 4),
  monomial_exp = create.monomial.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_exp = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  bspline_nharm2 = create.bspline.basis(rangeval = c(0, 1), nbasis = 5, norder = 4),
  bspline_nharm3 = create.bspline.basis(rangeval = c(0, 1), nbasis = 5, norder = 4),
  bspline_nharm4 = create.bspline.basis(rangeval = c(0, 1), nbasis = 6, norder = 4),
  monomial_nharm2 = create.monomial.basis(rangeval = c(0, 1), nbasis = 4),
  monomial_nharm3 = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  monomial_nharm4 = create.monomial.basis(rangeval = c(0, 1), nbasis = 6),
  fourier_nharm2 = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_nharm3 = create.fourier.basis(rangeval = c(0, 1), nbasis = 3),
  fourier_nharm4 = create.fourier.basis(rangeval = c(0, 1), nbasis = 5)
)

basis_functions_2_1 <- list(
  bspline_exp = create.bspline.basis(rangeval = c(0, 1), nbasis = 11, norder = 4),
  monomial_exp = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_exp = create.fourier.basis(rangeval = c(0, 1), nbasis = 9),
  bspline_nharm2 = create.bspline.basis(rangeval = c(0, 1), nbasis = 4, norder = 4),
  bspline_nharm3 = create.bspline.basis(rangeval = c(0, 1), nbasis = 6, norder = 4),
  bspline_nharm4 = create.bspline.basis(rangeval = c(0, 1), nbasis = 23, norder = 4),
  monomial_nharm2 = create.monomial.basis(rangeval = c(0, 1), nbasis = 10),
  monomial_nharm3 = create.monomial.basis(rangeval = c(0, 1), nbasis = 6),
  monomial_nharm4 = create.monomial.basis(rangeval = c(0, 1), nbasis = 8),
  fourier_nharm2 = create.fourier.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_nharm3 = create.fourier.basis(rangeval = c(0, 1), nbasis = 15),
  fourier_nharm4 = create.fourier.basis(rangeval = c(0, 1), nbasis = 7)
)

basis_functions_2_2 <- list(
  bspline_exp = create.bspline.basis(rangeval = c(0, 1), nbasis = 6, norder = 4),
  monomial_exp = create.monomial.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_exp = create.fourier.basis(rangeval = c(0, 1), nbasis = 7),
  bspline_nharm2 = create.bspline.basis(rangeval = c(0, 1), nbasis = 4, norder = 4),
  bspline_nharm3 = create.bspline.basis(rangeval = c(0, 1), nbasis = 6, norder = 4),
  bspline_nharm4 = create.bspline.basis(rangeval = c(0, 1), nbasis = 23, norder = 4),
  monomial_nharm2 = create.monomial.basis(rangeval = c(0, 1), nbasis = 10),
  monomial_nharm3 = create.monomial.basis(rangeval = c(0, 1), nbasis = 6),
  monomial_nharm4 = create.monomial.basis(rangeval = c(0, 1), nbasis = 8),
  fourier_nharm2 = create.fourier.basis(rangeval = c(0, 1), nbasis = 5),
  fourier_nharm3 = create.fourier.basis(rangeval = c(0, 1), nbasis = 15),
  fourier_nharm4 = create.fourier.basis(rangeval = c(0, 1), nbasis = 7)
)

##### evaluate bases at plot_grid #####

eval_bases_1_1 <- map(
  .x = basis_functions_1_1,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

eval_bases_1_2 <- map(
  .x = basis_functions_1_2,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

eval_bases_2_1 <- map(
  .x = basis_functions_2_1,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

eval_bases_2_2 <- map(
  .x = basis_functions_2_2,
  .f = function(bas) eval.basis(evalarg = plot_grid, basisobj = bas)
)

##### create betalist_list for basis expansion regression #####
betalist_list_1_1 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = betafdPar1,
      smooth_basis = fdPar(basis_functions_1_1[[i]])
    )
  }
)

betalist_list_1_2 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = betafdPar1,
      smooth_basis = fdPar(basis_functions_1_2[[i]])
    )
  }
)

betalist_list_2_1 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = betafdPar1,
      smooth_basis = fdPar(basis_functions_2_1[[i]])
    )
  }
)

betalist_list_2_2 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = betafdPar1,
      smooth_basis = fdPar(basis_functions_2_2[[i]])
    )
  }
)

##### create xfdlist_list for basis expansion regression #####
xfdlist_list_1_1 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = rep(x = 1, times = 200),
      smooth_basis = smooth.basis(
        argvals = grid, y = t(my_data), fdParobj = basis_functions_1_1[[i]]
      )$fd
    )
  }
)

xfdlist_list_1_2 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = rep(x = 1, times = 200),
      smooth_basis = smooth.basis(
        argvals = grid, y = t(my_data), fdParobj = basis_functions_1_2[[i]]
      )$fd
    )
  }
)

xfdlist_list_2_1 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = rep(x = 1, times = 200),
      smooth_basis = smooth.basis(
        argvals = grid, y = t(my_data), fdParobj = basis_functions_2_1[[i]]
      )$fd
    )
  }
)

xfdlist_list_2_2 <- map(
  .x = 1:3,
  .f = function(i) {
    list(
      const = rep(x = 1, times = 200),
      smooth_basis = smooth.basis(
        argvals = grid, y = t(my_data), fdParobj = basis_functions_2_2[[i]]
      )$fd
    )
  }
)

##### estimate basis expansion models #####
fReg_1_1 <- map(
  .x = 1:3,
  .f = function(i) {
    fRegress(
      y = Y1_1, xfdlist = xfdlist_list_1_1[[i]], betalist = betalist_list_1_1[[i]]
    )
  }
)

fReg_1_2 <- map(
  .x = 1:3,
  .f = function(i) {
    fRegress(
      y = Y1_2, xfdlist = xfdlist_list_1_2[[i]], betalist = betalist_list_1_2[[i]]
    )
  }
)

fReg_2_1 <- map(
  .x = 1:3,
  .f = function(i) {
    fRegress(
      y = Y2_1, xfdlist = xfdlist_list_2_1[[i]], betalist = betalist_list_2_1[[i]]
    )
  }
)

fReg_2_2 <- map(
  .x = 1:3,
  .f = function(i) {
    fRegress(
      y = Y2_2, xfdlist = xfdlist_list_2_2[[i]], betalist = betalist_list_2_2[[i]]
    )
  }
)

##### extract intercepts parameters from basis expansion models #####

fReg_interc_1_1 <- map(
  .x = fReg_1_1,
  .f = function(model) model$betaestlist$const$fd$coefs[1]
)

fReg_interc_1_2 <- map(
  .x = fReg_1_2,
  .f = function(model) model$betaestlist$const$fd$coefs[1]
)

fReg_interc_2_1 <- map(
  .x = fReg_2_1,
  .f = function(model) model$betaestlist$const$fd$coefs[1]
)

fReg_interc_2_2 <- map(
  .x = fReg_2_2,
  .f = function(model) model$betaestlist$const$fd$coefs[1]
)

##### extract coefficients from basis expansion models #####

fReg_coef_1_1 <- map(
  .x = fReg_1_1,
  .f = function(model) model$betaestlist$smooth_basis$fd$coefs
)

fReg_coef_1_2 <- map(
  .x = fReg_1_2,
  .f = function(model) model$betaestlist$smooth_basis$fd$coefs
)

fReg_coef_2_1 <- map(
  .x = fReg_2_1,
  .f = function(model) model$betaestlist$smooth_basis$fd$coefs
)

fReg_coef_2_2 <- map(
  .x = fReg_2_2,
  .f = function(model) model$betaestlist$smooth_basis$fd$coefs
)

###### calculate estimated coefficient functions for basis expansion models
##### and bind to data.frame #####
beta_hat_1_1 <- as_tibble(list.cbind(map(
  .x = 1:3,
  .f = function(i) (eval_bases_1_1[[i]] %*% fReg_coef_1_1[[i]]) / 401
)))

beta_hat_1_2 <- as_tibble(list.cbind(map(
  .x = 1:3,
  .f = function(i) (eval_bases_1_2[[i]] %*% fReg_coef_1_2[[i]]) / 401
)))

beta_hat_2_1 <- as_tibble(list.cbind(map(
  .x = 1:3,
  .f = function(i) (eval_bases_2_1[[i]] %*% fReg_coef_2_1[[i]]) / 401
)))

beta_hat_2_2 <- as_tibble(list.cbind(map(
  .x = 1:3,
  .f = function(i) (eval_bases_2_2[[i]] %*% fReg_coef_2_2[[i]]) / 401
)))

names(beta_hat_1_1) <- c("B-Spline", "Monomial", "Fourier")
names(beta_hat_1_2) <- c("B-Spline", "Monomial", "Fourier")
names(beta_hat_2_1) <- c("B-Spline", "Monomial", "Fourier")
names(beta_hat_2_2) <- c("B-Spline", "Monomial", "Fourier")

##### Construct fpcs with chosen specifications #####
pca_1_1 <- map(
  .x = 4:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(my_data), fdParobj = basis_functions_1_1[[i]])$fd,
      nharm = (i - 4) %% 3 + 2, centerfns = TRUE
    )
  }
)

pca_1_2 <- map(
  .x = 4:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(my_data), fdParobj = basis_functions_1_2[[i]])$fd,
      nharm = (i - 4) %% 3 + 2, centerfns = TRUE
    )
  }
)

pca_2_1 <- map(
  .x = 4:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(my_data), fdParobj = basis_functions_2_1[[i]])$fd,
      nharm = (i - 4) %% 3 + 2, centerfns = TRUE
    )
  }
)

pca_2_2 <- map(
  .x = 4:12,
  .f = function(i) {
    pca.fd(
      fdobj = smooth.basis(argvals = grid, y = t(my_data), fdParobj = basis_functions_2_2[[i]])$fd,
      nharm = (i - 4) %% 3 + 2, centerfns = TRUE
    )
  }
)

##### Extract Coefficients for fpcs #####
pca_basis_coef_1_1 <- map(
  .x = pca_1_1,
  .f = function(fpca_obj) fpca_obj$harmonics$coefs
)

pca_basis_coef_1_2 <- map(
  .x = pca_1_2,
  .f = function(fpca_obj) fpca_obj$harmonics$coefs
)

pca_basis_coef_2_1 <- map(
  .x = pca_2_1,
  .f = function(fpca_obj) fpca_obj$harmonics$coefs
)

pca_basis_coef_2_2 <- map(
  .x = pca_2_2,
  .f = function(fpca_obj) fpca_obj$harmonics$coefs
)

#### Multiply coefficients with basis_evaluations #####
pc_values_1_1 <- map(
  .x = 1:9,
  .f = function(i) eval_bases_1_1[[i + 3]] %*% pca_basis_coef_1_1[[i]]
)

pc_values_1_2 <- map(
  .x = 1:9,
  .f = function(i) eval_bases_1_2[[i + 3]] %*% pca_basis_coef_1_2[[i]]
)

pc_values_2_1 <- map(
  .x = 1:9,
  .f = function(i) eval_bases_2_1[[i + 3]] %*% pca_basis_coef_2_1[[i]]
)

pc_values_2_2 <- map(
  .x = 1:9,
  .f = function(i) eval_bases_2_2[[i + 3]] %*% pca_basis_coef_2_2[[i]]
)

##### Extract Scores from fpcs #####
pca_scores_coef_1_1 <- map(
  .x = pca_1_1,
  .f = function(fpca_obj) fpca_obj$scores
)

pca_scores_coef_1_2 <- map(
  .x = pca_1_2,
  .f = function(fpca_obj) fpca_obj$scores
)

pca_scores_coef_2_1 <- map(
  .x = pca_2_1,
  .f = function(fpca_obj) fpca_obj$scores
)

pca_scores_coef_2_2 <- map(
  .x = pca_2_2,
  .f = function(fpca_obj) fpca_obj$scores
)

##### Bind scores and responses to data frame for lm #####
pca_df_1_1 <- map(
  .x = pca_scores_coef_1_1,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y1_1, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_1_2 <- map(
  .x = pca_scores_coef_1_2,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y1_2, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_2_1 <- map(
  .x = pca_scores_coef_2_1,
  .f = function(scrs) {
    tmp <- as_tibble(cbind(Y2_1, scrs))
    colnames(tmp) <- c("Response", paste0("harm", 1:ncol(scrs)))
    return(tmp)
  }
)

pca_df_2_2 <- map(
  .x = pca_scores_coef_2_2,
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
  .x = 1:9,
  .f = function(i) t(pca_reg_coef_1_1[[i]] %*% t(pc_values_1_1[[i]])) / 401
)))

colnames(pca_beta_est_1_1) <- c(
  paste0("bspline_nharm", 2:4),
  paste0("monomial_nharm", 2:4),
  paste0("fourier_nharm", 2:4)
)

pca_beta_est_1_2 <- as_tibble(list.cbind(map(
  .x = 1:9,
  .f = function(i) t(pca_reg_coef_1_2[[i]] %*% t(pc_values_1_2[[i]])) / 401
)))

colnames(pca_beta_est_1_2) <- c(
  paste0("bspline_nharm", 2:4),
  paste0("monomial_nharm", 2:4),
  paste0("fourier_nharm", 2:4)
)

pca_beta_est_2_1 <- as_tibble(list.cbind(map(
  .x = 1:9,
  .f = function(i) t(pca_reg_coef_2_1[[i]] %*% t(pc_values_2_1[[i]])) / 401
)))

colnames(pca_beta_est_2_1) <- c(
  paste0("bspline_nharm", 2:4),
  paste0("monomial_nharm", 2:4),
  paste0("fourier_nharm", 2:4)
)

pca_beta_est_2_2 <- as_tibble(list.cbind(map(
  .x = 1:9,
  .f = function(i) t(pca_reg_coef_2_2[[i]] %*% t(pc_values_2_2[[i]])) / 401
)))

colnames(pca_beta_est_2_2) <- c(
  paste0("bspline_nharm", 2:4),
  paste0("monomial_nharm", 2:4),
  paste0("fourier_nharm", 2:4)
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
  beta_hat_1_1,
  pca_beta_est_1_1
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f1", "B-Spline", "Monomial", "Fourier",
    "bspline_nharm2", "bspline_nharm3", "bspline_nharm4",
    "monomial_nharm2", "monomial_nharm3", "monomial_nharm4",
    "fourier_nharm2", "fourier_nharm3", "fourier_nharm4"
  )))

plot_tibble_1_2 <- cbind(
  tibble(
    x = plot_grid,
    f1 = f_1(plot_grid)
  ),
  beta_hat_1_2,
  pca_beta_est_1_2
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f1", "B-Spline", "Monomial", "Fourier",
    "bspline_nharm2", "bspline_nharm3", "bspline_nharm4",
    "monomial_nharm2", "monomial_nharm3", "monomial_nharm4",
    "fourier_nharm2", "fourier_nharm3", "fourier_nharm4"
  )))

plot_tibble_2_1 <- cbind(
  tibble(
    x = plot_grid,
    f2 = f_2(plot_grid)
  ),
  beta_hat_2_1,
  pca_beta_est_2_1
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f2", "B-Spline", "Monomial", "Fourier",
    "bspline_nharm2", "bspline_nharm3", "bspline_nharm4",
    "monomial_nharm2", "monomial_nharm3", "monomial_nharm4",
    "fourier_nharm2", "fourier_nharm3", "fourier_nharm4"
  )))

plot_tibble_2_2 <- cbind(
  tibble(
    x = plot_grid,
    f2 = f_2(plot_grid)
  ),
  beta_hat_2_2,
  pca_beta_est_2_2
) %>%
  pivot_longer(cols = !x, names_to = "Curve", values_to = "y") %>%
  mutate(Curve = factor(Curve, levels = c(
    "f2", "B-Spline", "Monomial", "Fourier",
    "bspline_nharm2", "bspline_nharm3", "bspline_nharm4",
    "monomial_nharm2", "monomial_nharm3", "monomial_nharm4",
    "fourier_nharm2", "fourier_nharm3", "fourier_nharm4"
  )))

##### generate plots for basis expansion #####
basis_expansion_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
  filter(Curve %in% c("f1", "B-Spline", "Monomial", "Fourier"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/basis_expansion_1_1.pdf", plot = basis_expansion_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

basis_expansion_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
  filter(Curve %in% c("f1", "B-Spline", "Monomial", "Fourier"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/basis_expansion_1_2.pdf", plot = basis_expansion_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

basis_expansion_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
  filter(Curve %in% c("f2", "B-Spline", "Monomial", "Fourier"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/basis_expansion_2_1.pdf", plot = basis_expansion_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

basis_expansion_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
  filter(Curve %in% c("f2", "B-Spline", "Monomial", "Fourier"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/basis_expansion_2_2.pdf", plot = basis_expansion_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 2 #####
fpcr_nharm2_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
  filter(Curve %in% c("f1", "bspline_nharm2", "monomial_nharm2", "fourier_nharm2"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm2_1_1.pdf", plot = fpcr_nharm2_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm2", "monomial_nharm2", "fourier_nharm2"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm2_1_2.pdf", plot = fpcr_nharm2_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm2", "monomial_nharm2", "fourier_nharm2"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm2_2_1.pdf", plot = fpcr_nharm2_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm2_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm2", "monomial_nharm2", "fourier_nharm2"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm2_2_2.pdf", plot = fpcr_nharm2_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 3 #####
fpcr_nharm3_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm3", "monomial_nharm3", "fourier_nharm3"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm3_1_1.pdf", plot = fpcr_nharm3_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm3", "monomial_nharm3", "fourier_nharm3"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm3_1_2.pdf", plot = fpcr_nharm3_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm3", "monomial_nharm3", "fourier_nharm3"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm3_2_1.pdf", plot = fpcr_nharm3_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm3_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm3", "monomial_nharm3", "fourier_nharm3"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm3_2_2.pdf", plot = fpcr_nharm3_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)

##### generate plots for fpcr nharm = 3 #####
fpcr_nharm4_plot_1_1 <- ggplot(data = plot_tibble_1_1 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm4", "monomial_nharm4", "fourier_nharm4"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm4_1_1.pdf", plot = fpcr_nharm4_plot_1_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_1_2 <- ggplot(data = plot_tibble_1_2 %>%
                                 filter(Curve %in% c("f1", "bspline_nharm4", "monomial_nharm4", "fourier_nharm4"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm4_1_2.pdf", plot = fpcr_nharm4_plot_1_2,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_2_1 <- ggplot(data = plot_tibble_2_1 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm4", "monomial_nharm4", "fourier_nharm4"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm4_2_1.pdf", plot = fpcr_nharm4_plot_2_1,
  width = 20, height = 12, units = "in", dpi = 600
)

fpcr_nharm4_plot_2_2 <- ggplot(data = plot_tibble_2_2 %>%
                                 filter(Curve %in% c("f2", "bspline_nharm4", "monomial_nharm4", "fourier_nharm4"))) +
  geom_line(aes(x = x, y = y, col = Curve)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 48),
    legend.text = element_text(size = 40)
  ) +
  guides(color = guide_legend(override.aes = list(lwd = 5)))

ggsave(
  filename = "../Graphics/Curve_Estimates/fpcr_nharm4_2_2.pdf", plot = fpcr_nharm4_plot_2_2,
  width = 20, height = 12, units = "in", dpi = 600
)