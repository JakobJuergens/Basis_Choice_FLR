rm(list=ls())
data_generation <-function(fun){
    var1 = 1
    var2 = 2
    var3 = 0.2

    X <- matrix(0, nrow=n_obs, ncol=length(grid))
    for(i2 in 1:n_obs){
        X[i2,]=X[i2,]+rnorm(length(grid), 2, var1)
        X[i2,]=X[i2,]+runif(1, 0, var2)
        X[i2,]=X[i2,]+rnorm(1, 0, var3)*grid
        X[i2,]=X[i2,]*fun

        for(j2 in 1:5){
           e =abs(rnorm(2, 0, var1/j2^(2)))
           X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
           X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
        }
    }
    return(X)
}

data(gasoline)
octane <- (gasoline$octane)
NIR    <- as.matrix(gasoline$NIR)
###################################################################################################################
# set up "global" variables
###################################################################################################################
set.seed(100)
n_obs = 60
nharm = 4
n_var = 400
grid = seq(0, 1, length = n_var+1)
#smooth
f1 <- 2*sin(0.5*pi*grid) + 4*sin(1.5*pi*grid) + 5*sin(2.5*pi*grid)
#bumpy
f2 <- 1.5*exp(-0.5*(grid-0.3)^2/0.02^2) -  4*exp(-0.5*(grid-0.45)^2/0.015^2) +  8*exp(-0.5*(grid-0.6)^2/0.02^2) -  exp(-0.5*(grid-0.8)^2/0.03^2)
#two different variances of error
sigma_eps_squared1_1 = as.numeric((var(NIR %*% f1)/0.9) - var(NIR %*% f1))
sigma_eps_squared1_2 = as.numeric((var(NIR %*% f1)/0.6) - var(NIR %*% f1) )
sigma_eps_squared2_1 = as.numeric((var(NIR %*% f2)/0.9) - var(NIR %*% f2))
sigma_eps_squared2_2 = as.numeric((var(NIR %*% f2)/0.6) - var(NIR %*% f2) )



bspline_function <- function(rep, NIR, n_obs){
    CV_container_spline  <- c()
        for(j in seq(11,12,1)){

            CV_container  <- matrix(NaN, nrow = rep, ncol = 4)

            for(i in 1 : rep){

                #########
                # Divide in test / training!
                #########
                #each true beta and variance
                Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
                Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
                Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
                Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
        
                data = t(NIR)

                #print(dim(NIR))
                smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = 4)
                smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd
                xfdlist = list(smooth_basis=smooth_basis_fd)
                betabasis1 <- create.constant.basis(c(0, 60))
                betafd1    <- fd(0, betabasis1)
                betafdPar1 <- fdPar(betafd1)
                #betafd2    <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = as.numeric(order) )
                betafdPar2  <- fdPar(smallbasis)
                betalist <- list(smooth_basis_fd_data=betafdPar2)
    
                f_regress1_1 <- fRegress.CV(y = Y1_1, xfdlist, betalist)
                f_regress1_2 <- fRegress.CV(y = Y1_2, xfdlist, betalist)
                f_regress2_1 <- fRegress.CV(y = Y2_1, xfdlist, betalist)
                f_regress2_2 <- fRegress.CV(y = Y2_2, xfdlist, betalist)
        
                CV_container[i,1] <- f_regress1_1$SSE.CV
                CV_container[i,2] <- f_regress1_2$SSE.CV
                CV_container[i,3] <- f_regress2_1$SSE.CV
                CV_container[i,4] <- f_regress2_2$SSE.CV

            }
            scaled_MSE <- colMeans(CV_container)
            scaled_MSE[5] = j
            scaled_MSE[6] = order
            CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
        
    }
    colnames(CV_container_spline) = c("f1_e1_spline", "f1_e2_spline", "f2_e1_spline", "f2_e2_spline", "n_basis")
    return(CV_container_spline)
}

fpcr_function <- function(rep, NIR, n_obs){
    CV_container_spline  <- c()
    train.control <- trainControl(method = "cv", number = 10)
        for(j in seq(11,12,1)){

            CV_container  <- matrix(NaN, nrow = rep, ncol = 5)

            for(i in 1 : rep){

                #########
                # Divide in test / training!
                #########
                #each true beta and variance
                Y1_1 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_1)
                Y1_2 <- as.numeric(NIR %*% f1 + rnorm(n_obs, 0, 1) * sigma_eps_squared1_2)
                Y2_1 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_1)
                Y2_2 <- as.numeric(NIR %*% f2 + rnorm(n_obs, 0, 1) * sigma_eps_squared2_2)
        
                data = t(NIR)
                
                #print(dim(NIR))
                smallbasis      <- create.bspline.basis(rangeval = c(0, length(grid)), nbasis = as.numeric(j), norder = 4)
                smooth_basis_fd <- smooth.basis(y = data, fdParobj=smallbasis)$fd

                simulated_pcaObj = pca.fd(smooth_basis_fd, nharm = nharm, centerfns = TRUE)

         
                dataframe1_1 = data.frame(Y1_1, simulated_pcaObj$scores)
                CV_container[i,1] <- train(Y1_1 ~., data = dataframe1_1, method = "lm", trControl = train.control)$results[[2]]

                dataframe1_2 = data.frame(Y1_2, simulated_pcaObj$scores)
                CV_container[i,2] <- train(Y1_2 ~., data = dataframe1_2, method = "lm", trControl = train.control)$results[[2]]
                
                dataframe2_1 = data.frame(Y2_1, simulated_pcaObj$scores)
                CV_container[i,3] <- train(Y2_1 ~., data = dataframe2_1, method = "lm", trControl = train.control)$results[[2]]

                dataframe2_2 = data.frame(Y2_2, simulated_pcaObj$scores)
                CV_container[i,4] <- train(Y2_2 ~., data = dataframe2_2, method = "lm", trControl = train.control)$results[[2]]

                CV_container[i,5] <- sum(simulated_pcaObj$varprop)

            }
            scaled_MSE <- colMeans(CV_container)
            scaled_MSE[6] = j
            scaled_MSE[7] = order
            CV_container_spline <- rbind(CV_container_spline, scaled_MSE)
        
    }
    colnames(CV_container_spline) = c("f1_e1_fpcr", "f1_e2_fpcr", "f2_e1_fpcr", "f2_e2_fpcr","varprop" ,"n_basis")
    return(CV_container_spline)
}



################
# AIC and Mallows CP are not as expected!
# When using GCV / 10fcv, I got results as in validation set approach, so maybe we should stick with that!
################
aic_calculation = function(fregress_obj, params, gridlength){
#' Computes AIC for a fRegress object
#' Inputs: 
#' fregress_obj (obj): Object of class fda::fRegress, 
#' params (int): number of parameters of the model, 
#' gridlength (int): greidlength on which the function is evaluated
    mu_ml = mean(fregress_obj$yfdobj - fregress_obj$yhatfdobj)
    print(mu_ml)
    sig_sq = (1/length(fregress_obj$yfdobj))*sum(((fregress_obj$yfdobj - fregress_obj$yhatfdobj) - mu_ml)^2)
    #sig_sq = mean((f_regress_obj$yfdobj - mean(f_regress_obj$yhatfdobj))^2)
    AIC = 0
    for(i in 1:length(fregress_obj$yfdobj)){
        AIC = AIC -((gridlength/2)*log(2*pi) -(gridlength/2)*log(sig_sq) 
            -(1/2*sig_sq)*t(fregress_obj$yfdobj[i] - fregress_obj$yhatfdobj[i]) %*% 
            (fregress_obj$yfdobj[i] - fregress_obj$yhatfdobj[i]))
    }
    return(-AIC + params)
    }
#aic_calculation(f_regress, 15, 150)

mellow_cp = function(fregress_obj, fregress_obj_max_p, n_obs, p){
#' Computes MellowCP for a fRegress object
#' Inputs: 
#' fregress_obj (obj): Object to evaluate of class fda::fRegress, 
#' fregress_obj_max_p (obj): Object constructed using max number, default = 40, of bsplines
#' nobs (int): number of observations 
#' p (int): number of splines used to construct fregress_obj

    rsqr = sum((fregress_obj$yfdobj - fregress_obj$yhatfdobj)^2)
    mse40splines = mean((fregress_obj_max_p$yfdobj - fregress_obj_max_p$yhatfdobj)^2) #mse with max K

    MellowCP = (rsqr/mse40splines) - n_obs +2*(p+1)
    return(MellowCP)
}
#mellow_cp(f_regress, f_regress2, n_obs, p)