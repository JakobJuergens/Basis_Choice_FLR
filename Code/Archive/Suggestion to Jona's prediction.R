rm(list=ls())

library(refund);
library(ggplot2) ;
library(dplyr); 
library(reshape2); 
library(fda);

data(gasoline)

octane = (gasoline$octane)
waves = t(as.matrix(gasoline$NIR))

test = seq(1,30,1)

waves_test= waves[,test]
octane_test = octane[test]

waves_train= waves[,-test]
wavelengths = 1:401
octane_train = octane[-test]

smallbasis  <- create.bspline.basis(rangeval = c(0, 401), nbasis = 20, 5)

gasoline_fd_train <- smooth.basis(wavelengths, waves_train, smallbasis)$fd
gasoline_fd_test <- smooth.basis(wavelengths , waves_test, smallbasis)

f_regress2 <- fRegress(octane_train ~ gasoline_fd_train)

### I guess gasoline_fd_test$betaestlist[[2]] calls coefficients of basis functions.
### Does it, then, work only with coefficients without basis functions
### (maybe it doesn't matter because the same basis functions are shared and they are decided by coefficients)?
prec <- predict(f_regress2, gasoline_fd_test$betaestlist[[2]]$fd)

plot(prec,octane_test)
MSE <- mean((octane_test - prec)^2)
print(MSE)
abline(a=0,b=1,col="black",lty=6)
