# Clear Workspace
rm(list = ls())

# load libraries
library(fda)
library(refund)
library(tidyverse)

# load data
data(gasoline)
NIR    <- t(as.matrix(gasoline$NIR))

wavelengths <- seq(from = 900, to = 1700, by = 2)

NIRbasis5 <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                   nbasis = 5, norder = 5)
NIRbasis10 <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                  nbasis = 10, norder = 5)
NIRbasis30 <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                  nbasis = 30, norder = 5)
NIRbasis50 <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                  nbasis = 50, norder = 5)
NIRbasis100 <- create.bspline.basis(rangeval = c(0, max(wavelengths)),
                                  nbasis = 100, norder = 5)
NIRbasis5 <- fdPar(NIRbasis5)
NIRbasis10 <- fdPar(NIRbasis10)
NIRbasis30 <- fdPar(NIRbasis30)
NIRbasis50 <- fdPar(NIRbasis50)
NIRbasis100 <- fdPar(NIRbasis100)


gasoline_fd5 <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                            y = NIR, fdParobj = NIRbasis5)$fd
gasoline_fd10 <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                            y = NIR, fdParobj = NIRbasis10)$fd
gasoline_fd30 <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                            y = NIR, fdParobj = NIRbasis30)$fd
gasoline_fd50 <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                            y = NIR, fdParobj = NIRbasis50)$fd
gasoline_fd100 <- smooth.basis(argvals = seq(0, to = 1700, length.out = 401),
                            y = NIR, fdParobj = NIRbasis100)$fd

observation5 <- as.matrix(eval.fd(seq(0,1700, len = 401), gasoline_fd5, int2Lfd(0))[,1])
observation10 <- as.matrix(eval.fd(seq(0,1700, len = 401), gasoline_fd10, int2Lfd(0))[,1])
observation30 <- as.matrix(eval.fd(seq(0,1700, len = 401), gasoline_fd30, int2Lfd(0))[,1])
observation50 <- as.matrix(eval.fd(seq(0,1700, len = 401), gasoline_fd50, int2Lfd(0))[,1])
observation100 <- as.matrix(eval.fd(seq(0,1700, len = 401), gasoline_fd100, int2Lfd(0))[,1])

NIR_seq <- as.matrix(seq(900, 1700, length = 401))

NIR_plot5 <- cbind(observation5, NIR_seq)
NIR_plot10 <- cbind(observation10, NIR_seq)
NIR_plot30 <- cbind(observation30, NIR_seq)
NIR_plot50 <- cbind(observation50, NIR_seq)
NIR_plot100 <- cbind(observation100, NIR_seq)

colnames(NIR_plot5) = c('wavelength','time')
colnames(NIR_plot10) = c('wavelength','time')
colnames(NIR_plot30) = c('wavelength','time')
colnames(NIR_plot50) = c('wavelength','time')
colnames(NIR_plot100) = c('wavelength','time')

NIR_plot5 <- as.data.frame(NIR_plot5)
NIR_plot10 <- as.data.frame(NIR_plot10)
NIR_plot30 <- as.data.frame(NIR_plot30)
NIR_plot50 <- as.data.frame(NIR_plot50)
NIR_plot100 <- as.data.frame(NIR_plot100)

NIR <- cbind(as.data.frame(NIR[,1]),NIR_seq)
colnames(NIR) = c('wavelength','time')
plot(NIR[,wavelength])

NIR$'K' = 'K=0'
NIR_plot5$'K' = 'K=5'
NIR_plot10$'K' = 'K=10'
NIR_plot30$'K' = 'K=30'
NIR_plot50$'K' = 'K=50'
NIR_plot100$'K' = 'K=100'

NIR_data = rbind(NIR,NIR_plot5,NIR_plot10,NIR_plot30,NIR_plot50,NIR_plot100)

label <- paste0("K=",c(0,5,10,30,50,100))
label[1] <- "Real Line"

basis_exp_plot <- ggplot(data = NIR_data, aes(x=time, y=wavelength,color=K)) +
  geom_line()+
  xlab('Wavelength in nm') +
  ylab('Spectrum') +
  scale_color_discrete(breaks=c('K=0','K=5','K=10','K=30','K=50','K=100'), 
                       labels=label)+
  ggtitle('Basis expansions according to the different number of B-spline basis functions') +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 30),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20)
  )

# save plot in appropriate folder
ggsave(
  filename = "../Graphics/basis_expansions.pdf", plot = basis_exp_plot,
  width = 20, height = 12, units = "in", dpi = 600
)
  