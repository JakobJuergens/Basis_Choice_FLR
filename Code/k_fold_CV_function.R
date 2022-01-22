fRegress.CVk <- function(y, xfdlist, betalist, wt=NULL, CVobs=10,
                        returnMatrix=FALSE, ...)
{
  
  # FREGRESS.CV computes cross-validated error sum of squares
  # for scalar or functional responses.
  
  # This function is a modified version of the function by Jim Ramsay
  # Original function: https://github.com/cran/fda/blob/master/R/fRegress.CV.R
  # Modification was done to use kfold CV, with the default of 10 folds
  
  #  check the arguments
  
  argList  <- fRegressArgCheck(y, xfdlist, betalist, wt)
  
  yfdobj   <- argList$yfd
  xfdlist  <- argList$xfdlist
  betalist <- argList$betalist
  wt       <- argList$wt
  
  # extract dimensions of the data and analysis
  
  p <- length(xfdlist)
  N <- dim(xfdlist[[1]]$coef)[2]
  dataset_length <- length(y)
  
  #  branch to either scalar or functional dependent variable
  
  if (inherits(yfdobj, "numeric"))  {
    
    #  scalar dependent variable case
    
    yvec   <- yfdobj
    SSE.CV <- 0
    errfd  <- c()
    
    # check if dataset can be divided evenly by number of CV number
    if (dataset_length %% CVobs == 0)  {
      i <- dataset_length / CVobs
    } else {print("ERROR")}##############################################################  do smth
    start = 1
    stop = i
    # start CV
    for (m in 1:CVobs) {
      
      #  eliminate case (start:stop) from the weights
      wti <- wt[-(start:stop)]
      xfdlisti <- vector("list",p)
      for (j in 1:p) {
        xfdj          <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
          betafdParj <- betalist[[j]]
          betafdj    <- betafdParj$fd
          basisj     <- betafdj$basis
          betarangej <- basisj$rangeval
          conbasisj  <- create.constant.basis(betarangej)
          xfdj       <- fd(matrix(xfdj,1,N), conbasisj)
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-(start:stop)],1,N-1)
        else                    coefj <- as.matrix(coefj[,-(start:stop)])
        xfdlisti[[j]] <- fd(coefj,basisj)
      }
      yveci         <- yvec[-(start:stop)]
      fRegressListi <- fRegress(yveci, xfdlisti, betalist, wti)
      betaestlisti  <- fRegressListi$betaestlist
      yhati <- c()
      
  
      for (j in 1:p) {
        # use test data on trained model, get yhat for test data
        for (curve in start:stop){
          betafdParj <- betaestlisti[[j]]
          betafdj    <- betafdParj$fd
          xfdj       <- xfdlist[[j]]
          bbasisj    <- betafdj$basis
          rangej     <- bbasisj$rangeval
          nfine      <- max(501, bbasisj$nbasis*10+1)
          tfine      <- seq(rangej[1], rangej[2], len=nfine)
          delta      <- tfine[2]-tfine[1]
          betavec    <- eval.fd(tfine, betafdj, 0, returnMatrix)
          xveci      <- eval.fd(tfine, xfdj[curve], 0, returnMatrix)
          yhati     <- c(yhati, delta*(sum(xveci*betavec) -
                                         0.5*( xveci[1]    *betavec[1] +
                                                 xveci[nfine]*betavec[nfine] )))
        }
       
      }
      errfd[m] = mean(yvec[(start:stop)] - yhati);
      SSE.CV <- SSE.CV + errfd[m]^2
      start = start +i
      stop = stop +i
    }
  } 
  return(list(SSE.CV=SSE.CV, errfd.cv=errfd, fRegress_obj = fRegressListi))
}


# fRegressArgCheck is copied from the fda package, found here: https://github.com/cran/fda/blob/master/R/fRegressArgCheck.R
fRegressArgCheck <- function(yfd, xfdlist, betalist, wt=NULL) 
{
  #  FREGRESS_ARGCHECK checks the first four arguments for the functions
  #  for function regression, including FREGRESS.
  
  #  Last modified 16 December 2020 by Jim Ramsay
  
  #  --------------------  Check classes of arguments  --------------------
  
  #  check that YFD is of class either 'fd' or 'numeric' and compute sample size N
  
  if (!(is.fdPar(yfd) || is.fd(yfd) || is.numeric(yfd) || is.matrix(yfd))) stop(
    "First argument is not of class 'fdPar', 'fd', 'numeric' or 'matrix'.")
  
  #  As of 2020, if yfd is an fdPar object, it is converted to an fd object.
  #  The added structure of the fdPar class is not used in any of the fRegress codes.
  # The older versions of fda package used yfdPar as the name for the first member.
  
  if (is.fdPar(yfd)) yfd <- yfd$fd
  
  if (inherits(yfd, "fd")) {
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
  } else {
    N <- length(yfd)
  } 
  
  #  check that xfdlist is a list object and compute number of covariates p
  
  #  check XFDLIST
  
  if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)
  
  if (!inherits(xfdlist, "list")) stop(
    "Argument XFDLIST is not a list object.")
  
  #  get number of independent variables p
  
  p <- length(xfdlist)
  
  #  check BETALIST
  
  if (inherits(betalist, "fd")) betalist <- list(betalist)
  
  if (!inherits(betalist, "list")) stop(
    "Argument BETALIST is not a list object.")
  
  if (length(betalist) != p)  {
    cat(paste("\nNumber of regression coefficients does not match\n",
              "number of independent variables."))
    stop("")
  }
  
  #  extract the range if YFD is functional
  
  if (inherits(yfd, "fd")) {
    rangeval <- yfd$basis$rangeval
  } else {
    rangeval = c(0,1)
    #   allscalar <- TRUE
    #   for (j in 1:p) {
    #     if (inherits(xfdlist[[j]], "fd")) {
    #       rangeval <- xfdlist[[j]]$basis$rangeval            
    #       allscalar <- FALSE
    #       break
    #     }
    #   }
    # if (allscalar) stop(
    #   paste("The dependent variable and all the independent",   
    #         "variables are scalar."))
  }
  
  #  --------------------  check contents of XFDLIST  -------------------
  
  #  If the object is a vector of length N,
  #  it is converted to a functional data object with a
  #  constant basis
  
  onebasis <- create.constant.basis(rangeval)
  onesfd   <- fd(1,onebasis)
  
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xcoef <- xfdj$coefs
      if (length(dim(xcoef)) > 2) stop(
        paste("Covariate",j,"is not univariate."))
      #  check size of coefficient array
      Nj <- dim(xcoef)[2]
      if (Nj != N) {
        print(
          paste("Incorrect number of replications in XFDLIST",
                "for covariate",j))
        xerror = TRUE
      }
    } 
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != N && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
    } 
    if (!(inherits(xfdlist[[j]], "fd"     ) || 
          inherits(xfdlist[[j]], "numeric") ||
          inherits(xfdlist[[j]], "matrix" ))) {
      print(paste("XFDLIST[[",j,"]] is not an FD or numeric or matrix object."))
      xerror = TRUE
    }
  }
  
  #  --------------------  check contents of BETALIST  -------------------
  
  berror <- FALSE
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (inherits(betafdParj, "fd") || inherits(betafdParj, "basisfd")) {
      betafdParj    <- fdPar(betafdParj)
      betalist[[j]] <- betafdParj
    }
    if (!inherits(betafdParj, "fdPar")) {
      print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
      berror <- TRUE
    }
  }
  
  if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")
  
  #  --------------------  check contents of WEIGHTS  -------------------
  
  if (is.null(wt)) wt = rep(1,N)
  if (length(wt) != N) stop("Number of weights not equal to N.")
  if (any(wt < 0))     stop("Negative weights found.")
  
  #  ---------------------  return the argument list  --------------------
  
  # The older versions of fda package used yfdPar as the name for the first member.
  
  return(list(yfd=yfd, xfdlist=xfdlist, betalist=betalist, wt=wt))
  
}
############################

# This is the origional function from here: https://github.com/cran/fda/blob/master/R/fRegress.CV.R
# We use it for testing the modified function
# May be removed in the future
fRegress.CVoriginal <- function(y, xfdlist, betalist, wt=NULL, CVobs=1:N,
                                returnMatrix=FALSE, ...)
{
  
  # FREGRESS.CV computes cross-validated error sum of squares
  # for scalar or functional responses. NOTE: ordinary and
  # generalized cross validation scores are now returned by fRegress
  # when scalar responses are used.
  
  # last modified 16 December 2020 by Jim Ramsay
  
  #  check the arguments
  
  argList  <- fRegressArgCheck(y, xfdlist, betalist, wt)
  
  yfdobj   <- argList$yfd
  xfdlist  <- argList$xfdlist
  betalist <- argList$betalist
  wt       <- argList$wt
  
  # extract dimensions of the data and analysis
  
  p <- length(xfdlist)
  N <- dim(xfdlist[[1]]$coef)[2]
  M <- length(CVobs)
  
  #  branch to either scalar or functional dependent variable
  
  if (inherits(yfdobj, "numeric"))  {
    
    #  scalar dependent variable case
    
    yvec   <- yfdobj
    SSE.CV <- 0
    errfd  <- c()
    for (m in 1:M) {
      i        <- CVobs[m]
      #  eliminate case i from the weights
      wti <- wt[-i]
      xfdlisti <- vector("list",p)
      for (j in 1:p) {
        xfdj          <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
          betafdParj <- betalist[[j]]
          betafdj    <- betafdParj$fd
          basisj     <- betafdj$basis
          betarangej <- basisj$rangeval
          conbasisj  <- create.constant.basis(betarangej)
          xfdj       <- fd(matrix(xfdj,1,N), conbasisj)
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-i],1,N-1)
        else                    coefj <- as.matrix(coefj[,-i])
        xfdlisti[[j]] <- fd(coefj,basisj)
      }
      yveci         <- yvec[-i]
      fRegressListi <- fRegress(yveci, xfdlisti, betalist, wti)
      betaestlisti  <- fRegressListi$betaestlist
      yhati <- 0
      for (j in 1:p) {
        betafdParj <- betaestlisti[[j]]
        betafdj    <- betafdParj$fd
        xfdj       <- xfdlist[[j]]
        bbasisj    <- betafdj$basis
        rangej     <- bbasisj$rangeval
        nfine      <- max(501, bbasisj$nbasis*10+1)
        tfine      <- seq(rangej[1], rangej[2], len=nfine)
        delta      <- tfine[2]-tfine[1]
        betavec    <- eval.fd(tfine, betafdj, 0, returnMatrix)
        xveci      <- eval.fd(tfine, xfdj[i], 0, returnMatrix)
        yhati      <- yhati + delta*(sum(xveci*betavec) -
                                       0.5*( xveci[1]    *betavec[1] +
                                               xveci[nfine]*betavec[nfine] ))
      }
      errfd[i] = yvec[i] - yhati;
      SSE.CV <- SSE.CV + errfd[i]^2
      print(errfd)
     
    }
  } else {
    
    #  functional dependent variable case
    
    yfd      <- yfdobj
    SSE.CV   <- 0
    errcoefs <- c()
    for(m in 1:length(CVobs)){
      # index of case to eliminate
      i <-  CVobs[m]
      # eliminate case i from the weights
      wti <- wt[-i]
      # eliminate case i from covariates
      txfdlist <- xfdlist
      for(k in 1:p){
        txfdlist[[k]] <- xfdlist[[k]][-i]
      }
      # eliminate case i from dependent variable
      yfdi <- yfd[-i]
      # carry out the functional regression analysis
      tres <- fRegress(yfdi,txfdlist,betalist,wti)
      #  extract the regression coefficient functions
      betaestlisti <- tres$betaestlist
      #  compute the fit to the data for case i
      yhatfdi <- 0
      for(k in 1:p){
        betafdPark = betaestlisti[[k]]
        betafdk    = betafdPark$fd
        xfdk       = xfdlist[[k]]
        xfdik      = xfdk[i]
        tempfd     = xfdik*betafdk
        yhatfdi <- yhatfdi + tempfd
      }
      #  compute the residual function
      errfdi   <- yfd[i] - yhatfdi
      #  increment the error sum of squares by the integral of the
      #  square of the residual function
      SSE.CV   <- SSE.CV + inprod(errfdi,errfdi)
      #  add the coefficients for the residual function
      errcoefs <- cbind(errcoefs,errfdi$coefs)
    }
    #  set up the functional data object for the residual fns
    errfd <- fd(errcoefs,errfdi$basis)
    names(errfd$fdnames)[[3]] <- "Xval Errors"
  }
  return(list(SSE.CV=SSE.CV, errfd.cv=errfd))
}


