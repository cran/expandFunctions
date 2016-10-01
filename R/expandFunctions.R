# expandFunctions.R
# roxygen2::roxygenize(getwd()) # save file first!

#' Reset annoyingly persistent warning messages.
#' @return Returns TRUE invisibly.  Used for side effects only.
#' @export
#' @references
#' This function is built around the snippet found here:
#' \url{http://stackoverflow.com/questions/5725106/r-how-to-clear-all-warnings}
#' @examples
#' ## reset.warnings()
# This might be better as an add-in
reset.warnings <- function() {
  assign("last.warning", NULL, envir = baseenv())
  invisible(TRUE)
}

#' Select and fit sparse linear model with LASSO
#'
#' @description The purpose of this function is to make the process
#'              of LASSO modelling as simple as possible.
#'
#'              This is a simple wrapper on two glmnet functions
#'              which, when given input matrix X and response vector
#'              y, and a criterion for model selection, will
#'              estimate the lambda parameter, and return the
#'              LASSO results as a glmnet model.  This model
#'              can then be used to find coefficients and predictions.
#'
#' @param X         Predictor matrix, nXp, with n observations and p
#'                  features.
#' @param y         Response vector, or column or row matrix.  Must
#'                  have length n.
#' @param criterion String describing which lambda criterion to
#'                  use in selecting a LASSO model.  Choices
#'                  currently are c("lambda.1se","lambda.min").
#' @return a glmnet model
#' @export
#' @seealso \code{\link[glmnet]{glmnet}} and
#'          \code{\link[glmnet]{cv.glmnet}}
#' @examples
#' set.seed(1)
#' nObs <- 100
#' X <- distMat(nObs,6)
#' A <- cbind(c(1,0,-1,rep(0,3)))
#'   # Y will only depend on X[,1] and X[,3]
#' Y <- X %*% A + 0.1*rnorm(nObs)
#' lassoObj <- easyLASSO(X=X,y=Y) # LASSO fitting
#' Yhat <- predict(lassoObj,newx=X)
#' yyHatPlot(Y,Yhat)
#' coef( lassoObj ) # Sparse coefficients
#' coefPlot( lassoObj )
easyLASSO <- function(X,y,criterion="lambda.1se") {

  # Cross-validate a choice of lambda
  cvglmnetObj <- glmnet::cv.glmnet(x=X,y=y)

  # Fit LASSO using the estimated lambda and criterion
  lassoObj <- glmnet::glmnet(x=X,y=y,lambda=cvglmnetObj[[criterion]])

  return( lassoObj )
}

#' @title Extend outer product.
#'
#' @description
#' Extends outer \{base\} \code{outer(x,y,FUN)} to include functions
#' \code{FUN(x,y,\dots)} where the first argument of \code{FUN}
#' is a vector but the second argument must be a scalar.
#'
#' @param x     Vector, with the same function as
#'              in outer \{base\}.  Each value will
#'              correspond to a row in the return matrix.
#' @param y     Vector.  Each
#'              element in the vector corresponds
#'              to a column in the return matrix.
#' @param FUN   Function. x and y will
#'              be the first and second arguments.  Unlike
#'              \code{outer}, however, while a vector can be the
#'              first argument, FUN might only allow
#'               \emph{one value} as the
#'              second argument.  This means eOuter
#'              can use lagshift, for instance, as FUN.
#' @param ...   Additional parameters for FUN.
#'
#' @details
#' outer has limitations; it only works with functions which
#' can take vector inputs for \emph{both}
#' the first and second arguments, such as "^".  As a result,
#' many functions cannot be used for FUN.  The function eOuter
#' gets around this limitation by additionally allowing functions
#' which accept a vector for the first argument, but only scalars
#' for the second argument.  It can be used everywhere
#' that \code{outer} can be used, but also when FUN is
#' limited in this way.
#'
#' @return A matrix \code{Z} of size
#' \code{length(x) X length(y)}
#' containing \code{Z[,i]} with values \code{FUN(x,y[i],\dots)}.
#' @export
#' @seealso \code{\link{outer}} and \code{\link{ePow}}
#' @examples
#' # This implements a function similar to ePow
#' # FIXME: change ePow to use eOuter!!!
#' eOuter(1:6,0:2,FUN = `^`)
#' # Other functions of columns
#' eOuter(1:10,0:3,FUN = lagshift,lagMax=3,pad=NA)
#' # FIXME: Make function to allow polynomials to be used:
#' # eOuter(1:10,1:3,FUN = glaguerre.polynomials, alpha=0.5)
eOuter <- function(x,y,FUN,...) {
  temp <- lapply(y,function(w) FUN(x,w,...))
  return( do.call(cbind, temp ) )
}

#' Convert vector x into a matrix \eqn{X_{ij} = {x_i}^j}
#'
#' @param x              Data vector.
#' @param colParamVector Vector of column powers.
#'
#' @return A matrix X of size length(x) X length(colParamVector)
#'         \deqn{X_{ij} = {x_i}^j}
#' @export
#' @examples
#' x <- 1:6
#' ePow(x,0:2)
ePow <- function(x,colParamVector) {
  return( eOuter(x,colParamVector,FUN=`^`) )
}

#' Helper function for eLag.
#'
#' Generates shifted columns.
#' @param x       Input vector
#' @param i       Shift (integer)
#' @param lagMax  Maximum lag that will be needed
#' @param pad     Scalar used for padding.
#' @return vector padded front and back with padding appropriate
#' for generating lag.
#' @export
#' @examples
#' lagshift(1:3,0,1,NA)
lagshift <- function(x,i,lagMax,pad) {
  return( c(rep(pad,i),x,rep(pad,lagMax - i)) )
}

#' Convert vector into a matrix of lag columns
#'
#' @param x               Data \emph{vector}
#' @param colParamVector  Vector of lags for embedding
#' @param pad             Scalar for padding embedding
#' @return A matrix whose columns are x lagged by the
#'         corresponding values in colParamVector.
#' @export
#' @seealso \code{\link[stats]{embed}} and
#' \code{\link[tseriesChaos]{embedd}}, which
#'  are related functions.
#' @examples
#' eLag(1:6, 0:2)
#' eLag(1:6, 0:2, pad=0)
eLag <- function(x,colParamVector,pad=NA) {
  colParamVector <- colParamVector - min(colParamVector)
  lagMax <- max(colParamVector)
  X <- eOuter(x,colParamVector,FUN=lagshift,
              lagMax=lagMax,pad=pad)
  return( X )
}

#' Multivariate second order polynomial expansion.
#'
#' Expand matrix columns into linear, square, and unique product columns.
#'
#' @param X   vector or matrix.  If a vector, it will be converted to
#'       a column matrix.  If it is desired that the squares
#'       and products of a \emph{vector} are computed, pass rbind(X)
#'       instead of X, and thereby pass a row matrix.
#' @param FUN Binary function which forms the products of the columns.
#'       By default, this is `*`, but other \emph{commuting} operators
#'       or kernels can be used if desired.
#' @param ... Options for FUN.  Not needed if FUN doesn't have options.
#'
#' @details Form a matrix with columns composed of  into linear, square, and
#' product columns:
#'
#' \deqn{[X | FUN(X[,i], X[,j])]}
#'
#' where \eqn{i, j} are the unique combinations of \eqn{i} and \eqn{j},
#' including \eqn{i=j}.
#'
#' By default, the function used to form the squares and
#' products, FUN, is just conventional multiplication = `*`, but any
#' \emph{commuting} binary operator can be used.
#'
#' This particular expansion is often applied in
#'
#' \itemize{
#'   \item{General Method of Data Handling (GMDH).}
#'   \item{Nonlinear Slow Feature Analysis (SFA).  Performing
#'   a multivariate polynomial of second degree expansion
#'   in all the features, then
#'   performing \emph{linear} SFA on the resulting expanded
#'   feature matrix, is a very common approach, and in fact
#'   is the default method in \code{sfa2 \{rSFA\}}.}
#' }
#'
#' @return \eqn{[X,X^2,unique products of columns of X]}.  The unique
#'       products are in row major upper right triangular order.
#'       Thus, for X with columns 1:3, the order is
#'       \deqn{X[,1]^2, X[,2]^2, X[,3]^2,
#'       X[,1]*X[,2], X[,1]*X[,3], X[,2]*X[,3]}
#' @export
#' @seealso \code{\link[rSFA]{sfa2}}
#' @examples
#' # # Examples
#' # eQuad(1:5)
#' # eQuad(matrix(1:12,ncol=3),FUN=`+`)
eQuad <- function(X,FUN=`*`,...) {

  # If X is a vector, make a column matrix
  if (is.vector(X)) X <- cbind(X)

  nVar <- NCOL(X) # Number of columns

  # If X only has one column, there are no product columns,
  # so handle as a special case
  if (nVar==1) return( cbind(X,FUN(X,X,...)) )

  a <- combn(1:nVar,2) # Matrix of all unique (transposes excluded)
  # 2 element combinations of 1:nVar.  Each row is one
  # pair

  # Form X, X^2 and unique products of columns of X and return
  return( cbind( X, FUN(X,X,...), FUN(X[,a[1,]], X[,a[2,]], ...) ) )
}

# wrap polynomial functions into a format of
# FUN(x,degree,...); this will allow eOuter to be used
# to build matrices from vectors.
# Additionally, it would be good to have a "matricize"
# function that allows this to be applied to an X matrix
# rather than an x vector, generating blocks corresponding
# to each column.

#' Extends eOuter to allow a matrix for the first argument
#'
#' @param X              R object coercible to a matrix
#'                       the columns of this will be the
#'                       argument of FUN (see below).
#' @param colParamVector Vector input which will be the second
#'                       argument of FUN (see below).
#' @param FUN            Function which will be applied to
#'                       FUN(X[,i],colParamVector[j],...)
#' @param ...            Additional arguments to FUN.
#' @details This function is a simple extension of eOuter which allows
#' the function eOuter(X[,i],colParamVector,FUN,...) for
#' i in the columns of X.
#' @return  Returns a matrix with the matrics generated by eOuter
#' for each column column bound together.  This means that each
#' row of the returned matrix represents single observations (at
#' least as long as no lags are used).
#' @export
#' @examples
#' A <- matrix(1:6,ncol=2)
#' temp <- eMatrixOuter(A,0:2,FUN=`^`)
eMatrixOuter <- function(X,colParamVector,FUN,...) {
  return(
    do.call(cbind,
      alply(X,2,function(y) eOuter(y,colParamVector,FUN,...)))
  )
}

#' Matrix size-preserving diff function
#'
#' Returns a matrix, the same size as the
#' input matrix X, containing the back difference.
#' @param X     R object coercible to matrix
#' @param pad   Pad the first row with this value;
#'              the default is NA. 0 would be another value often
#'              used in signal processing.
#' @return Returns a matrix, the same size as the
#'         input matrix X, containing the back difference by column.
#'         The first row is filled with copies of pad.
#' @export
#' @examples
#' eDiff( 1:8 )
#' eDiff( as.data.frame(1:8) )
#' eDiff( matrix(1:8,ncol=2) )
# The pattern used here makes sense for any operator
# that operates by columns.
eDiff <- function(X,pad=NA) {
  temp <- diff( as.matrix(X) )
  return( rbind(rep(pad,NCOL(temp)), temp ) )
}

#' Remove padded rows from matrix X
#'
#' @param X          R object coercible to matrix
#' @param pad        Value representing padded elements.  By
#'                   default it is NA, but could be any value.
#' @return  A matrix.
#' @export
#' @examples
#' n <- 10
#' x <- rnorm(n)    # x vector
#' X <- eLag(x,0:1) # X matrix
#' t <- 1:n         # time vector
#' T <- eLag(t,0:1) # time matrix; the column corresponding
#'                  # to 0 is the time for each row,
#'                  # even after trimming
#' matplot(X,type="l",lty=1)
#' X <- eTrim(X)
#' T <- eTrim(T)
#' matplot(x=T[,1],y=X,type="l",lty=1,
#'   xlab="Time")
#
# The reference seems clunky; it would be nice if the
# reference information was carried along from operation to
# operation (say, as an attribute of the matrix object)
# TODO: Make this work for time series, and have the
# time series time values inherited.
eTrim <- function(X,pad=NA) {

  # Coerce X to a matrix; vectors become column matrices.
  X <- as.matrix(X)

  # Find and keep rows w/o pad
  cleanRow <- apply(X,1,function(y) !any(pad %in% y))
  X <- X[cleanRow,,drop=FALSE]

  return(X)
}

#' Replace values in an R object coerible to a matrix
#'
#' @description Replace values in an R object coerible to a matrix.
#'   It is useful for replacing NA with other values, etc.,
#'   such as with padding values.
#'
#' @param X R object coercible to a matrix
#' @param a Value to be replaced
#' @param b	Value to replace
#' @return X with all a's replaced with b's.  a may be NA
#' @export
#' @seealso \code{\link[base]{replace}}, which performs the same
#' operation on vectors, and on which this operation is based.
#' @examples
#' A <- matrix(1:6,ncol=2)
#' A <- eReplace(A,1,NA)
#' A <- eReplace(A,NA,-9999)
#' A <- eReplace(A,-9999,0)
eReplace <- function(X,a,b) {
  if ((length(a)!=1) || (length(b)!=1)) {
    stop("Use scalar values for a and/or b")
  }
  X <- as.matrix(X)

  if (is.na(a)) {
    X[is.na(X)] <- b
  } else {
    X[X==a]     <- b
  }
  return(X)
}

# Notes:
# purrr may be useful.  It is based on data.frames, not matrices,
# but much of the functionality is useful for both, or
# a matrix version easily made.  Most of purrr is actually
# oriented around manipulating lists, and making
# lists easier to use.
# purrr::prepend, append look useful.
# purrr::rbernouli - have written elsewhere.
# reduce, map, etc. - generally useful.
# split_by
# transpose for lists

# ExpandList
# 1) Uses a list of closure?  Function + parameter vector?
#    pipes w/data in function calls?  For optimization
#    needs to have parameters some place that optimization
#    routines can get at them, or has a function for extracting
#    the parameters, modifying them, and reinserting them.
# stft (with NA padding) is another expanding function
# random projection is another expanding closure (need
# to store the random projection matrix)
# recursive random projection (with NA padding)
# Jacobian blocks
# Layers (which is nearly there based on lists of expansions)
# Transformations - Gaussianization, Haar-Fisz, returns.
# Several "activation functions" can be used in expanding.
# AR values
# column differences
# last error(s)
# Blocking - takes vector and converts to frame-type matrix.
#   Note this is different from regular embedding, in that
#   the time offset can be different from row to row.
# Find time index of a given embedding, blocking, or descendent.

# Random projectors - store seed in attr? - probably
# not, since this may reset the seed for other routines.

# FIXME: make wrapper functions for polynomial functions
# to return simple functions of the right degree.
# (Tchebyshev) polynomial expansion
# arbitrary embedding would be good in this package as well.
# random projection would be good as well.
# polynom
# orthopolynom
# MonoPoly
# moments - might be useful
# mpoly - symbolic and more multivariate polynomials

#' Generate special functions using orthonormal functions
#'
#' orthopolynom can be used to generate special functions,
#' but for expansion they should be modified.  As of this
#' writing, orthopolynom generates polynomials for
#' Chebyshev, Hermite, Legendre and many other functions,
#' their integrals and derivatives, and more.
#'
#' The function polywrapper does 2 things:
#'
#' \itemize{
#'   \item{Generate functions from polynomial coefficients.}
#'   \item{Uses x as the 1st argument, and the order as
#'   the second; this means the generated functions can be
#'   used in eOuter and eMatrixOuter.}
#' }
#'
#' The functions so generated can be used as simple special functions,
#' as well as being useful in feature building.
#'
#' @param basePoly A polynomial list from orthopoly
#' @param kMax     Integer.  The maximum order of the function
#'                 generated.
#'
#' @details Since the coefficients from orthopolynom are generated
#'          by recursion, an upper limit of the function order
#'          needs to be set when calling polywrapper.  This is the
#'          main limitation of polywrapper.  Fortunately, since
#'          the functions are compactly stored, kMax can be set
#'          quite high if desired.  Note that usually the kMax
#'          is known, and is relatively small.
#'
#'          NB: The input x may need to be normalized.  orthopolynom
#'          has the function scaleX for just such a purpose.
#'
#' @return Function which is compatible with eOuter and eMatrixOuter
#' @export
#' @examples
#' # Generate a Chebyshev function of the form
#' # chebyFUN(x,k), where x is the input and k is the order.
#' # In this case, k must be no more than 5 (since that
#' # is the value passed to kMax), although it is
#' # easy to set this to a higher order if desired.
#' chebyFUN <- polywrapper(basePoly=orthopolynom::chebyshev.t.polynomials,
#'   kMax=5)
#' # Now the function chebyFUN
#' # can be used as any other function:
#' x <- seq(-1,+1,0.01)
#' plot(x,chebyFUN(x,5),type="l")
#' eOuter(seq(-1,+1,0.01),0:3,chebyFUN)
polywrapper <-
  function(basePoly=orthopolynom::chebyshev.t.polynomials,kMax=0) {

  # Get a list of the polynomials; this will
  # be encapsulated into the function returned
  polyList <- basePoly(kMax)

  # Build a function whose 2nd argument chooses the
  # function of the appropriate degree
  FUN <- function(x,k) {
    if (k>kMax) stop("Order must be no greater than kMax.")
    return( polynomial.values( list(polyList[[k+1]]), x )[[1]] )
  }

  return(FUN)
}

#' Make a matrix with coefficients distributed as dist
#'
#' @description Generate a pXq matrix with coefficients drawn
#'              from the univariate distribution dist with
#'              options distOpt.
#'
#' @param p         Number of rows
#' @param q         Number of columns
#' @param dist      distribution of coefficients to draw from;
#'                  default is rnorm.
#' @param distOpt   Named list of additional parameters for dist.
#'                  \emph{Always omit the first parameter,n, of the
#'                  distribution sampling function}. Defaults may
#'                  be omitted if desired (see examples).
#' @details The user may provide their own distribution function,
#'          but note that the number of values to return, n,
#'          \emph{must} be the first argument, just as with
#'          the built-in distributions.  The first argument does
#'          not have to be named.
#' @return A pXq matrix with coefficients distributed
#'         as dist with parameters defined in `...`
#' @export
#' @examples
#' X <- distMat(10,2)
#' X <- distMat(10,2,distOpt=list(mean=1,sd=2))
#' X <- distMat(5,3,rnorm,list(mean=1,sd=2))
#' X <- distMat(5,3,rnorm,list(sd=2))
#' X <- distMat(50,3,rt,list(df=3))
distMat <- function(p,q,dist=rnorm,distOpt=NULL) {
  w <- do.call(dist,c(p*q,distOpt,list()))
  return( matrix(w,nrow=p,ncol=q) )
}

#' @title Define a Random Affine Projection Transformation (RAPT) object
#'
#' @description
#' Create a Random Affine Projection Transformation (RAPT) object.
#' Such objects use random affine projection transformation to the
#' resulting matrix.  This allows RAPT objects serve as a basis
#' for a large number of kinds of expansions.
#'
#' @param   p        Number of input features (columns of \code{X}).
#' @param   q        Number of output features,
#' @param   Wdist    W distribution function.  Coefficients for
#'                   the random projection matrix W are drawn from
#'                   this distribution.  The default is rnorm.
#' @param   WdistOpt List of optional parameters for Wdist.
#'                   If this is NULL (default),
#'                   then only defaults of the distribution
#'                   are used.
#' @param   bDist    b distribution function.  Coefficients for
#'                   the offset vector b are drawn from this
#'                   distribution.   The default is runif.
#' @param   bDistOpt List of optional parameters for bDist.
#'                   If this is NULL
#'                   then only defaults of the distribution
#'                   are used.  The default is
#'                   \code{bDistOpt=list(min=0,max=0)},
#'                   which results in b = 0, with no offset.
#' @details This initializes a eRAPTobj, which holds all the
#'          parameters needed to perform a random projection
#'          transformation expansion (RAPT).
#'
#'          An RAPT of X is defined as
#'
#'          \deqn{X W + b}
#'
#'          where
#'
#'          X is the input matrix
#'
#'          W is a matrix of coefficients drawn from Wdist with
#'            options WdistOpt
#'
#'          b is a column matrix of coefficients drawn from bDist
#'            with options bDistOpt
#'
#'          If there is a need for multiple W or b distributions,
#'          then make multiple raptObj.  This makes
#'          each raptObj fairly simple, while allowing arbitrary
#'          complexity through multiple expansion and composition.
#'
#'          A simple way to get a linear feature, in addition
#'          to the RAPT features, is to simply cbind the
#'          original matrix X in with the raptObj matrix.
#'
#' @return An expand object, which defines the following fields:
#'   W       Input weighting matrix
#'   b       Input offset matrix
#' @export
#' @examples
#' raptObj <- raptMake(21,210,bDistOpt=list(min=-pi,max=+pi))
raptMake <- function(p, q,
                     Wdist=rnorm, WdistOpt=NULL,
                     bDist=runif, bDistOpt=list(min=0,max=0) ) {

  # Generate affine transform coefficients.
  W <- distMat(p,q,Wdist,WdistOpt)
  b <- distMat(q,1,bDist,bDistOpt)

  # Assemble components into eRAPTobj structure
  raptObj <- list(W = W, b = b)

  # Return object
  return( raptObj )
}

#' @title Expand an input matrix X using raptObj.
#'
#' @description Expand an input matrix X using
#' a Random Affine Projection Transformation (RAPT) object.
#' Such objects use random affine projection transformation to the
#' resulting matrix.  This allows RAPT objects serve as a basis
#' for a large number of kinds of expansions.  If p are the
#' number of features of X, and q are number of expanded features,
#' the applications fall into two broad categories:
#'
#' \itemize{
#'   \item{p > q using the Johnson-Lindenstrauss theorem:
#'     \itemize{
#'       \item{Compressed sensing.}
#'       \item{Manifold learning.}
#'       \item{Dimension reduction.}
#'       \item{Graph embedding.}
#'       \item{...}
#'     }
#'   }
#'   \item{p < q using Bochner's theorem:
#'     \itemize{
#'       \item{Approximate kernel projection.}
#'       \item{Fast approximate SVD.}
#'       \item{Estimation of dependence.}
#'       \item{...}
#'     }
#'   }
#' }
#'
#' @param   X       Input data matrix
#' @param   raptObj raptObj generated by raptMake
#'
#' @details Computes
#'
#'          \deqn{X W + b}
#'
#'          where
#'
#'          W = raptObj$W
#'
#'          b = raptObj$b
#'
#' @return  A matrix of randomly (but repeatable) features.
#' @export
#' @seealso Details of how the rapt object is built
#'  are in \code{\link{raptMake}}.
#' @references
#' \url{https://en.wikipedia.org/wiki/Johnson\%E2\%80\%93Lindenstrauss_lemma},
#' \url{https://en.wikipedia.org/wiki/Bochner\%27s_theorem}
#' @examples
#' # Toy problem
#' set.seed(1)
#' nObs <- 100 # Number of observations
#' X <- matrix(seq(-4,+4,length.out = nObs),ncol=1)
#' Ytrue <- sin(5*X) + 2*cos(2*X) # True value Ytrue = g(X)
#' Y <- Ytrue + rnorm(nObs) # Noisy measurement Y
#'
#' # Standardize X
#' Xstd <- scale(X)
#' attributes(Xstd) <- attributes(X)
#'
#' # Bochner (random Fourier) projection object
#' nDim <- NCOL(X)
#' h <- 10 # Estimated by goodness of fit Adj R^2.
#'   # Normally this would be fit by cross validation.
#' raptObj <- raptMake(nDim,nDim*200,WdistOpt = list(sd=h),
#'                     bDistOpt=list(min=-pi,max=+pi))
#'
#' # Apply raptObj to Xstd to generate features,
#' # keeping unaltered features Xstd as well.
#' Xrapt <- cbind( Xstd, cos( rapt(Xstd,raptObj) ) )
#'
#' # Standardize results
#' XraptStd <- scale(Xrapt)
#' attributes(XraptStd) <- attributes(Xrapt)
#'
#' # A linear fitting of Y to the features XraptStd
#' lmObj <- lm(Y ~ XraptStd)
#' summary(lmObj)
#'
#' # Plot measurements (Y), predictions (Yhat),
#' # Kernel smoothing with Gaussian kernel and same bandwidth,
#' # true Y without noise.
#' Yhat <- predict(lmObj)
#' plot (X,Y   ,main="Linear Fitting", ylim=c(-6,10))
#' lines(X,Yhat,col="red",lty=1,lwd=2)
#' grid(col="darkgray")
#' kFit <- ksmooth(X,Y,kernel="normal",bandwidth=1/h)
#' lines(kFit$x,kFit$y,lty=1,col="green",lwd=2)
#' lines(X,Ytrue,lty=1,col="blue",lwd=2)
#' legend("topleft",
#'         legend=c("Noisy measurements",
#'                  "Estimated Y from RAPT",
#'                  "Estimated Y from Kernel Smooth",
#'                  "True Y"),
#'         col=1:4,
#'         pch=c( 1,NA,NA,NA),
#'         lty=c(NA, 1, 1, 1),
#'         lwd=2,
#'         bty="n")
#'
#' # Fit sparse model w/LASSO and
#' # lambda criteria = 1 standard deviation.
#' # This avoids overgeneralization errors usually
#' # associated with fitting large numbers of features
#' # to relatively few data points.  It also improves
#' # the end effects, which are of paramount importance
#' # in high dimensional problems (since by the curse
#' # of dimensionality, almost all points are close an edge
#' # in high dimensional problems).
#' lassoObj <- easyLASSO(XraptStd,Y)
#' Yhat <- predict(lassoObj, newx = XraptStd)
#' # Use linear fit of prediction Yhat as goodness of fit.
#' summary(lm(Y ~ Yhat))
#'
#' # Plot results of LASSO fitting
#' # These show LASSO does a better job fitting edges.
#' plot(X,Y,main="LASSO Fitting",ylim=c(-6,10))
#' lines(X,Yhat,col="red",lty=1,lwd=2)
#' grid(col="darkgray")
#' kFit <- ksmooth(X,Y,kernel="normal",bandwidth=1/h)
#' lines(kFit$x,kFit$y,lty=1,col="green",lwd=2)
#' lines(X,Ytrue,lty=1,col="blue",lwd=2)
#' legend("topleft",
#'         legend=c("Noisy measurements",
#'                  "Estimated Y from RAPT",
#'                  "Estimated Y from Kernel Smooth",
#'                  "True Y"),
#'         col=1:4,
#'         pch=c( 1,NA,NA,NA),
#'         lty=c(NA, 1, 1, 1),
#'         lwd=2,
#'         bty="n")
rapt <- function(X,raptObj) {
  return( sweep(X %*% raptObj$W, 2, raptObj$b, `+`) )
}

#' Freeman-Tukey transform
#'
#' This transform takes Poisson (count) information and
#' makes it more Gaussian, then z-scales (standardizes
#' by centering and scaling to var = 1) the results.
#'
#' @param x Data vector
#' @return The transformed vector
#' @export
#' @references
#' Taken from
#' \url{https://en.wikipedia.org/wiki/Anscombe_transform}
#' @examples
#' x <- freemanTukey(sunspot.month)
freemanTukey <- function(x) {
  x <- as.vector(x)
  x <- sqrt(x) + sqrt(x+1)
  x <- as.vector( scale(x) )
  return(x)
}

#' Plot y and yHat on the same scale w/reference line
#'
#' @description Plots y and yHat on the same scale as a
#'              scatterplot with a 1:1 reference line in red.
#'              This is useful for visually comparing actual
#'              data y with estimates yHat, determining
#'              outliers, etc.
#'
#' @param y    Vector or matrix coercible to vector. Typically
#'             will be the quantity to be predicted.
#' @param yHat Vector or matrix coercible to vector, same
#'             length as y.  Typically will be the prediction.
#' @param ...  Optional additional graph parameters.
#' @details  Normally only makes sense with vectors, column matrices,
#'           or row matrices.
#' @return Returns invisibly - only used for graphic side effects.
#' @export
#' @examples
#' set.seed(1)
#' nObs <- 80
#' X <- distMat(nObs,2)
#' A <- cbind(c(1,-1))
#' Y <- X %*% A + rnorm(nObs) # Response data
#' lmObj <- lm(Y ~ X)
#' Yhat <- predict(lmObj) # Estimated response
#' yyHatPlot(Y,Yhat)
yyHatPlot <- function(y,yHat,...) {

  plot(y,yHat,asp=1,...)
  abline(a=0,b=1,col="red")
  grid(col="darkgray")

  invisible(TRUE)
}

#' Informative plots for Y and Yhat
#'
#' This function presents diagnostic plots of estimate Yhat
#' and response Y.
#'
#' @param   Y    R object representing response,
#'               coercible to a vector.
#' @param   Yhat R object representing estimate,
#'               coercible to a vector.
#'               The length of Y and Yhat must be equal.
#' @param   ...  Options for \code{\link[stats]{cor}} function.
#'               The defaults are use = "everything" and
#'               method = "pearson".
#'
#' @details The plots shown are:
#'
#' \itemize{
#'   \item{Y vs Yhat.  Under a perfect noise-free fitting,
#'          this would be a straight line with
#'          the points lined up on the red line, and the
#'          correlation wpuld be 1.0000.}
#'   \item{Y, Yhat and Y-Yhat (residual) time domain plots.
#'          The time steps are in samples.}
#'   \item{These show the ACF for the original Y, the residual,
#'          and |residual|.  The latter helps identify
#'          nonlinearity in the residual.}
#' }
#'
#' @return  Invisibly returns TRUE; this routine
#'          is only used for its graphical side effects
#'          described in Details.
#'
#' @export
#' @seealso \code{\link[stats]{cor}}
#' @examples
#' # The order here looks backwards, but is chosen to
#' # simulate a typical pair - Yhat will normally have
#' # a smaller range than Y.
#' set.seed(2)
#' nObs <- 100 # Number of observations
#' x <- stats::filter(rnorm(nObs),c(-0.99),
#'      method="recursive")
#' x <- x + (x^2) # Nonlinear component
#' myLags <- 0:2
#' X <- eTrim(eLag(x,myLags))
#' Y <- X[,+1,drop=FALSE]
#' X <- X[,-1,drop=FALSE]
#' lmObj <- lm(Y ~ X)
#' Yhat <- predict(lmObj)
#' Ydiagnostics(Y,Yhat)
Ydiagnostics <- function(Y,Yhat,...) {

  # A useful metric
  myMetric <- cor(Y,Yhat,...)

  # Symmetrical limits on y plots
  temp <- max(abs(range(Y,Yhat)))
  ylim <- c(-temp,temp)

  par(mfrow=c(1,1))
  yyHatPlot(Y, Yhat, xlim=ylim, xlab="Y", ylab="Yhat",
            main=sprintf("Correlation = %4.4f",myMetric))

  myResidual <- Y - Yhat

  par(mfrow=c(2,1))
  matplot(cbind(Y,Yhat),type="l",lty=1,ylim=ylim,
          ylab="Y, Yhat")
  grid(col="darkgray")
  plot(myResidual, type="l",ylab="Residual",ylim=ylim)
  grid(col="darkgray")
  par(mfrow=c(1,1))

  lag.max <- length(Y)

  par(mfrow=c(3,1))
  # Save and set margins
  marOld <- par(mar=c(0,5.1,0,0))
  acf(Y              , lag.max = lag.max, ylab="Y",
      ylim=c(-1,+1), main="ACF(Y)")
  acf(myResidual     , lag.max = lag.max, ylab="Residual",
      ylim=c(-1,+1), main="ACF(Residual)")
  acf(abs(myResidual), lag.max = lag.max, ylab="|Residual|",
      ylim=c(-1,+1), main="ACF(|Residual|)")
  # Restore margins to original values.
  par(mar=marOld)
  par(mfrow=c(1,1))

  invisible(TRUE)
}

#' Plots coefficients in an impulse response format
#'
#' @description Given a model xObj for which coef(xObj)
#'     returns a set of coefficients, plot the coefficients.
#'
#'     The plots make it easier to compare which features are large,
#'     which are set to zero, and how features change from run
#'     to run in a graphical manner.
#'
#'     If the fitting process is linear (e.g. lm, glmnet, etc.)
#'     and the original features are appropriately ordered lags,
#'     this will generate an impulse response.
#'
#'     Any coefficients that are \emph{exactly} zero (for instance,
#'     set that way by LASSO) will appear as red X's; non-zero
#'     points will be black O's.
#'
#' @param xObj             Output of a fitting model.
#' @param includeIntercept Should the 1st coefficient be plotted?
#'                         Default is FALSE.
#' @param type             Graphics type.  Default is "h", which
#'                         results in an impulse-like plot.
#' @param main             "main" title; default is the relative
#'                         number of non-zero coefficients,
#'                         a measure of sparsity.
#' @param ...              Optional additional graphical parameters,
#'                         for instance to set ylim to a fixed value.
#' @details If includeIntercept==TRUE, the intercept of the model
#'          will be plotted as index 0.
#'
#'          Changing the type using \code{type="b"}
#'          will result in a parallel coordinate-like plot rather
#'          than an impulse-like plot.  It is sometimes easier to
#'          see the differences in coefficients with type="b"
#'          rather than type="h".
#'
#' @return Invisibly returns TRUE.  Used for its
#'         graphic side effects only.
#' @export
#' @examples
#' set.seed(1)
#' nObs <- 100
#' X <- distMat(nObs,6)
#' A <- cbind(c(1,0,-1,rep(0,3))) # Y will only depend on X[,1] and X[,3]
#' Y <- X %*% A + 0.1*rnorm(nObs)
#' lassoObj <- easyLASSO(X,Y)
#' Yhat <- predict(lassoObj,newx=X)
#' yyHatPlot(Y,Yhat)
#' coef( lassoObj ) # Sparse coefficients
#' coefPlot( lassoObj )
#' coefPlot( lassoObj, includeIntercept=TRUE )
#' coefPlot( lassoObj, type="b" )
coefPlot <- function(xObj, includeIntercept=FALSE,
                     type = "h",
                     main = NULL, ...) {

  # Get coefficients,including intercept as needed.
  if (includeIntercept) {
    xCoef <- coef(xObj)
    myIndex <- seq_along(xCoef) - 1
  } else {
    xCoef <- coef(xObj)[-1]
    myIndex <- seq_along(xCoef)
  }

  # Which xCoef are exactly zero?
  myZeros <- xCoef==0

  # Default main is 0 (all xCoef are exactly zero) to
  # 1 (all xCoef are used)
  if (is.null(main)) main <- mean(!myZeros)

  # Colors and symbols for zeros and non-zeros
  myCol   <- ifelse(myZeros,"red","black")
  myPch   <- ifelse(myZeros,    4,      1)

  # Plot
  plot  (myIndex, xCoef, type=type, main=mean(xCoef!=0), ...)
  points(myIndex, xCoef, pch= myPch, col=myCol)
  grid(col="darkgray")

  invisible(TRUE)
}

# For STFT and wavelets, note that not only
# does padding have to occur up front, but
# depending on the frame rate used, features
# need to be repeated or interpolated.
# e1071::stft is magnitude only (still useful, though).
# The interface is simple, but also isn't very
# flexible - no ability to specify window, for instance.
# seewave::spectro is the stft, istft is the inverse
# (uses OLA).
# complex can be returned using complex=TRUE.
# temp <- spectro(rnorm(16000),f=1,complex=TRUE,
# norm=FALSE,dB=NULL,plot=FALSE)

# tuneR::melfcc - calculates mel-frequency cepstral coefficients.
# seewave::ceps - cepstrum or real cepstrum.

# Linear SFA can be applied to resulting features.
# NB: SFA that only saves centering and projection
# information w/o saving the input data is much more
# compact.  SFA that has a predict.sfa method is
# also a useful approach.  plot.sfa also can be useful.
# Since SFA is an instantaneous projection method,
# the inputs can contain NAs for predict.sfa
# (although they must be removed for the actual
# fitting process).  summary.sfa is also useful.

# purrr may be able to perform interaction expansion
# by converting formulas to functions.





