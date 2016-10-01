#' expandFunctions: a feature matrix builder
#'
#' @description
#' A variety of  functions for conversion of
#' vectors and matrices to other matrices to use as
#' features.  This allows one to quickly build feature structures
#' and apply various machine learning methods to those features
#' for exploration and pedantic purposes.
#'
#' @details
#' The \strong{expandFunctions} package contains functions
#' that can be used to expand feature vectors and matrices into
#' larger feature matrices.  These functions include lag
#' embedding, special function univariate exansion, quadratic
#' expansion, and random vector projection.
#'
#' The general steps for feature generation for time domain data
#' (which subsumes multivariate data via lags) are:
#'
#' \itemize{
#'   \item{Preprocess data - remove mean, transform, etc., to a
#'   useful vector or matrix.}
#'   \item{If not a matrix, functionally expand vector into a matrix.
#'   This is
#'   typically done by lag embedding, but may also include STFT, wavelet
#'   transforms, etc.}
#'   \item{Functionally expand matrices generated.}
#'   \item{Combine resulting matrices into a single feature matrix.}
#'   \item{Dimensional reduction, feature selection, and/or
#'   feature extraction to reduce the number of features.}
#'   \item{Use machine learning method(s) on the resulting feature
#'   matrix.}
#' }
#'
#' Most of the steps above are well supported in \R on CRAN, but the
#' expansion steps tend to be scattered in a variety of packages,
#' poorly represented, or custom built by the user.  The
#' \strong{expandFunction} package is intended
#' to collect many of these functions together in one place.
#'
#' Preprocessing almost always should include centering and scaling the
#' data.  However, it may also include a variety of transformations,
#' such as Freeman-Tukey, in order to make the vector fit
#' more closely to a given model (say, a linear model with Gaussian
#' noise).
#'
#' If the input isn't a time domain vector, but is instead already
#' in tabular form (for instance, Boston Housing Data),
#' the embedding step can be skipped.
#'
#' Dimension reduction is outside the scope of this package, but
#' is normally performed to reduce the variables that need handling,
#' reducing the memory used and speeding up the analysis.
#'
#' The package functions are "magrittr-friendly", that is,
#' built so that they can be directly pipelined since X, the data,
#' is the first argument.
#'
#' Most functions are prefixed with "e" to help distinguish them
#' from being confused with similarly named functions.
#'
#' @examples
#' \dontrun{
#' # Sunspot counts can be somewhat Gaussianized by the
#' # Freeman-Tukey transform.
#' x <- freemanTukey(sunspot.month)
#' par(mfrow=c(1,1)) # just in case multiplots were left over.
#' plot(x,type="l")
#'
#' # Embed x using eLag
#' # Since the base period of sunspots is 11*12 months,
#' # pick the lags to be fractions of this.
#' myLags <- seq(from=0,to=200,by=1)
#' X <- eTrim(eLag(x,myLags))
#' Y <- X[,+1,drop=FALSE]
#' X <- X[,-1,drop=FALSE]
#' # matplot(X,type="l",lty=1)
#'
#' # OLS fitting on the lag feature set
#' lmObj <- lm(y ~ .,data = data.frame(x=X,y=Y))
#' coefPlot(lmObj,type="b")
#' summary(lmObj)
#' Yhat <- predict(lmObj, newdata = data.frame(x=X))
#' Ydiagnostics(Y,Yhat)
#'
#' # LASSO fitting on the lag feature set
#' lassoObj <- easyLASSO(X,Y,criterion="lambda.min")
#' coefPlot(lassoObj,type="b")
#' Yhat <- predict(lassoObj,newx = X)
#' Ydiagnostics(Y,Yhat)
#'
#' # Reduce the lag feature set using non-zero
#' # LASSO coefficients
#' useCoef <- coef(lassoObj)[-1]!=0
#' myLags <- seq(from=0,to=200,by=1)[c(TRUE,useCoef)]
#' X <- eTrim(eLag(x,myLags))
#' Y <- X[,+1,drop=FALSE]
#' X <- X[,-1,drop=FALSE]
#'
#' # OLS fitting on the reduced lag feature set
#' lmObj <- lm(y ~ .,data = data.frame(x=X,y=Y))
#' summary(lmObj)
#' coefPlot(lmObj)
#' Yhat <- predict(lmObj, newdata = data.frame(x=X))
#' Ydiagnostics(Y,Yhat)
#'
#' # 1st nonlinear feature set
#' # Apply a few Chebyshev functions to the columns of the
#' # lag matrix. Note these exclude the constant values,
#' # but include the linear.
#' chebyFUN <- polywrapper(basePoly=orthopolynom::chebyshev.t.polynomials,
#'                         kMax=5)
#' Z <- eMatrixOuter(X,1:5,chebyFUN)
#'
#' # OLS fitting on the 1st nonlinear feature set
#' lmObj <- lm(y ~ .,data = data.frame(z=Z,y=Y))
#' summary(lmObj)
#' Yhat <- predict(lmObj, newdata = data.frame(z=Z))
#' Ydiagnostics(Y,Yhat)
#'
#' # LASSO fitting on the 1st nonlinear feature set
#' lassoObj <- easyLASSO(Z,Y)
#' coefPlot(lassoObj)
#' Yhat <- predict(lassoObj,newx = Z)
#' Ydiagnostics(Y,Yhat)
#'
#' # 2nd nonlinear feature set
#' # Use eQuad as an alternative expansion of the lags
#' Z <- eQuad(X)
#'
#' # OLS fitting on the 2nd nonlinear feature set
#' lmObj <- lm(y ~ .,data = data.frame(z=Z,y=Y))
#' summary(lmObj)
#' Yhat <- predict(lmObj, newdata = data.frame(z=Z))
#' Ydiagnostics(Y,Yhat)
#'
#' # LASSO fitting on the 2nd nonlinear feature set
#' lassoObj <- easyLASSO(Z,Y)
#' coefPlot(lassoObj)
#' Yhat <- predict(lassoObj,newx = Z)
#' Ydiagnostics(Y,Yhat)
#' }
#'
#' @name    expandFunctions
#' @docType package
#' @author  Scott Miller <sam3CRAN@gmail.com>
#' @import  utils stats graphics orthopolynom plyr
NULL
