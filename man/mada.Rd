\name{mada}
\alias{mada}
\title{ Multi-class AdaBoost }
\description{ One-vs-all multi-class AdaBoost }
\usage{
mada(xtr, ytr, xte=NULL, yte=NULL, mstop=50, nu=0.1, interaction.depth=1)
}
\arguments{
  \item{xtr}{ training data matrix containing the predictor variables in the model.}
  \item{ytr}{ training vector of responses. \code{ytr} must be integers from 1 to C, for C class problem. }
  \item{xte}{ test data matrix containing the predictor variables in the model.}
  \item{yte}{ test vector of responses. \code{yte} must be integers from 1 to C, for C class problem. }
  \item{mstop}{ number of boosting iteration. }
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{interaction.depth}{ used in gbm to specify the depth of trees. }
}
\details{
For a C-class problem (C > 2), each class is separately compared against all other classes with AdaBoost, and C functions are estimated to represent confidence for each class. The classification rule is to assign the class with the largest estimate.
}
\value{
  A list contains variable selected \code{xselect} and training and testing error \code{err.tr, err.te}. }  
 
\seealso{\code{\link{cv.mada}} for cross-validated stopping iteration. }

\author{ Zhu Wang }
\keyword{classification}
