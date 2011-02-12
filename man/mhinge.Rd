\name{mhinge}
\alias{mhinge}
\alias{print.mhinge}
\title{ Multi-class HingeBoost}
\description{
  Multi-class HingeBoost with one-against-all which optimizes the hinge loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
mhinge(xtr, ytr, xte=NULL, yte=NULL, cost = NULL, nu=0.1, 
learner="ls", m1=200, K=10, cv1=FALSE, twin = FALSE, m2=200, cv2=FALSE)
\method{print}{mhinge}(x, ...)
}
\arguments{
  \item{xtr}{ training data containing the predictor variables.}
  \item{ytr}{ vector of training data responses. \code{ytr} must be in \{1,2,...,k\}.}
  \item{xte}{ test data containing the predictor variables.}
  \item{yte}{ vector of test data responses. \code{yte} must be in \{1,2,...,k\}.}
  \item{cost}{ default is NULL for equal cost; otherwise a numeric vector indicating price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{m1}{ number of boosting iteration }
  \item{K}{ K-fold cross-validation }
  \item{cv1}{ logical. if TRUE, then tuning parameter is selected by cross-validation in the first round of boosting.} 
  \item{twin}{ logical: twin boosting? }
  \item{m2}{ number of twin boosting iteration }
  \item{cv2}{ logical. If TRUE, then tuning parameter is selected by cross-validation in twin boosting.} 
  \item{x}{ class of \code{\link{mhinge}}. }
  \item{\dots}{ additional arguments. }
}
\details{
 A linear or nonlinear multi-class HingeBoost classifier is fitted using a boosting algorithm based on one-against component-wise
  base learners for +1/-1 responses, with possible cost-sensitive hinge loss function.
}
\value{
  An object of class \code{mhinge} with \code{\link{print}} method being available.
}
\references{

  Peter Buhlmann and Torsten Hothorn (2007),
  Boosting algorithms: regularization, prediction and model fitting.
  \emph{Statistical Science}, \bold{22}(4), 477--505.
  
  Peter Buhlmann and Torsten Hothorn (2010),
  Twin Boosting: improved feature selection and prediction.
  \emph{Statistics and Computing}, \bold{20}(2), 119--138.
 
  Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.
}
\seealso{\code{\link{bst}} for HingeBoost binary classification. Furthermore see \code{\link{cv.bst}} for stopping iteration selection by cross-validation, and \code{\link{bst_control}} for control parameters.}

\author{ Zhu Wang }
\examples{
\dontrun{
dat1 <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/
thyroid-disease/ann-train.data")
dat2 <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/
thyroid-disease/ann-test.data")
res <- mhinge(xtr=dat1[,-22], ytr=dat1[,22], xte=dat2[,-22], yte=dat2[,22], 
cost=c(2/3, 0.5, 0.5), nu=0.5, learner="ls", m1=100, K=5, cv1=FALSE, 
twin=TRUE, m2= 200, cv2=FALSE)
res <- mhinge(xtr=dat1[,-22], ytr=dat1[,22], xte=dat2[,-22], yte=dat2[,22], 
cost=c(2/3, 0.5, 0.5), nu=0.5, learner="ls", m1=100, K=5, cv1=FALSE, 
twin=TRUE, m2= 200, cv2=TRUE)
}
}
\keyword{classification}
