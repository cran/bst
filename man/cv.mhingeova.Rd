\name{cv.mhingeova}
\alias{cv.mhingeova}
\title{ Cross-Validation for one-vs-all HingeBoost with multi-class problem }
\description{
  Cross-validated estimation of the empirical misclassification error for boosting parameter selection.
}
\usage{
cv.mhingeova(x, y, balance=FALSE, K=10, cost = NULL, nu=0.1, 
learner=c("tree", "ls", "sm"), maxdepth=1, m1=200, twinboost = FALSE, 
m2=200, trace=FALSE, plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of multi class responses. \code{y} must be an integer vector from 1 to C for C class problem. }
  \item{balance}{ logical value. If TRUE, The K
parts were roughly balanced, ensuring that the classes were distributed
proportionally among each of the K parts.}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{maxdepth}{ tree depth used in \code{learner=tree}}
  \item{m1}{ number of boosting iteration }
  \item{twinboost}{ logical: twin boosting? }
  \item{m2}{ number of twin boosting iteration }
  \item{trace}{ if TRUE, iteration results printed out }
  \item{plot.it}{ a logical value, to plot the estimated risks if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{\dots}{ additional arguments. }
}

\value{
  object with
  \item{residmat }{ empirical risks in each cross-validation at boosting iterations }
  \item{fraction}{ abscissa values at which CV curve should be computed. }
  \item{cv}{The CV curve at each value of fraction}
  \item{cv.error}{The standard error of the CV curve}
  ...
}
\seealso{ \code{\link{mhingeova}} }
\note{ The functions for balanced cross validation were from R package pmar. } 
