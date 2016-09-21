\name{cv.mbst}
\alias{cv.mbst}
\title{ Cross-Validation for Multi-class Boosting}
\description{
  Cross-validated estimation of the empirical multi-class loss
  for boosting parameter selection.
}
\usage{
cv.mbst(x, y, balance=FALSE, K = 10, cost = NULL, 
family = c("hinge","hinge2","thingeDC", "closs", "clossMM"), 
learner = c("tree", "ls", "sm"), ctrl = bst_control(), 
type = c("loss","error"), plot.it = TRUE, se = TRUE, n.cores=2, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be integers from 1 to C for C class problem. }
  \item{balance}{ logical value. If TRUE, The K
parts were roughly balanced, ensuring that the classes were distributed
proportionally among each of the K parts.}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ \code{family} = "hinge" for hinge loss. "hinge2" is a different hinge loss}
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ for \code{family="hinge"}, \code{type="loss"} is hinge risk. For \code{family="thingeDC"}, \code{type="loss"}}
  \item{plot.it}{ a logical value, to plot the estimated risks if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{n.cores}{The number of CPU cores to use. The cross-validation loop
  will attempt to send different CV folds off to different cores.}
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
\seealso{ \code{\link{mbst}} }
