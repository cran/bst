\name{cv.mhingebst}
\alias{cv.mhingebst}
\title{ Cross-Validation for Multi-class Hinge Boosting}
\description{
  Cross-validated estimation of the empirical multi-class hinge loss
  for boosting parameter selection.
}
\usage{
cv.mhingebst(x, y, balance=FALSE, K = 10, cost = NULL, family = "hinge", 
learner = c("tree", "ls", "sm"), ctrl = bst_control(), 
type = c("risk","misc"), plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be integers from 1 to C for C class problem. }
  \item{balance}{ logical value. If TRUE, The K
parts were roughly balanced, ensuring that the classes were distributed
proportionally among each of the K parts.}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ \code{family} = "hinge" for hinge loss.}
Implementing the negative gradient corresponding
                to the loss function to be minimized. 
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ for \code{family="hinge"}, \code{type="risk"} is hinge risk.}
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
\seealso{ \code{\link{mhingebst}} }
