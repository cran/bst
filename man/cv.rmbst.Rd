\name{cv.rmbst}
\alias{cv.rmbst}
\title{ Cross-Validation for Truncated Multi-class Loss Boosting}
\description{
  Cross-validated estimation of the empirical multi-class loss
  for boosting parameter selection.
}
\usage{
cv.rmbst(x, y, balance=FALSE, K = 10, cost = NULL, 
rfamily = "thinge", learner = c("tree", "ls", "sm"), ctrl = bst_control(), 
type = c("loss","error"), plot.it = TRUE, main = NULL, se = TRUE, n.cores=2, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be integers from 1 to C for C class problem. }
  \item{balance}{ logical value. If TRUE, The K
parts were roughly balanced, ensuring that the classes were distributed
proportionally among each of the K parts.}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{rfamily}{ \code{rfamily} = "thinge" for truncated multi-class hinge loss.}
Implementing the negative gradient corresponding
                to the loss function to be minimized. 
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ loss value or misclassification error. }
  \item{plot.it}{ a logical value, to plot the estimated loss or error with cross validation if \code{TRUE}. }
  \item{main}{ title of plot}
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
