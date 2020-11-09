\name{cv.mada}
\alias{cv.mada}
\title{ Cross-Validation for one-vs-all AdaBoost with multi-class problem }
\description{
  Cross-validated estimation of the empirical misclassification error for boosting parameter selection.
}
\usage{
cv.mada(x, y, balance=FALSE, K=10, nu=0.1, mstop=200, interaction.depth=1, 
trace=FALSE, plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{ a data matrix containing the variables in the model.}
  \item{y}{ vector of multi class responses. \code{y} must be an integer vector from 1 to C for C class problem.}
  \item{balance}{ logical value. If TRUE, The K parts were roughly balanced, ensuring that the classes were distributed proportionally among each of the K parts.}
  \item{K}{ K-fold cross-validation }
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{mstop}{ number of boosting iteration. }
  \item{interaction.depth}{ used in gbm to specify the depth of trees. }
  \item{trace}{ if TRUE, iteration results printed out. }
  \item{plot.it}{ a logical value, to plot the cross-validation error if \code{TRUE}. }
  \item{se}{ a logical value, to plot with 1 standard deviation curves. }
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
\seealso{ \code{\link{mada}} }
