\name{bst_control}
\alias{bst_control}
\title{ Control Parameters for Boosting }
\description{
  Specification of the number of boosting iterations, step size
  and other parameters for boosting algorithms.
}

\usage{
bst_control(mstop = 50, nu = 0.1, twinboost = FALSE, 
f.init = NULL, xselect.init = NULL, center = FALSE, trace = FALSE, 
numsample = 50, df=4)
}
\arguments{
  \item{mstop}{ an integer giving the number of boosting iterations. }
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{twinboost}{ a logical value: \code{TRUE} for twin boosting. }
  \item{f.init}{ the estimate from the first round of twin boosting. Only useful when \code{twinboost=TRUE}. }
  \item{xselect.init}{ the variable selected from the first round of twin boosting. Only useful when \code{twinboost=TRUE}. }
  \item{center}{ a logical value: \code{TRUE} to center covariates with mean. }
  \item{trace}{ a logical value for printout of more details of information during
                the fitting process. }
  \item{numsample}{ number of random sample variable selected in the first round of twin boosting. This is potentially useful in the future implementation. }
  \item{df}{ degree of freedom used in smoothing splines. }
}
\details{
  Objects to specify parameters of the boosting algorithms 
  implemented in \code{\link{bst}}, via the \code{ctrl} argument.
}
\value{
  An object of class \code{bst_control}, a list.
}
\seealso{ \code{\link{bst}}}
