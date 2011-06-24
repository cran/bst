\name{cv.msvm}
\alias{cv.msvm}
\title{ Cross-Validation }
\description{
  Cross-validated estimation of the empirical multi-class hinge loss
  for boosting parameter selection.
}
\usage{
cv.msvm(x, y, K = 10, cost = NULL, family = "hinge", 
learner = c("tree", "ls", "sm"), ctrl = bst_control(), 
type = c("risk","misc"), plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for \code{family} = "hinge".}
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
  \item{fraction }{ number of boosting iterations }
  \item{cv }{ empirical risks with cross-validation at boosting iteration \code{fraction}. }
  \item{cv.err }{ empirical risks standard errors with cross-validation at boosting iteration \code{fraction}. }
  ...
}
\seealso{ \code{\link{bst}} }
\examples{
x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
x <- as.data.frame(x)
cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls", type="misc")
}
