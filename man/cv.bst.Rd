\name{cv.bst}
\alias{cv.bst}
\title{ Cross-Validation for Binary HingeBoost}
\description{
  Cross-validated estimation of the empirical risk
  for boosting parameter selection.
}
\usage{
cv.bst(x, y, K = 10, cost = 0.5, family = c("hinge", "gaussian"), 
learner = c("ls", "sm", "tree"), ctrl = bst_control(), 
type = c("risk", "misc"), plot.it = TRUE, se = TRUE, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for \code{family} = "hinge".}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ \code{family} = "hinge" for hinge loss and \code{family}="gaussian" for squared error loss.
Implementing the negative gradient corresponding
                to the loss function to be minimized. By default, hinge loss
                for +1/-1 binary responses is used.}
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ cross-validation criteria. For \code{family="hinge"}, \code{type="risk"} is hinge risk and \code{type="misc"} is misclassification error. For \code{family="gaussian"}, only empirical risks. }
  \item{plot.it}{ a logical value, to plot the estimated risks if \code{TRUE}. }
  \item{se}{ a logical value, to plot with standard errors. }
  \item{\dots}{ additional arguments. }
}

\value{
  object with
  \item{residmat }{ empirical risks in each cross-validation at boosting iterations }
  \item{mstop}{ boosting iteration steps at which CV curve should be computed. }
  \item{cv}{The CV curve at each value of mstop}
  \item{cv.error}{The standard error of the CV curve}
  \item{family}{ \code{family} = "hinge" for hinge loss and \code{family}="gaussian" for squared error loss.}
  ...
}
\seealso{ \code{\link{bst}} }
\examples{
\dontrun{
x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
x <- as.data.frame(x)
cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls", type="misc")
}
}
