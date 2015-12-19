\name{rbstpath}
\alias{rbstpath}
\title{ Robust Boosting Path for Truncated Loss Functions}
\description{
  Gradient boosting path for optimizing robust loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
rbstpath(x, y, rmstop=seq(40, 400, by=20), ctrl=bst_control(), del=1e-16, ...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\}. }
  \item{rmstop}{ vector of boosting iterations}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{del}{convergency critera}
  \item{...}{arguments passed to \code{rbst}}
}
\details{
 This function invokes \code{rbst} with \code{mstop} being each element of vector \code{rmstop}. It can provide different paths. Thus \code{rmstop} serves as another hyper-parameter. However, the most important hyper-parameter is the loss truncation point.  
}
\value{
  A length \code{rmstop} vector of lists with each element being an object of class \code{rbst}.
}
\seealso{\code{\link{rbst}}}

\author{ Zhu Wang }
\examples{
x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
y[1:10] <- -y[1:10]
x <- as.data.frame(x)
dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
coefir=coef(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rbst(x, y, ctrl = bst_control(mstop=50, s=0, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)
rmstop <- seq(10, 40, by=10)
dat.m3 <- rbstpath(x, y, rmstop, rfamily = "thinge", learner = "ls")
}
\keyword{classification}
