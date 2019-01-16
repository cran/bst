\name{rmbst}
\alias{rmbst}
\title{ Robust Boosting for Multi-class Robust Loss Functions}
\description{
  MM (majorization/minization) based gradient boosting for optimizing nonconvex robust loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
rmbst(x, y, cost = 0.5, rfamily = c("thinge", "closs"), ctrl=bst_control(),
control.tree=list(maxdepth = 1),learner=c("ls","sm","tree"),del=1e-10)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, 2, ..., k\}. }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{rfamily}{ \code{family} = "thinge" is currently implemented.}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{del}{convergency critera}
}
\details{
An MM algorithm operates by creating a convex surrogate function that majorizes the nonconvex objective function. When the surrogate function is minimized with gradient boosting algorithm, the desired objective function is decreased. The MM algorithm contains difference of convex (DC) for \code{rfamily="thinge"}, and quadratic majorization boosting algorithm (QMBA) for \code{rfamily="closs"}. 
}
\value{
  An object of class \code{bst} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{plot}} and \code{\link{predict}} methods are available for linear models.
  For nonlinear models, methods \code{\link{print}} and \code{\link{predict}} are available.
\item{x, y, cost, rfamily, learner, control.tree, maxdepth}{These are input variables and parameters}
\item{ctrl}{the input \code{ctrl} with possible updated \code{fk} if \code{type="adaptive"}}
\item{yhat}{predicted function estimates}
\item{ens}{a list of length \code{mstop}. Each element is a fitted model to the psedo residuals, defined as negative gradient of loss function at the current estimated function}
\item{ml.fit}{the last element of \code{ens}}
\item{ensemble}{a vector of length \code{mstop}. Each element is the variable selected in each boosting step when applicable}
\item{xselect}{selected variables in \code{mstop}}
\item{coef}{estimated coefficients in \code{mstop}}
}
\references{

  Zhu Wang (2018),
  Quadratic Majorization for Nonconvex Loss with
                Applications to the Boosting Algorithm, \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3), 491-502, \url{http://dx.doi.org/10.1080/10618600.2018.1424635}
  
  Zhu Wang (2018), Robust boosting with truncated loss functions, \emph{Electronic Journal of Statistics}, \bold{12}(1), 599-650, \url{https://doi.org/10.1214/18-EJS1404}

}

\seealso{\code{\link{cv.rmbst}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

\author{ Zhu Wang }
\examples{
x <- matrix(rnorm(100*5),ncol=5)
c <- quantile(x[,1], prob=c(0.33, 0.67))
y <- rep(1, 100)
y[x[,1] > c[1] & x[,1] < c[2] ] <- 2
y[x[,1] > c[2]] <- 3
x <- as.data.frame(x)
x <- as.data.frame(x)
dat.m <- mbst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- mbst(x, y, ctrl = bst_control(twinboost=TRUE, 
f.init=predict(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rmbst(x, y, ctrl = bst_control(mstop=50, s=1, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)
}
\keyword{classification}
