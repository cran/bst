\name{rbst}
\alias{rbst}
\title{ Robust Boosting for Robust Loss Functions}
\description{
  MM (majorization/minimization) algorithm based gradient boosting for optimizing nonconvex robust loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
rbst(x, y, cost = 0.5, rfamily = c("tgaussian", "thuber","thinge", "tbinom", "binomd", 
"texpo", "tpoisson", "clossR", "closs", "gloss", "qloss"), ctrl=bst_control(), 
control.tree=list(maxdepth = 1), learner=c("ls","sm","tree"),del=1e-10)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for classification. }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{rfamily}{ robust loss function, see details. }
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{del}{convergency criteria}
}
\details{
An MM algorithm operates by creating a convex surrogate function that majorizes the nonconvex objective function. When the surrogate function is minimized with gradient boosting algorithm, the desired objective function is decreased. The MM algorithm contains difference of convex (DC) algorithm for \code{rfamily=c("tgaussian", "thuber","thinge", "tbinom", "binomd", "texpo", "tpoisson")} and quadratic majorization boosting algorithm (QMBA) for \code{rfamily=c("clossR", "closs", "gloss", "qloss")}. 
  
\code{rfamily} = "tgaussian" for truncated square error loss, "thuber" for truncated Huber loss, "thinge" for truncated hinge loss, "tbinom" for truncated logistic loss, "binomd" for logistic difference loss, "texpo" for truncated exponential loss, "tpoisson" for truncated Poisson loss, "clossR" for C-loss in regression, "closs" for C-loss in classification, "gloss" for G-loss, "qloss" for Q-loss.

\code{s} must be a numeric value to be specified in \code{bst_control}. For \code{rfamily="thinge", "tbinom", "texpo"}  \code{s < 0}. For \code{rfamily="binomd", "tpoisson", "closs", "qloss", "clossR"} , \code{s > 0} and for \code{rfamily="gloss"}, \code{s > 1}. Some suggested \code{s} values: "thinge"= -1, "tbinom"= -log(3), "binomd"= log(4), "texpo"= log(0.5), "closs"=1, "gloss"=1.5, "qloss"=2, "clossR"=1.
}
\value{
  An object of class \code{bst} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{plot}} and \code{\link{predict}} methods are available for linear models.
  For nonlinear models, methods \code{\link{print}} and \code{\link{predict}} are available.
\item{x, y, cost, rfamily, learner, control.tree, maxdepth}{These are input variables and parameters}
\item{ctrl}{the input \code{ctrl} with possible updated \code{fk} if \code{family="tgaussian"}, \code{"thingeDC"}, \code{"tbinomDC"}, \code{"binomdDC"} or \code{"tpoisson"}.}
\item{yhat}{predicted function estimates}
\item{ens}{a list of length \code{mstop}. Each element is a fitted model to the pseudo residuals, defined as negative gradient of loss function at the current estimated function}
\item{ml.fit}{the last element of \code{ens}}
\item{ensemble}{a vector of length \code{mstop}. Each element is the variable selected in each boosting step when applicable}
\item{xselect}{selected variables in \code{mstop}}
\item{coef}{estimated coefficients in \code{mstop}}
}

\references{

  Zhu Wang (2018),
  Quadratic Majorization for Nonconvex Loss with
                Applications to the Boosting Algorithm, \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3), 491-502, \url{https://doi.org/10.1080/10618600.2018.1424635}

  Zhu Wang (2018), Robust boosting with truncated loss functions, \emph{Electronic Journal of Statistics}, \bold{12}(1), 599-650, \url{https://doi.org/10.1214/18-EJS1404}

}

\seealso{\code{\link{cv.rbst}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

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
}
\keyword{classification}
