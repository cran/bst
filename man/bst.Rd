\name{bst}
\alias{bst}
\alias{print.bst}
\alias{predict.bst}
\alias{plot.bst}
\alias{coef.bst}
\alias{fpartial.bst}
\title{ Boosting for Classification and Regression}
\description{
  Gradient boosting for optimizing hinge or squared error loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
bst(x, y, cost = 0.5, family = c("hinge", "gaussian"), ctrl = bst_control(), 
control.tree = list(maxdepth = 1), learner = c("ls", "sm", "tree"))
\method{print}{bst}(x, ...)
\method{predict}{bst}(object, newdata=NULL, newy=NULL, mstop=NULL, 
type=c("response", "all.res", "class", "loss", "error"), ...)
\method{plot}{bst}(x, type = c("step", "norm"),...)
\method{coef}{bst}(object, ...)
\method{fpartial}{bst}(object, mstop=NULL, newdata=NULL)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for \code{family} = "hinge".}
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ \code{family} = "hinge" for hinge loss and \code{family}="gaussian" for squared error loss. 
Implementing the negative gradient corresponding
                to the loss function to be minimized. By default, hinge loss
                for +1/-1 binary responses is used.}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{type}{ in \code{predict} a character indicating whether the response, all responses across the boosting iterations, classes, loss or classification errors should be predicted in case of \code{hinge} 
               problems. in \code{plot}, plot of boosting iteration or $L_1$ norm. }
  \item{object}{ class of \code{\link{bst}}. }
  \item{newdata}{ new data for prediction with the same number of columns as \code{x}. }
  \item{newy}{ new response. }
  \item{mstop}{ boosting iteration for prediction. }
  \item{\dots}{ additional arguments. }
}
\details{
 A linear or nonlinear classifier is fitted using a boosting algorithm for +1/-1 responses. 
}
\value{
  An object of class \code{bst} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{plot}} and \code{\link{predict}} methods are available for linear models.
  For nonlinear models, methods \code{\link{print}} and \code{\link{predict}} are available.
}
\references{

  Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.
}
\seealso{\code{\link{cv.bst}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

\author{ Zhu Wang }
\examples{
x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
x <- as.data.frame(x)
dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
f.init=predict(dat.m), xselect.init = dat.m$xselect, mstop=50))
}
\keyword{classification}
