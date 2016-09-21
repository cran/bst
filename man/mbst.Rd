\name{mbst}
\alias{mbst}
\alias{print.mbst}
\alias{predict.mbst}
\alias{fpartial.mbst}
\title{ Boosting for Multi-Classification}
\description{
  Gradient boosting for optimizing multi-class loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
mbst(x, y, cost = NULL, family = c("hinge", "hinge2", "thingeDC", "closs", "clossMM"), 
ctrl = bst_control(), control.tree=list(fixed.depth=TRUE, 
n.term.node=6, maxdepth = 1), learner = c("ls", "sm", "tree"))
\method{print}{mbst}(x, ...)
\method{predict}{mbst}(object, newdata=NULL, newy=NULL, mstop=NULL, 
type=c("response", "class", "loss", "error"), ...)
\method{fpartial}{mbst}(object, mstop=NULL, newdata=NULL)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be 1, 2, ..., k for a k classification problem}
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ \code{family} = "hinge" for hinge loss, \code{family}="hinge2" for hinge loss but the response is not recoded (see details). \code{family="thingeDC"} for DCB loss function, see \code{rmbst}.} 
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{type}{ in \code{predict} a character indicating whether the response, all responses across the boosting iterations, classes, loss or classification errors should be predicted in case of \code{hinge} 
               problems. in \code{plot}, plot of boosting iteration or $L_1$ norm. }
  \item{object}{ class of \code{\link{mbst}}. }
  \item{newdata}{ new data for prediction with the same number of columns as \code{x}. }
  \item{newy}{ new response. }
  \item{mstop}{ boosting iteration for prediction. }
  \item{\dots}{ additional arguments. }
}
\details{
 A linear or nonlinear classifier is fitted using a boosting algorithm for multi-class responses. This function is different from \code{mhingebst} on how to deal with zero-to-sume constraint and loss functions. If \code{family="hinge"}, the loss function is the same as in \code{mhingebst} but the boosting algorithm is different. If \code{family="hinge2"}, the loss function is different from \code{family="hinge"}: the response is not recoded as in Wang (2012). In this case, the loss function is  
\deqn{\sum{I(y_i \neq j)(f_j+1)_+}.} \code{family="thingeDC"} for robust loss function used in the DCB algorithm.
}
\value{
  An object of class \code{mbst} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{plot}} and \code{\link{predict}} methods are available for linear models.
  For nonlinear models, methods \code{\link{print}} and \code{\link{predict}} are available.
\item{x, y, cost, family, learner, control.tree, maxdepth}{These are input variables and parameters}
\item{ctrl}{the input \code{ctrl} with possible updated \code{fk} if \code{family="thingeDC"}}
\item{yhat}{predicted function estimates}
\item{ens}{a list of length \code{mstop}. Each element is a fitted model to the psedo residuals, defined as negative gradient of loss function at the current estimated function}
\item{ml.fit}{the last element of \code{ens}}
\item{ensemble}{a vector of length \code{mstop}. Each element is the variable selected in each boosting step when applicable}
\item{xselect}{selected variables in \code{mstop}}
\item{coef}{estimated coefficients in each iteration. Used internally only}
}
\references{
  Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.

  Zhu Wang (2012), Multi-class HingeBoost: Method and Application to the Classification of Cancer Types Using Gene Expression Data. \emph{Methods of Information in Medicine}, \bold{51}(2), 162--7.

}
\seealso{\code{\link{cv.mbst}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

\author{ Zhu Wang }
\examples{
x <- matrix(rnorm(100*5),ncol=5)
c <- quantile(x[,1], prob=c(0.33, 0.67))
y <- rep(1, 100)
y[x[,1] > c[1] & x[,1] < c[2] ] <- 2
y[x[,1] > c[2]] <- 3
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
