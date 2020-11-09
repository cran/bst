\name{bst}
\alias{bst}
\alias{print.bst}
\alias{predict.bst}
\alias{plot.bst}
\alias{coef.bst}
\alias{fpartial.bst}
\title{ Boosting for Classification and Regression}
\description{
  Gradient boosting for optimizing loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
bst(x, y, cost = 0.5, family = c("gaussian", "hinge", "hinge2", "binom", "expo", 
"poisson", "tgaussianDC", "thingeDC", "tbinomDC", "binomdDC", "texpoDC", "tpoissonDC",
 "huber", "thuberDC", "clossR", "clossRMM", "closs", "gloss", "qloss", "clossMM",
"glossMM", "qlossMM", "lar"), ctrl = bst_control(), control.tree = list(maxdepth = 1), 
learner = c("ls", "sm", "tree"))
\method{print}{bst}(x, ...)
\method{predict}{bst}(object, newdata=NULL, newy=NULL, mstop=NULL, 
type=c("response", "all.res", "class", "loss", "error"), ...)
\method{plot}{bst}(x, type = c("step", "norm"),...)
\method{coef}{bst}(object, which=object$ctrl$mstop, ...)
\method{fpartial}{bst}(object, mstop=NULL, newdata=NULL)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for \code{family} = "hinge".}
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{family}{ A variety of loss functions.
\code{family} = "hinge" for hinge loss and \code{family}="gaussian" for squared error loss. 
Implementing the negative gradient corresponding
                to the loss function to be minimized. For hinge loss, +1/-1 binary responses is used.}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ type of prediction or plot, see \code{\link{predict}}, \code{\link{plot}}}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{object}{ class of \code{\link{bst}}. }
  \item{newdata}{ new data for prediction with the same number of columns as \code{x}. }
  \item{newy}{ new response. }
  \item{mstop}{ boosting iteration for prediction. }
  \item{which}{ at which boosting \code{mstop} to extract coefficients. }
  \item{\dots}{ additional arguments. }
}
\details{
 Boosting algorithms for classification and regression problems. In a classification problem, suppose \eqn{f}{f} is a classifier for a response \eqn{y}{y}. A cost-sensitive or weighted loss function is \deqn{L(y,f,cost)=l(y,f,cost)\max(0, (1-yf))}{L(y,f,cost)=l(y,f,cost)max(0, (1-yf)).} For \code{family="hinge"}, 
\deqn{l(y,f,cost)= 
1-cost, if \, y= +1;
\quad cost, if \, y= -1}{l(y,f,cost)= 1-cost, if y = +1; = cost, if y = -1.} For \code{family="hinge2"}, 
l(y,f,cost)= 1, if y = +1 and f > 0 ; = 1-cost, if y = +1 and f < 0; = cost, if y = -1 and f > 0; = 1, if y = -1 and f < 0. 

For twin boosting if \code{twinboost=TRUE}, there are two types of adaptive boosting if \code{learner="ls"}: for \code{twintype=1}, weights are based on coefficients in the first round of boosting; for \code{twintype=2}, weights are based on predictions in the first round of boosting. See Buehlmann and Hothorn (2010).
}
\value{
  An object of class \code{bst} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{plot}} and \code{\link{predict}} methods are available for linear models.
  For nonlinear models, methods \code{\link{print}} and \code{\link{predict}} are available.
\item{x, y, cost, family, learner, control.tree, maxdepth}{These are input variables and parameters}
\item{ctrl}{the input \code{ctrl} with possible updated \code{fk} if \code{family="thingeDC", "tbinomDC", "binomdDC"}}
\item{yhat}{predicted function estimates}
\item{ens}{a list of length \code{mstop}. Each element is a fitted model to the pseudo residuals, defined as negative gradient of loss function at the current estimated function}
\item{ml.fit}{the last element of \code{ens}}
\item{ensemble}{a vector of length \code{mstop}. Each element is the variable selected in each boosting step when applicable}
\item{xselect}{selected variables in \code{mstop}}
\item{coef}{estimated coefficients in each iteration. Used internally only}
}
\references{

  Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.

  Peter Buehlmann and Torsten Hothorn (2010),
  Twin Boosting: improved feature selection and prediction,
  \emph{Statistics and Computing}, \bold{20}, {119-138}.
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
predict(dat.m)
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
coefir=coef(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rbst(x, y, ctrl = bst_control(mstop=50, s=0, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)
}
\keyword{classification}
