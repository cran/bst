\name{msvm}
\alias{msvm}
\alias{print.msvm}
\alias{predict.msvm}
\alias{fpartial.msvm}
\title{ Boosting for Multi-class Classification}
\description{
  Gradient boosting for optimizing multi-class hinge loss functions with componentwise
  linear, smoothing splines, tree models as base learners.
}

\usage{
msvm(x, y, cost = NULL, family = c("hinge"), ctrl = bst_control(), 
control.tree = list(maxdepth = 1), learner = c("ls", "sm", "tree"))
\method{print}{msvm}(x, ...)
\method{predict}{msvm}(object, newdata=NULL, newy=NULL, mstop=NULL, type=c("response", "class", "loss", "error"), ...)
\method{fpartial}{msvm}(object, mstop=NULL, newdata=NULL)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for \code{family} = "hinge".}
  \item{cost}{ equal costs for now and unequal costs will be implemented in the future.}
  \item{family}{ \code{family} = "hinge" for multi-class hinge loss. }
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{control.tree}{ control parameters of rpart. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models, 
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{type}{ in \code{predict} a character indicating whether the response, classes, loss or classification errors should be predicted in case of \code{hinge}} 
               problems. 
  \item{object}{ class of \code{\link{msvm}}. }
  \item{newdata}{ new data for prediction with the same number of columns as \code{x}. }
  \item{newy}{ new response. }
  \item{mstop}{ boosting iteration for prediction. }
  \item{\dots}{ additional arguments. }
}
\details{
 A linear or nonlinear classifier is fitted using a boosting algorithm based on component-wise
  base learners for multi-class responses. 
}
\value{
  An object of class \code{msvm} with \code{\link{print}} and \code{\link{predict}} methods being available for fitted models.
}
\references{
 
  Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.

  Zhu Wang (2011), Multi-class HingeBoost: Method and Application to the Classification of Cancer Types Using Gene Expression Data. Manuscript.

}
\seealso{\code{\link{cv.msvm}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

\author{ Zhu Wang }
\examples{
dat <- ex1data(200)
res <- msvm(x=dat$x, y=dat$y)
predict(res)
}
\keyword{classification}
