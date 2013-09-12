\name{mhingebst}
\alias{mhingebst}
\alias{print.mhingebst}
\alias{predict.mhingebst}
\alias{fpartial.mhingebst}
\title{ Boosting for Multi-class Classification}
\description{
  Gradient boosting for optimizing multi-class hinge loss functions with componentwise linear least squares, smoothing splines and trees as base learners.
}

\usage{
mhingebst(x, y, cost = NULL, family = c("hinge"), ctrl = bst_control(), 
control.tree = list(fixed.depth=TRUE, n.term.node=6, maxdepth = 1), 
learner = c("ls", "sm", "tree"))
\method{print}{mhingebst}(x, ...)
\method{predict}{mhingebst}(object, newdata=NULL, newy=NULL, mstop=NULL, 
type=c("response", "class", "loss", "error"), ...)
\method{fpartial}{mhingebst}(object, mstop=NULL, newdata=NULL)
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
  \item{object}{ class of \code{\link{mhingebst}}. }
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
  An object of class \code{mhingebst} with \code{\link{print}} and \code{\link{predict}} methods being available for fitted models.
}
\references{
 
 Zhu Wang (2011),
  HingeBoost: ROC-Based Boost for Classification and Variable Selection.
  \emph{The International Journal of Biostatistics}, \bold{7}(1), Article 13.

 Zhu Wang (2012), Multi-class HingeBoost: Method and Application to the Classification of Cancer Types Using Gene Expression Data. \emph{Methods of Information in Medicine}, \bold{51}(2), 162--7.

}
\seealso{\code{\link{cv.mhingebst}} for cross-validated stopping iteration. Furthermore see
  \code{\link{bst_control}}}

\author{ Zhu Wang }
\examples{
\dontrun{
dat <- ex1data(100)
res <- mhingebst(x=dat$x, y=dat$y)
}
}
\keyword{classification}
