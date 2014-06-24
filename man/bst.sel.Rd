\name{bst.sel}
\alias{bst.sel}
\title{Function to select number of predictors} 
\description{Function to determine the first q predictors in the boosting path, or perform (10-fold) cross-validation and determine the optimal set of parameters}
\usage{
bst.sel(x, y, q, type=c("firstq", "cv"), ...)
}
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Continuous response vector for linear regression}
  \item{q}{Maximum number of predictors that should be selected.}
  \item{type}{if \code{type="firstq"}, return the first \code{q} predictors in the boosting path. if \code{type="cv"}, perform (10-fold) cross-validation and determine the optimal set of parameters}
  \item{...}{Further arguments to be passed to \code{\link{cv.bst}.}}
}
\details{Function to determine the first q predictors in the boosting path, or perform (10-fold) cross-validation and determine the optimal set of parameters. This may be used for p-value calculation. See below.}
\value{Vector of selected predictors.}
\author{Zhu Wang}

\examples{
\dontrun{
x <- matrix(rnorm(100*100), nrow = 100, ncol = 100)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
sel <- bst.sel(x, y, q=10)
library("hdi")
fit.multi <- hdi(x, y, method = "multi.split",
model.selector =bst.sel,
args.model.selector=list(type="firstq", q=10))
fit.multi
fit.multi$pval[1:10] ## the first 10 p-values
fit.multi <- hdi(x, y, method = "multi.split",
model.selector =bst.sel,
args.model.selector=list(type="cv"))
fit.multi
fit.multi$pval[1:10] ## the first 10 p-values
}
}
\keyword{models}
\keyword{regression}
