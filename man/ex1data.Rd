\name{ex1data}
\alias{ex1data}
\title{ Generating Three-class Data with 50 Predictors}
\description{
  Randomly generate data for a three-class model. 
}

\usage{
ex1data(n.data, p=50)
}
\arguments{
  \item{n.data}{ number of data samples.}
  \item{p}{ number of predictors. }
}
\details{
 The data is generated based on Example 1 described in Wang (2012).
}
\value{
  A list with n.data by p predictor matrix \code{x}, three-class response \code{y} and conditional probabilities.
}
\references{

 Zhu Wang (2012), Multi-class HingeBoost: Method and Application to the Classification of Cancer Types Using Gene Expression Data. \emph{Methods of Information in Medicine}, \bold{51}(2), 162--7.

}

\author{ Zhu Wang }
\examples{
\dontrun{
dat <- ex1data(100, p=5)
mhingebst(x=dat$x, y=dat$y)
}
}
\keyword{classification}
