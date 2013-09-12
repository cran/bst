\name{ex1data}
\alias{ex1data}
\title{ Generating Three-class Data}
\description{
  Randomly generate data for a three-class model. 

  Modified code from http://www.stat.osu.edu/~yklee/software.html.
}

\usage{
ex1data(n.data)
}
\arguments{
  \item{n.data}{ number of data samples.}
}
\details{
 The data is generated based on Example 1 described in Wang (2011).
}
\value{
  A list with n.data by 50 predictor matrix \code{x}, three-class response \code{y} and conditional probabilities \code{p}.
}
\references{

 Zhu Wang (2012), Multi-class HingeBoost: Method and Application to the Classification of Cancer Types Using Gene Expression Data. \emph{Methods of Information in Medicine}, \bold{51}(2), 162--7.

}

\author{ Zhu Wang }
\examples{
\dontrun{
dat <- ex1data(200)
mhingebst(x=dat$x, y=dat$y)
}
}
\keyword{classification}
