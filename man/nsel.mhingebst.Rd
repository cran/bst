\name{nsel}
\alias{nsel}

\title{ Find Number of Variables In Multi-class Boosting Iterations}
\description{ Find Number of Variables In Multi-class Boosting Iterations}

\usage{
nsel(object, mstop)
}
\arguments{
  \item{object}{an object of \code{\link{mhingebst}}, \code{\link{mbst}}, or \code{\link{rmbst}}}
  \item{mstop}{boosting iteration number}
}
\value{a vector of length \code{mstop} indicating number of variables selected in each boosting iteration}
\author{ Zhu Wang }

