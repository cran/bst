\name{bfunc}
\alias{bfunc}
\title{ Compute upper bound of second derivative of loss}
\description{
  Compute upper bound of second derivative of loss.
}

\usage{
bfunc(family, s)
}
\arguments{
  \item{family}{ a family from "closs", "gloss", "qloss" for classification and "clossR" for robust regression. }
  \item{s}{ a parameter related to robustness. }
}
\details{
A finite upper bound is required in quadratic majorization.
}
\value{
A positive number.
}
\author{ Zhu Wang }
\keyword{classification, regression}
