\name{evalerr}
\alias{evalerr}
\title{ Compute prediction errors }
\description{
  Compute prediction errors for classification and regression problems.
}

\usage{
evalerr(family, y, yhat)
}
\arguments{
  \item{family}{ a family used in \code{bst}. Classification or regression family. }
  \item{y}{ response variable. For classification problems, y must be 1/-1. }
  \item{yhat}{ predicted values. }
}
\details{
For classification, returns misclassification error. For regression, returns mean squared error. 
}
\value{
For classification, returns misclassification error. For regression, returns mean squared error. 
}
\author{ Zhu Wang }
\keyword{internal}
