\name{cv.rbst}
\alias{cv.rbst}
\title{ Cross-Validation for Nonconvex Loss Boosting}
\description{
  Cross-validated estimation of the empirical risk/error, can be used for tuning parameter selection.
}
\usage{
cv.rbst(x, y, K = 10, cost = 0.5, rfamily = c("tgaussian", "thuber", "thinge", 
"tbinom", "binomd", "texpo", "tpoisson", "clossR", "closs", "gloss", "qloss"), 
learner = c("ls", "sm", "tree"), ctrl = bst_control(), type = c("loss", "error"), 
plot.it = TRUE, main = NULL, se = TRUE, n.cores=2,...)
}
\arguments{
  \item{x}{ a data frame containing the variables in the model.}
  \item{y}{ vector of responses. \code{y} must be in \{1, -1\} for binary classification}
  \item{K}{ K-fold cross-validation }
  \item{cost}{ price to pay for false positive, 0 < \code{cost} < 1; price of false negative is 1-\code{cost}.}
  \item{rfamily}{ nonconvex loss function types. }
  \item{learner}{ a character specifying the component-wise base learner to be used:
      \code{ls} linear models,
      \code{sm} smoothing splines,
      \code{tree} regression trees.
}
  \item{ctrl}{ an object of class \code{\link{bst_control}}.}
  \item{type}{ cross-validation criteria. For \code{type="loss"}, loss function values and \code{type="error"} is misclassification error. }
  \item{plot.it}{ a logical value, to plot the estimated loss or error with cross validation if \code{TRUE}. }
  \item{main}{ title of plot}
  \item{se}{ a logical value, to plot with standard errors. }
 \item{n.cores}{The number of CPU cores to use. The cross-validation loop
  will attempt to send different CV folds off to different cores.}
  \item{\dots}{ additional arguments. }
}

\value{
  object with
  \item{residmat }{ empirical risks in each cross-validation at boosting iterations }
  \item{mstop}{ boosting iteration steps at which CV curve should be computed. }
  \item{cv}{The CV curve at each value of mstop}
  \item{cv.error}{The standard error of the CV curve}
  \item{rfamily}{ nonconvex loss function types. }
  ...
}
\seealso{ \code{\link{rbst}} }
\author{ Zhu Wang }
\examples{
\dontrun{
x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
x <- as.data.frame(x)
cv.rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls", type="lose")
cv.rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls", type="error")
dat.m <- rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls")
dat.m1 <- cv.rbst(x, y, ctrl = bst_control(twinboost=TRUE, coefir=coef(dat.m), 
xselect.init = dat.m$xselect, mstop=50), family = "thinge", learner="ls")
}
}
