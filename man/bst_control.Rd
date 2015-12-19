\name{bst_control}
\alias{bst_control}
\title{ Control Parameters for Boosting }
\description{
  Specification of the number of boosting iterations, step size
  and other parameters for boosting algorithms.
}

\usage{
bst_control(mstop = 50, nu = 0.1, twinboost = FALSE, twintype=1, 
threshold=c("standard", "adaptive"), f.init = NULL, coefir = NULL,
xselect.init = NULL, center = FALSE, trace = FALSE, numsample = 50, df = 4, 
s = NULL, sh = NULL, q = NULL, qh = NULL, fk = NULL, iter = 10, intercept = FALSE)
}
\arguments{
  \item{mstop}{ an integer giving the number of boosting iterations. }
  \item{nu}{ a small number (between 0 and 1) defining the step size or shrinkage parameter. }
  \item{twinboost}{ a logical value: \code{TRUE} for twin boosting. }
  \item{twintype}{ for \code{twinboost=TRUE} only. For \code{learner="ls"}, if \code{twintype=1}, twin boosting with weights from magnitude of coefficients in the first round of boosting. If \code{twintype=2}, weights are correlations between predicted values in the first round of boosting and current predicted values. For learners not componentwise least squares, \code{twintype=2}. }
  \item{threshold}{ if \code{threshold="adaptive"}, the estimated function \code{ctrl$fk} is updated in every boosting step. Otherwise, no update for \code{ctrl$fk} in boosting steps. Only used if in robust loss functions with the difference convex loss. }
  \item{f.init}{ the estimate from the first round of twin boosting. Only useful when \code{twinboost=TRUE} and \code{learner="sm" or "tree"}. }
  \item{coefir}{ the estimated coefficients from the first round of twin boosting. Only useful when \code{twinboost=TRUE} and \code{learner="ls"}. }
  \item{xselect.init}{ the variable selected from the first round of twin boosting. Only useful when \code{twinboost=TRUE}. }
  \item{center}{ a logical value: \code{TRUE} to center covariates with mean. }
  \item{trace}{ a logical value for printout of more details of information during
                the fitting process. }
  \item{numsample}{ number of random sample variable selected in the first round of twin boosting. This is potentially useful in the future implementation. }
  \item{df}{ degree of freedom used in smoothing splines. }
  \item{s,q}{ truncation parameter \code{s} or frequency \code{q} of outliers for robust regression and classification. If \code{s} is missing but \code{q} is available, \code{s} may be computed as the \code{1-q} quantile of robust loss values using conventional software.}
  \item{sh, qh}{ threshold value or frequency \code{qh} of outliers for Huber regressionn \code{family="huber"} or \code{family="rhuberDC"}.
For \code{family="huber"}, if \code{sh} is not provided, \code{sh} is then updated adaptively with the median of \code{y-yhat} where \code{yhat} is the estimated \code{y} in the last boosting iteration. For \code{family="rhuberDC"}, if \code{sh} is missing but \code{qh} is available, \code{sh} may be computed as the \code{1-qh} quantile of robust loss values using conventional software.}
  \item{fk}{ used for robust classification. A function estimate used in difference of convex algorithm}
  \item{iter}{number of iteration in difference of convex algorithm}
  \item{intercept}{ logical value, if TRUE, estimation of intercept with linear predictor model}
}

\details{
 Objects to specify parameters of the boosting algorithms
  implemented in \code{\link{bst}}, via the \code{ctrl} argument. The default value of \code{s} is -1 if \code{family="thinge"}, -log(3) if \code{family="tbinom"}, and 4 if \code{family="binomd"}
}
\value{
  An object of class \code{bst_control}, a list. Note \code{fk} may be updated for robust boosting.
}
\seealso{ \code{\link{bst}}}
