mhinge <- function(xtr, ytr, xte=NULL, yte=NULL, cost = NULL, nu=0.1, learner="ls", m1=200, K=10, cv1=FALSE, twin = FALSE, m2=200, cv2=FALSE){
  call <- match.call()
  if(m1 < 1 || m2 < 1) stop("Too small m1 or m2\n")
  n1 <- length(ytr)
  ncla <- length(unique(ytr))
  if(min(ytr)!=1 || max(ytr)!=ncla)
    stop("Response variable must be 1, 2, ...\n")
  if(!is.null(xte) && dim(xtr)[2] != dim(xte)[2])
    stop("The training data and test data have different dimensions\n")  
  if(is.null(cost))
    cost.mul <- rep(0.5, ncla)
  else {
    cost.mul <- cost
    if(any(cost <= 0) || any(cost >= 1))
    stop("The costs must in (0, 1)\n")
  }
  z <- matrix(-1, ncol=ncla, nrow=n1)
  fit.tr <- fit.te <- xsel <- NULL
  for(k in 1:ncla){
    z[which(ytr==k), k] <- 1
    cost1 <- cost.mul[k]
    m1up <- m1
    if(cv1){
      eval(parse(text = paste("tr",k, " <- cv.bst(x = xtr, y = z[,k], K=K, cost=cost.mul[k], learner=learner, ctrl = bst_control(mstop=m1, nu=nu))", sep=""))) 
      m1up <- eval(parse(text = paste("which.min(tr",k, "$cv)", sep="")))
    } 
    eval(parse(text = paste("tr",k, " <- bst(x = xtr, y = z[,k],  cost=cost.mul[k], learner=learner, ctrl = bst_control(mstop=m1up, nu=nu))", sep=""))) 
    if(twin){
     m2up <- m2
     if(cv2){
      eval(parse(text = paste("cvtr <- cv.bst(x = xtr, y = z[,k], K=K, cost=cost.mul[k], learner=learner, ctrl = bst_control(twinboost=TRUE, f.init=predict(tr", k, "), xselect.init = tr", k, "$xselect, mstop=m2, nu=nu))", sep=""))) 
      m2up <- eval(parse(text = paste("which.min(cvtr$cv)", sep="")))
     }
      eval(parse(text = paste("tr",k, " <- bst(x = xtr, y = z[,k],  cost=cost.mul[k], learner=learner, ctrl = bst_control(twinboost=TRUE, f.init=predict(tr", k, "), xselect.init = tr", k, "$xselect, mstop=m2up, nu=nu))", sep=""))) 
    } 
    fit.tr <- cbind(fit.tr, eval(parse(text=paste("predict(tr",k,")", sep=""))))           
    xsel <- c(xsel, eval(parse(text=paste("tr",k,"$xselect", sep=""))))           
    if(!is.null(xte)){
      fit.te <- cbind(fit.te, eval(parse(text=paste("predict(tr",k,", newdata=xte)", sep=""))))           
    }
  }
  err.tr <- sum(ytr != apply(fit.tr, 1, which.max))/length(ytr)
  if(!is.null(yte) && !is.null(xte))
    err.te <- sum(yte != apply(fit.te, 1, which.max))/length(yte)
  RET=list(call=call, learner=learner, nu=nu, twin=twin, m1=m1, m2=m2, cv1=cv1, cv2=cv2, err.tr = err.tr, err.te = err.te, xsel=round(table(xsel)/ncla,2))
  class(RET) = "mhinge"
  return(RET)
       }

print.mhinge <- function(x, ...) {

  cat("\n")
  cat("\t Multi-class HingeBoost Fitted with One-against-All\n")
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  cat("Base learner: ", x$learner, "\n")
  cat("Step size: ", x$nu, "\n")
  if(!cv1)
  cat("Number of boosting iterations: mstop =", x$m1, "\n")
  if(x$twin)
    cat("Twin boosting", "\n")
  if(!cv2)
  cat("Number of twin boosting iterations: mstop =", x$m2, "\n")
  if(x$learner=="sm")
  cat("Degree of freedom used is: ", x$ctrl$df, "\n")
  invisible(x)
}


