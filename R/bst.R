require(rpart)
### a: cost of false positive, 0 < a < 1
### If the cost of false negative is twice the cost of false positve, then (1-a)/a=2, thus, a=1/3
loss <- function(y, f, cost = NULL, family = "hinge"){
  if(family == "hinge"){
    if(is.null(cost)) cost <- 0.5
    hingeloss(y, f, a = cost)
  }
  else if(family == "gaussian")
    gaussloss(y, f)
}

gaussloss <- function(y, f) 1/2*(y - f)^2

hingeloss <- function(y, f, w=1, a=0.5){
  if(any(!y %in% c(-1, 1)))
    stop("y must be either 1 or -1\n")
  if(a >=1 || a <= 0)
    stop("misclassification error cost must be between 0 and 1\n")
  tmp <- 1-y*f
  res <- rep(0, length(tmp))
  for(i in 1:length(tmp)){
    if(tmp[i] >0)
      res[i] <- tmp[i]
    if(y[i]==1 && sign(f[i])==-1)
      res[i] <- (1-a) * res[i]
    if(y[i]==-1 && sign(f[i])==1)
      res[i] <- a * res[i]
  }
  res
}

ngradient <- function(y, f, cost = NULL, family = "hinge"){
  if(family == "hinge"){
    if(is.null(cost)) cost <- 0.5
    hingengra(y, f, a = cost)
  }
  else if(family == "gaussian")
    gaussngra(y, f)
}
###negative gradient w.r.t f
gaussngra <- function(y, f) y - f

### cost must be the same as in hinge loss
hingengra <- function(y, f, a=0.5){
  tmp <- 1-y*f
  res <- rep(0, length(tmp))
  for(i in 1:length(tmp)){
    if(tmp[i] >0)
      res[i] <- y[i]
    if(length(f)==length(y)){ ### excluding the initial iteration for which f=0, not a vector
      if(y[i]==1 && sign(f[i])==-1)
        res[i] <- (1-a) * res[i]
      if(y[i]==-1 && sign(f[i])==1)
        res[i] <- a * res[i]
    }
  }
  res
}

bst_control <- function(mstop = 50, nu = 0.1, twinboost = FALSE, f.init = NULL, xselect.init = NULL, center = FALSE, trace = FALSE, numsample = 50, df = 4) {

  RET <- list(mstop = mstop, nu = nu, 
              center = center, df = df,
              trace = trace, twinboost = twinboost, f.init = f.init, xselect.init = xselect.init)
  class(RET) <- c("bst_control")
  RET
}

#######################################################################################################################################################
bst <- function(x,y, cost=0.5, family = c("hinge", "gaussian"), ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree")){
  call <- match.call()
  family <- match.arg(family)
  learner <- match.arg(learner)
  if(learner == "tree" && is.null(colnames(x)))
  colnames(x) <- paste("x", 1:ncol(x), sep = "")
  mstop <- ctrl$mstop
  nu <- ctrl$nu
  twinboost <- ctrl$twinboost
  f.init <- ctrl$f.init
  xselect.init <- ctrl$xselect.init
  center <- ctrl$center
  trace <- ctrl$trace
  numsample <- ctrl$numsample
  df <- ctrl$df
  if(twinboost && (is.null(f.init) | is.null(xselect.init)))
    stop("Twin boosting requires initial function estimates and variable selected in the first round\n")
  nsample <- dim(x)[1]
  p <- dim(x)[2]
  if(learner == "tree" && p > 10 && twinboost && control.tree$maxdepth >= 2 && is.null(numsample))
    stop("for large p and degree >=2, random sample is suggested\n")   
  if(center){
    one <- rep(1,nrow(x))
    meanx <- drop(one %*% as.matrix(x))/length(y)
    x <- scale(x, meanx, FALSE) # centers x
  }
  oldx <- x; one <- rep(1,length(y))
  ens <- vector("list",mstop)
  Fboost <- offset <- pred.val <- 0
  if(family == "gaussian"){
    offset <- mean(y)
    y <- y - offset
  }
  Fboost <- rep(Fboost,length(y))
  ystar <- res <- matrix(NA,length(y),ncol(x))
  m <- 1
  coef0 <- sse <- minid <- rep(NA,ncol(x))
  sse <- minid <- rep(NA,ncol(x))
  risk <- xselect <- coef <- rep(NA, mstop)

  if(learner=="tree"){
    maxdepth <- control.tree$maxdepth
    if(maxdepth==1){
      p1 <- ifelse(!twinboost, p,length(xselect.init))
      xselect <- rep(NA,mstop)
      if(twinboost){
        xselect.new <- xselect.init
        inde <- as.matrix(1:p1, ncol=1)
      }
    }
    else if(twinboost && maxdepth==2){
      if(missing(vimp.init)) vimp.init <- rep(1,length(xselect.init))
      if(p > 10){
        inde <- NULL
        for (i in 1:numsample)
          inde <- rbind(inde, sample(xselect.init,maxdepth,prob=vimp.init[vimp.init>0]))
      }  
      else
        inde <- t(combn(xselect.init,2)) #generate interactions
#      P <- dim(inde)[1]
      xselect <- matrix(NA, ncol=2, nrow=mstop)
      xselect.new <- xselect.init[-(length(xselect.init))]
    }
  }
  while (m <= mstop){
    cor.w <- mse.w <- coef.w <- rep(NA,p)
    u <- ngradient(y, Fboost, cost = cost, family = family)
    if(!twinboost) xselect.init <- 1:p
    if(learner=="ls"){
      for (j in xselect.init){
        coef0[j] <- 1/sum(x[,j]^2)*sum(x[,j] * u)
        pred.tr <- x[,j] * coef0[j] 
        ss <- sum(pred.tr^2)
        if(twinboost){
          cor.w[j] <- cov(f.init,pred.tr)/sqrt(sum(pred.tr^2))
          mse.w[j] <- cor.w[j]^2 * (2*sum(u*pred.tr) - ss)
        }
        else mse.w[j] <- 2*sum(u*pred.tr) - ss
      }
      ind <- which.max(mse.w)
      ml.fit <- lm(u~x[, ind]-1)
      coef[m] <- coef0[ind]
      xselect[m] <- ind
    }
    else
      if(learner=="sm"){
        tree.twin <- vector("list",length(xselect.init))
        for(j in xselect.init){
###Twin L2 Boosting with genral weak learner, Buhlmann, page 8, step 4, in Twin boosting, improved feature selection and prediction
          tree.twin[[j]] <- smooth.spline(x=x[,j],y=u,df=df)
          pred.tr <- fitted(tree.twin[[j]])
          ss <- sum(pred.tr^2)
          if(twinboost){   
            cor.w[j] <- cov(f.init,pred.tr)/sqrt(sum(pred.tr^2))
            mse.w[j] <- cor.w[j]^2 * (2*sum(u*pred.tr) - ss)
          }
          else mse.w[j] <- 2*sum(u*pred.tr) - ss
        }   
        ind <- which.max(mse.w)
        ml.fit <- tree.twin[[ind]]
        xselect[m] <- ind
      }
      else
        if(learner=="tree"){
          cntrl <- rpart.control(maxdepth = maxdepth, #minsplit = nsample-1, #minbucket = 1,
                                 maxsurrogate = 0, maxcompete = 0, #usesurrogate=0
                                 cp = 0, xval = 0)
          if(!twinboost){
            ml.fit <- rpart(u~.,data=data.frame(cbind(u,x)),method="anova",control=cntrl)
            labs <- rownames(ml.fit[["splits"]])
            xselect[m] <- which(colnames(x) %in% labs)
          }
          else{
            tree.twin <- vector("list",nrow(inde))  ### check if correct for maxdepth=1
#            for(j in 1:nrow(inde)){
            for(j in 1:nrow(inde)){
              if(maxdepth==1){
                data.tr <- as.data.frame(cbind(u,x[,xselect.new[j]])); 
                colnames(data.tr) <- c("u",colnames(x)[xselect.new[j]])
              }
              else{
                data.tr <- as.data.frame(cbind(u,x[,inde[j,]])); 
                colnames(data.tr) <- c("u",colnames(x)[inde[j,]])
              }
###Twin L2 Boosting with genral weak learner, Buhlmann, page 8, step 4, in Twin boosting, improved feature selection and prediction
### cf weakboostadapt.rpart in twin.R by Buhlmann, received from Horton in 2009
              tree.twin[[j]] <- rpart(u~.,data=data.tr,method="anova",control=cntrl)
              pred.tr <- predict(tree.twin[[j]])
              ss <- sum(pred.tr^2)
              cor.w[j] <- cov(f.init,pred.tr)/sqrt(ss)
              mse.w[j] <- cor.w[j]^2 * (2*sum(u*pred.tr) - ss)
            }
            ml.fit <- tree.twin[[which.max(mse.w)]]
            if(maxdepth==1) xselect[m] <- which.max(mse.w)
            else {
              tmp <- ml.fit$frame$var[ml.fit$frame$var%in%colnames(x)]
              tmp <- unique(tmp)
              if(length(tmp)!=0)
                xselect[m,] <- as.character(tmp)
            }
          }
        }
### update prediction
    if(learner=="sm")
      Fboost <- Fboost + nu * fitted(ml.fit)
    else
      Fboost <- Fboost + nu * predict(ml.fit)
    risk[m] <- mean(loss(y, Fboost, cost = cost, family = family))
    if(trace){
      if(m %% 10==0) cat("\nm=", m, "  risk = ", risk[m])
    } 
    ens[[m]] <- ml.fit
    m <- m + 1
  }

  ensemble <- xselect
#  if(!twinboost) xselect <- sort(unique(xselect))
  xselect <- sort(unique(xselect))
  RET <- list(y=y,x=oldx, cost=cost, family = family, learner=learner, yhat=Fboost,offset=offset, ens=ens, control.tree=control.tree, risk=risk, ctrl = list(center=center, mstop=mstop,nu=nu, df=df), xselect=xselect, coef = coef, ensemble=ensemble)
  RET$call <- call
  class(RET) <- "bst"
  return(RET)
}

predict.bst <- function(object, newdata=NULL, newy=NULL, mstop=NULL, type=c("response", "all.res", "class", "loss", "error"), ...){
  if(is.null(mstop))
    mstop <- object$ctrl$mstop
  else if(mstop > object$ctrl$mstop)
      stop("mstop must be equal or smaller than the one used for estimation ", object$ctrl$mstop)
#  if((type=="loss" || type=="error") && (is.null(newdata) || is.null(newy)))
#    stop("For estimation of loss or error, both newdata and newy are needed\n")
  if (!is.null(newdata)) {
    if (is.null(colnames(newdata)))
      stop("missing column names for ", sQuote("newdata"))
  }
  type <- match.arg(type)
  one <- rep(1,nrow(object$x))
  x <- object$x
  y <- object$y
  if(is.null(newdata) && is.null(newy))
  ynow <- y
  else ynow <- newy
  if(!missing(newdata)){
    if(object$ctrl$center){
      meanx <- drop(one %*% as.matrix(x))/nrow(x)
      newdata <- scale(newdata, meanx, FALSE) # centers x
    }
  }
  ens <- object$ens
  lp <- object$offset
  nu <- object$ctrl$nu
  cost <- object$cost
  family <- object$family
  if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
  risk <- rep(NA, mstop)
  if(missing(newdata)) res <- matrix(NA, ncol=mstop, nrow=dim(x)[1])
  else res <- matrix(NA, ncol=mstop, nrow=dim(newdata)[1])
  for(m in 1:mstop){
    if(missing(newdata)){
      if(object$learner=="tree") 
        lp <- lp + nu*predict(ens[[m]])
      else lp <- lp + nu*fitted(ens[[m]])
    }
    else
      if(object$learner=="tree") 
        lp <- lp + nu*predict(ens[[m]], newdata = newdata)
      else if(object$learner=="sm")
        lp <- lp + nu * predict(object$ens[[m]], newdata[, object$ensemble[m]])$y
      else if(object$learner=="ls")
        lp <- lp + nu * object$coef[m] * newdata[, object$ensemble[m]]
    if(type=="all.res")
     res[,m] <- lp
     else if(type=="loss"){
#     risk[m] <- mean(loss(newy, lp, cost = cost, family = family))
     risk[m] <- mean(loss(ynow, lp, cost = cost, family = family))
     }
     else if(type == "error"){
      tmp <- sign(lp)
#      risk[m] <- (mean(newy != tmp))
      risk[m] <- (mean(ynow != tmp))
    }
  }

  if(type == "all.res")
   return(res)
  else 
  if(type == "class")
  lp <- sign(lp)
  else if(type %in% c("loss", "error")) lp <- risk
  return(drop(lp))
}

"cv.bst" <-
  function(x, y, K = 10, cost = 0.5, family = c("hinge", "gaussian"), learner = c("tree","ls", "sm"), ctrl = bst_control(), type = c("risk", "misc"), plot.it = TRUE, se = TRUE, ...)
{
  call <- match.call()
  family <- match.arg(family)
  learner <- match.arg(learner)
  type <- match.arg(type)
  family <- match.arg(family)
  learner <- match.arg(learner)
  mstop <- ctrl$mstop
  nu <- ctrl$nu
  df <- ctrl$df
  twinboost <- ctrl$twinboost
  trace <- ctrl$trace
  ctrl.cv <- ctrl
  if(family == "gaussian" && type =="misc") stop("Only risk option is implemented for gaussion family\n")
  all.folds <- cv.folds(length(y), K)
  fraction <- seq(from = 1, to = mstop, by=5)
  m1 <- 1:length(fraction)
  residmat <- matrix(0, length(fraction), K)
  for(i in seq(K)) {
    if(trace)
      cat("\n CV Fold", i, "\n\n")
    omit <- all.folds[[i]]
    if(ctrl$twinboost)
    ctrl.cv$f.init <- ctrl$f.init[ - omit]
    fit <- bst(x[ - omit,,drop=FALSE  ], y[ - omit], cost = cost, family = family, learner = learner, ctrl = ctrl.cv, ...)
    fit <- sapply(fraction, function(m) predict(fit, newdata = x[omit,  ,drop=FALSE], mstop = m))
    if(length(omit)==1)fit<-matrix(fit,nrow=1)
### hinge risk or misclassification error
    if(family == "hinge"){
      if(type == "risk"){
        tmp <- sapply(m1, function(m) loss(y[omit], fit[,m], cost = cost, family = family))
        residmat[, i] <- apply(tmp, 2, mean)
      }
      else{
        tmp <- sapply(m1, function(m) {
               tmp1 <- 0
              for (s in 1:length(omit)){
               if(y[omit[s]] == 1 && fit[s,m] <= 0)
               tmp1 <- tmp1 + (1-cost)
               else if(y[omit[s]] == -1 && fit[s,m] > 0)
               tmp1 <- tmp1 + cost
                }
               tmp1 <- tmp1/length(omit)
               }
               )
        residmat[, i] <- tmp
      }
    }
    else{
      tmp <- sapply(m1, function(m) loss(y[omit], fit[,m], cost = cost, family = family))
      residmat[, i] <- apply(tmp, 2, mean)
    }
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(residmat = residmat, fraction = fraction, cv = cv, cv.error = cv.error)
  if(plot.it) plotCVbst(object,se=se)
  invisible(object)
}

print.bst <- function(x, ...) {

  cat("\n")
  cat("\t Models Fitted with Gradient Boosting\n")
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  show(x$family)
  cat("\n")
  if(!is.null(x$ctrl$twinboost))
    cat("Twin boosting", "\n")
  cat("Base learner: ", x$learner, "\n")
  cat("Number of boosting iterations: mstop =", x$ctrl$mstop, "\n")
  cat("Step size: ", x$ctrl$nu, "\n")
  cat("Offset: ", x$offset, "\n")
  cat("\n")
  if(x$learner=="ls"){
    cat("Coefficients: \n")
    cf <- coef(x)
    print(cf)
    cat("\n")
  }
  if(x$learner=="sm")
  cat("Degree of freedom used is: ", x$ctrl$df, "\n")
  invisible(x)
}


### methods: coefficients
coef.bst <- function(object, ...) {
  if(object$learner != "ls")
    stop("Coefficients only implemented for linear least squares\n")
  ret <- numeric(NCOL(object$x))
  xselect <- object$ensemble
  for (j in unique(xselect)){
    tmp <- which(xselect == j)
    ret[j] <- sum(object$coef[tmp])
  }
  names(ret) <- colnames(object$x)
  RET <- ret * object$ctrl$nu
  attr(RET, "offset") <- object$offset
  RET
}

coefpath.bst <- function(object, ...) {
  if(object$learner != "ls")
    stop("Coefficients only implemented for linear least squares\n")

  vars <- colnames(object$x)
  xselect <- object$ensemble
  svars <- vars[tabulate(xselect, nbins = length(vars)) > 0]
  ret <- matrix(0, nrow = object$ctrl$mstop, ncol = length(svars))
  colnames(ret) <- svars
  for (j in unique(xselect)) {
    indx <- which(xselect == j)
    ret[indx, svars[svars == vars[j]]] <-
      object$coef[indx]
  }
  RET <- ret * object$ctrl$nu
  apply(RET, 2, cumsum)
}

plot.bst <- function(x, type=c("step","norm"), ...) {
  type <- match.arg(type)
  if(x$learner != "ls")
    stop("Coefficients path only implemented for linear least squares\n")
  cp <- coefpath.bst(x)
  cf <- cp[nrow(cp),]
  p <- dim(cp)[2]
  if(p > 1 && is.null(col)) col=2:(p+1)
  else col = 1
  x <- apply(cp, 1, function(x) sum(abs(x)))
  if(type == "step")
    matplot(cp, type = "l", col = col, xlab = "Number of boosting iterations",
            ylab = "Coefficients", ...)
  else matplot(x, cp, type = "l", col = col, xlab = "L_1 norm",
               ylab = "Coefficients", ...)
  axis(4, at = cp[nrow(cp),],labels = colnames(cp))
}

