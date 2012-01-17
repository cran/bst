mada <- function(xtr, ytr, xte=NULL, yte=NULL, mstop=50, nu=0.1, interaction.depth=1){
  call <- match.call()
  ncla <- length(unique(ytr))
  z <- matrix(0, ncol=ncla, nrow=length(ytr))
  if(!is.null(xte) && !is.null(yte))
    z.te <- matrix(0, ncol=ncla, nrow=length(yte))
  fit.tr <- fit.te <- fpar <- vector("list", ncla)
  err.tr <- err.te <- risk.te <- xsel <- NULL
  msel <- vector("list", ncla)
  for(k in 1:ncla){
    z[which(ytr==k), k] <- 1
    if(!is.null(xte) && !is.null(yte))
      z.te[which(yte==k), k] <- 1 
    dat <- data.frame(z=z[,k], xtr=xtr)
    colnames(dat) <- c("z", paste("x", 1:ncol(xtr), sep = ""))
    colnames(xte) <- paste("x", 1:ncol(xte), sep="")
    eval(parse(text = paste("tr",k, " <- gbm(z ~ ., data=dat, distribution='adaboost', shrinkage=0.1, n.trees=mstop, interaction.depth=interaction.depth, verbose=FALSE)", sep="")))
    msel[[k]] <- matrix(NA, nrow=dim(xtr)[2], ncol=mstop)
    for(m in 1:mstop){
    res <- gbm(z ~ ., data=dat, distribution='adaboost', shrinkage=0.1, n.trees=m, interaction.depth=interaction.depth, verbose=FALSE)
    xselect <- summary(res, order=FALSE,plotit=FALSE)[,2]
    msel[[k]][, m] <- ifelse(xselect > 0, 1, 0)
      #variable selected if xselect=1, o.w. 0.
}
      eval(parse(text = paste("fit.tr[[",k, "]] <- predict(tr", k, ",newdata=xtr, n.trees=1:mstop)", sep="")))
    if(mstop > 1)
    err.tr <- unlist(lapply(1:mstop, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.tr[[i]][,j]); 1/length(ytr) * sum(ytr != apply(tmp, 1, which.max)) }))
    else err.tr <- unlist(lapply(1:mstop, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.tr[[i]][j]); 1/length(ytr) * sum(ytr != apply(tmp, 1, which.max)) }))
    if(!is.null(yte) && !is.null(xte)){
      eval(parse(text = paste("fit.te[[",k, "]] <- predict(tr", k, ",newdata=xte, n.trees=1:mstop)", sep="")))
  }
    if(mstop > 1)
    err.te <- unlist(lapply(1:mstop, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.te[[i]][,j]); 1/length(yte) * sum(yte != apply(tmp, 1, which.max)) }))
    else err.te <- unlist(lapply(1:mstop, function(j){tmp <- NULL; for(i in 1:ncla) tmp <- cbind(tmp, fit.te[[i]][j]); 1/length(yte) * sum(yte != apply(tmp, 1, which.max)) }))
}
  xsel <- rep(NA, mstop)
  tmp <- NULL
  for(m in 1:mstop){
   for(k in 1:ncla)
    tmp <- c(tmp, which(msel[[k]][, m] == 1))
   xsel[m] <- length(unique(tmp))
  }
 RET <- list(xselect=xsel, err.tr = err.tr, err.te = err.te)
}

cv.mada <- function(x, y, balance=FALSE, K=10, nu=0.1, mstop=200, interaction.depth=1, trace=FALSE, plot.it = TRUE, se = TRUE, ...)
{
  call <- match.call()
  if(balance)
  all.folds <- balanced.folds(y, K)
  else
  all.folds <- cv.folds(length(y), K)
  fraction <- 1:mstop
#  fraction <- seq(from = 1, to = m1, by=5)
  residmat <- matrix(0, length(fraction), K)
  for(i in seq(K)) {
    if(trace)
      cat("\n CV Fold", i, "\n\n")
    omit <- all.folds[[i]]
    fit <- mada(xtr = x[ - omit,,drop=FALSE  ], ytr = y[ - omit], xte = x[ omit,,drop=FALSE ], yte=y[ omit ], mstop = mstop, nu=nu, interaction.depth=interaction.depth)
    residmat[,i] <- fit$err.te
  if(trace && i==K)
  cat("End of cross-validation\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(residmat = residmat, fraction = fraction, cv = cv, cv.error = cv.error)
  if(plot.it) plotCVbst(object,se=se)
  invisible(object)
}


