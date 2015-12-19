rmbst <- function(x,y, cost=0.5, rfamily="thinge", threshold=c("adaptive", "standard"), ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree"), del=1e-10){
    call <- match.call()
    learner <- match.arg(learner)
    rfamily <- match.arg(rfamily)
    threshold <- match.arg(threshold)
    s <- ctrl$s
    if(!is.null(s)){
            if(s < 0) stop("s must be >= 0\n")
        }
    fk <- ctrl$fk
    if(is.null(s)){
            s <- switch(rfamily,       
                        "thinge"= 1)
    }
    famtype <- switch(rfamily,
                      "thinge"="thingeDC"
                      )
    ctrl$s <- s
    iter <- ctrl$iter
    trace <- ctrl$trace
    if(trace) cat("\ngenerate initial values\n") 
### initiate values are important, best with nonrobust intercept models
### may need to upgrade for other nonrobust methods
    if(is.null(fk)){
	    #RET <- NULL
	    #RET$k <- length(table(y))
	    #RET$yhat <- matrix(0, ncol=RET$k, nrow=length(y))
	    bsttype <- switch(rfamily,
		    	                  "thinge"="hinge2"
		          )
		    RET <- mbst(x, y, cost=cost, family=bsttype, ctrl = bst_control(mstop=1), control.tree=control.tree, learner=learner)
    }
    else {
        RET <- NULL
        RET$yhat <- fk
    }
    los <- loss.mbst(y, f=RET$yhat, fk=fk, s=ctrl$s, k=RET$k, family = rfamily, cost=cost)
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\nrobust boosting ...\n")
        cat("\ninitial loss", mean(los), "\n")
    }
    los <- rep(NA, iter)
	while(d1 > del && k <= iter){
        ctrl$fk <- RET$yhat
        RET <- mbst(x, y, cost=cost, family=famtype, ctrl = ctrl, control.tree=control.tree, learner=learner)
	los[k] <- mean(loss.mbst(y, f=RET$yhat, fk=NULL, s=ctrl$s, k=RET$k, family = rfamily, cost=cost))
	#difference of convex function is linearly majorized, thus tmp1 - los[k] >= 0. cf Wang (2015)
	if(trace){
	tmp <- matrix(NA, nrow=length(y), ncol=RET$k)
	f <- RET$yhat; fk <- ctrl$fk
	for(j in 1:RET$k)
	   tmp[,j] <- (y!=j)*(mapply(function(x) max(x, 0), f[,j]+1) - mapply(function(x) max(x, 0), fk[,j]-s)- (f[,j]-fk[,j])*(fk[,j] >= s))
        tmp1 <- sum(tmp)/length(y)
 	cat("difference of convex function is linearly majorized, returning a non-negative number", tmp1-los[k], "\n")
        }
	d1 <- sum((RET$yhat - ctrl$fk)^2)/sum(ctrl$fk^2)
        if(trace) cat("\niteration", k, ": relative change of fk", d1, ", robust loss value", los[k], "\n") 
        if(k > 1){
        if(los[k] > los[k-1])
        k <- iter
        }
        k <- k + 1
    }
    RET$x <- x
    RET$y <- y
    RET$call <- call
    RET$cost <- cost
    RET$rfamily <- RET$family <- rfamily
    RET
}

