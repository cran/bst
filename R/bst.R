loss <- function(y, f, cost = 0.5, family = c("gaussian", "hinge", "hinge2", "binom", "expo", "poisson", "thingeDC", "tgaussianDC", "tpoissonDC", "tbinomDC", "binomdDC", "texpoDC", "thinge", "tgaussian", "tbinom", "binomd", "texpo", "tpoisson", "huber", "thuber", "thuberDC"), s=-1, sh=NULL, fk=NULL){
    family <- match.arg(family)
    if(family!="gaussian") ly <- ifelse(y==1, 1-cost, cost)
    if(family == "hinge"){
        hingeloss(y, f, a = cost)
    }
    else if(family == "hinge2"){
        hingeloss2(y, f, a = cost)
    }
    else if(family == "thingeDC"){
        hingeloss(y, f, a = cost) + ly*y*f*(s-y*fk > 0)
    }
    else if(family == "binom"){
        ly*log(1+exp(-y*f)) 
    }
    else if(family == "expo"){
        ly*exp(-y*f)
    }
    else if(family == "tbinomDC"){
        res <- log(1+exp(-y*f))+y*f*exp(-y*fk)/(1+exp(-y*fk))*(y*fk < s)  
        ly*res
    }
    else if(family == "texpoDC"){
        res <- exp(-y*f)+y*f*exp(-y*fk)*(y*fk < s)  
        ly*res
    }
    else if(family == "texpo"){
        ly*pmin(exp(-y*f), exp(-s))
    }
    else if(family == "binomdDC"){
        res <- log(1+exp(-y*f))+y*f*exp(-y*fk-s)/(1+exp(-y*fk-s))
        ly*res
    }
    else if(family=="thinge"){
        tmp <- 1-y*f
        tmp <- tmp*(tmp > 0)
        tmp1 <- s - y*f
        tmp1 <- tmp1*(tmp1 > 0)
        ly*(tmp - tmp1)
    }
    else if(family == "tbinom"){
        ly*pmin(log(1+exp(-y*f)), log(1+exp(-s)))
    }
    else if(family == "binomd"){
        res <- log(1+exp(-y*f))- log(1+exp(-y*f-s))
        ly*res
    }
    else if(family == "tbinomDC"){
        res <- log(1+exp(-y*f))+y*f*exp(-y*fk)/(1+exp(-y*fk))*(y*fk < s)  
        ly*res
    }
    else if(family == "binomdDC"){
        res <- log(1+exp(-y*f))+y*f*exp(-y*fk-s)/(1+exp(-y*fk-s))
        ly*res
    }
    else if(family == "gaussian"){
        gaussloss(y, f)
    }
    else if(family == "tgaussianDC"){
                                        #        if (is.null(s) || s < 0){
### fk is the previous fitted f
                                        #            if(is.null(fk)) fk <- 0
                                        #            s <- quantile(gaussloss(y, fk), 0.8)   ### test
                                        #        }
        gaussloss(y, f)+(y-fk)*(gaussloss(y,fk) > s)*f
    }
    else if(family == "tgaussian"){
                                        #        if (is.null(s) || s < 0){
### fk is the previous fitted f
                                        #            if(is.null(fk)) fk <- 0
                                        #            s <- quantile(gaussloss(y, fk), 0.8)   ### test
                                        #        }
        pmin(gaussloss(y, f), s)
    }
    else if(family == "poisson")
        -y*f+exp(f)
                                        #-dpois(y, exp(f), log = TRUE)
    else if(family == "tpoissonDC")
        -y*f+exp(f) - (-y+exp(fk))*(-y*fk+exp(fk)>s)*f
    else if(family == "tpoisson")
        pmin(-y*f+exp(f), s)
    else if(family == "huber"){
        if (is.null(sh) || sh < 0){
            if(is.null(fk)) fk <- 0
            sh <- median(abs(y - fk)) ### fk was updated in each boosting iteration with adaptive s
        }
        ifelse((a <- abs(y - f)) <= sh, a^2/2, sh*(a - sh/2))
    }
    else if(family == "thuberDC"){
        s2 <- sh
        f1 <- ifelse((a <- abs(y - f)) <= s2, a^2/2, s2*(a - s2/2))
        f2 <- ifelse(abs(y - fk) <= s2, y - fk, s2 * sign(y - fk))
        f3 <- ifelse((a <- abs(y - fk)) <= s2, a^2/2, s2*(a - s2/2))
        f1 + f2*(f3 > s)*f
    }
    else if(family == "thuber"){
### Todo: what if no s provided?
                                        #if (is.null(s) || s < 0){
                                        #    if(is.null(fk)) fk <- 0
                                        #    s <- median(abs(y - fk)) ### fk was updated in each boosting iteration with adaptive s
                                        #}
        b <- ifelse((a <- abs(y - f)) <= sh, a^2/2, sh*(a - sh/2))
        pmin(b, s)
    }

}

gaussloss <- function(y, f) 1/2*(y - f)^2

hingeloss <- function(y, f, w=1, a=0.5){
    if(any(!y %in% c(-1, 1)))
        stop("y must be either 1 or -1\n")
    if(a >=1 || a <= 0)
        stop("misclassification error cost must be between 0 and 1\n")
    tmp <- 1-y*f
    tmp <- tmp*(tmp > 0)
    ly <- ifelse(y==1, 1-a, a)
    ly*tmp 
}

hingeloss2 <- function(y, f, w=1, a=0.5){
    if(any(!y %in% c(-1, 1)))
        stop("y must be either 1 or 1\n")
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


ngradient <- function(y, f, cost = 0.5, family = c("gaussian", "hinge", "hinge2", "binom", "expo", "poisson", "thingeDC", "tgaussianDC", "tpoissonDC", "tbinomDC", "texpoDC", "binomdDC", "huber", "thuberDC"), s=-1, sh=NULL, fk=NULL){
    family <- match.arg(family)
    if(family!="gaussian") ly <- ifelse(y==1, 1-cost, cost)
    if(family == "hinge"){
        hingengra(y, f, a = cost)
    }
    else if(family == "hinge2")
        hingengra2(y, f, a = cost)
    else if(family == "thingeDC"){
        hingengra(y, f, a = cost) - ly*y*(s-y*fk > 0)
    }
    else if(family == "binom")
        ly*y*exp(-y*f)/(1+exp(-y*f))
    else if(family == "expo")
        ly*y*exp(-y*f)
    else if(family == "texpoDC")
        ly*y*exp(-y*f) - ly*y*exp(-y*fk)*(y*fk < s)
    else if(family == "tbinomDC")
        ly*y*exp(-y*f)/(1+exp(-y*f)) - ly*y*exp(-y*fk)/(1+exp(-y*fk))*(y*fk < s)
    else if(family == "binomdDC")
        ly*y*exp(-y*f)/(1+exp(-y*f)) - ly*y*exp(-y*fk-s)/(1+exp(-y*fk-s))
    else if(family == "gaussian")
        gaussngra(y, f)
    else if(family == "tgaussianDC"){
                                        #        if (is.null(s) || s < 0){
### fk is the previous fitted f
                                        #            if(is.null(fk)) fk <- 0
                                        #            s <- quantile(gaussloss(y, fk), 0.8)   ### test
                                        #        }
        gaussngra(y, f) - (y-fk)*(gaussloss(y, fk) > s)
    }
    else if(family == "poisson")
        y-exp(f)
    else if(family == "tpoissonDC")
        (y-exp(f))-(y-exp(fk))*(-y*fk+exp(fk) > s)
    else if(family == "huber"){
        if (is.null(sh) || sh < 0){
### fk is the previous fitted f
            if(is.null(fk)) fk <- 0
            sh <- median(abs(y - fk))
        }
        ifelse(abs(y - f) < sh, y - f, sh * sign(y - f))
    }
    else if(family == "thuberDC"){
        f1 <- ifelse(abs(y - f) < sh, y - f, sh * sign(y - f))
        f2 <- ifelse(abs(y - fk) < sh, y - fk, sh * sign(y - fk))
        f1 - f2*(loss(y, fk, family="huber", sh=sh) > s)
    }
}
###negative gradient w.r.t f
gaussngra <- function(y, f) y - f

### cost must be the same as in hinge loss
hingengra <- function(y, f, a=0.5){
    tmp <- ifelse(1-y*f > 0, 1, 0)
    ly <- ifelse(y==1, 1-a, a)
    ly*y*tmp 
}

### cost must be the same as in hinge loss
hingengra2 <- function(y, f, a=0.5){
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

bst_control <- function(mstop = 50, nu = 0.1, twinboost = FALSE, twintype=1, threshold=c("standard", "adaptive"), f.init = NULL, coefir = NULL, xselect.init = NULL, center = FALSE, trace = FALSE, numsample = 50, df = 4, s=NULL, sh=NULL, q=NULL, qh=NULL, fk=NULL, iter=10, intercept=FALSE) {
    threshold <- match.arg(threshold)
    RET <- list(mstop = mstop, nu = nu, 
                center = center, df = df, threshold=threshold,
                trace = trace, twinboost = twinboost, twintype=twintype, f.init = f.init, coefir = coefir, xselect.init = xselect.init, s=s, sh=sh, q=q, qh=qh, fk=fk, iter=iter, intercept=intercept)
    class(RET) <- c("bst_control")
    RET
}

#######################################################################################################################################################

rbstpath <- function(x, y, rmstop=seq(40, 400, by=20), ctrl=bst_control(), del=1e-16, ...){
    fit <- vector("list", length(rmstop))
    for(jk in 1:length(rmstop)){
        if(jk==1) fk <- NULL
        else fk <- fit[[jk-1]]$yhat
        ctrl$mstop <- rmstop[jk]
        ctrl$fk <- fk
        fit[[jk]] <- rbst(x, y, ctrl = ctrl, del=del, ...)
    }
    fit
}

rbst <- function(x,y, cost=0.5, rfamily=c("tgaussian", "thuber", "thinge", "tbinom", "binomd", "texpo", "tpoisson"), ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree"), del=1e-10){
    call <- match.call()
    learner <- match.arg(learner)
    rfamily <- match.arg(rfamily)
    s <- ctrl$s
    sh <- ctrl$sh
    q <- ctrl$q
    qh <- ctrl$qh
    if(is.null(ctrl$s) && rfamily=="tgaussian")
        sa <- TRUE  ### adaptive s
    else sa <- FALSE
                                        #if(is.null(s) && is.null(q)) stop("s or q must be provided\n")
    if(!is.null(q)){
        if(q < 0 || q > 1)
            stop("proportion of outliers must between 0 and 1\n")
        if(rfamily=="thuber"){
            if(is.null(qh)) stop("qh must be provided\n")
            if(q <= qh) stop("q should be larger than qh\n")
        }
    }
    else if(!is.null(s)){
        if(rfamily=="tbinom")
            if(s > 0) stop("s must be non-negative for rf='tbinom'\n")
        if(rfamily=="binomd")
            if(s <= 0) stop("s must be positive for rf='binomd'\n")
        if(rfamily=="thuber"){
            if(is.null(sh)) stop("sh must be provided\n")
            if(s <= sh) stop("s should be larger than sh\n")
        }
    }
    fk <- ctrl$fk
    if(is.null(s)){
        if(rfamily=="tgaussian" || rfamily=="thuber"){
            stop("how to find s is not implemented\n")
        }
        else
            s <- switch(rfamily,       
                        "thinge"= -1,
                        "tbinom"= -log(3),
                        "binomd"= log(4),
                        "texpo"= log(0.5),
                        "tpoisson"= 5*mean(y))
    }
    famtype <- switch(rfamily,
                      "tgaussian"="tgaussianDC",
                      "thuber"="thuberDC",
                      "thinge"="thingeDC",
                      "tbinom"="tbinomDC",
                      "binomd"="binomdDC",
                      "texpo"="texpoDC",
                      "tpoisson"="tpoissonDC",
                      )
    ctrl$s <- s
    ctrl$sh <- sh
    iter <- ctrl$iter
    trace <- ctrl$trace
    if(trace) cat("\ngenerate initial values\n") 
    ly <- ifelse(y==1, 1-cost, cost)
### initiate values are important, best with nonrobust intercept models
### may need to upgrade for other nonrobust methods
    if(is.null(fk)){
        bsttype <- switch(rfamily,
                                        # "tgaussian"="gaussian",
                          "tgaussian"="huber",
                          "thuber"="huber",
                          "thinge"="hinge",
                          "tbinom"="binom",
                          "binomd"="binom",
                          "texpo"="expo",
                          "tpoisson"="poisson",
                          )
        RET <- bst(x, y, cost=cost, family=bsttype, ctrl = bst_control(mstop=1), control.tree=control.tree, learner=learner)
    }
    else {
        RET <- NULL
        RET$yhat <- fk
    }
    los <- loss(y, f=RET$yhat, cost, family = rfamily, s=ctrl$s, sh=ctrl$sh, fk=fk)
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\nrobust boosting ...\n")
        cat("\ninitial loss", mean(los), "\n")
    }
    los <- rep(NA, iter)
    while(d1 > del && k <= iter){
        ctrl$fk <- RET$yhat
        if(sa){
### fk is the previous fitted f
            if(is.null(ctrl$fk)) fk <- 0
            ctrl$s <- quantile(gaussloss(y, ctrl$fk), 0.5)   ### adaptive s,  test program
        }
        RET <- bst(x, y, cost=cost, family=famtype, ctrl = ctrl, control.tree=control.tree, learner=learner)
        los[k] <- mean(loss(y, f=RET$yhat, cost, family = rfamily, s=ctrl$s, sh=ctrl$sh, fk=NULL))
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

bst <- function(x,y, cost=0.5, family = c("gaussian", "hinge", "hinge2", "binom", "expo", "poisson", "tgaussianDC", "thingeDC", "tbinomDC", "binomdDC", "texpoDC", "tpoissonDC", "huber", "thuberDC"), 
                ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree")){
    call <- match.call()
    family <- match.arg(family)
    threshold <- ctrl$threshold
    learner <- match.arg(learner)
    if(learner == "tree" && is.null(colnames(x)))
        colnames(x) <- paste("x", 1:ncol(x), sep = "")
    mstop <- ctrl$mstop
    nu <- ctrl$nu
    twinboost <- ctrl$twinboost
    twintype <- ctrl$twintype
    f.init <- ctrl$f.init
    coefir <- ctrl$coefir
    xselect.init <- ctrl$xselect.init
    center <- ctrl$center
    trace <- ctrl$trace
    numsample <- ctrl$numsample
    df <- ctrl$df
    s <- ctrl$s
    sh <- ctrl$sh
    fk <- ctrl$fk
    maxdepth <- control.tree$maxdepth
    if(twinboost){
        if(is.null(xselect.init))
            stop("Twin boosting requires initial variable selected in the first round\n")
	else if(learner %in% c("sm", "tree") && is.null(f.init))
            stop("Twin boosting requires initial function estimates in the first round\n")
	else if(learner=="ls" && is.null(coefir))
            stop("Twin boosting requires initial coefficients estimates in the first round\n")
    }
    nsample <- dim(x)[1]
    p <- dim(x)[2]
    if(learner == "tree" && p > 10 && twinboost && control.tree$maxdepth >= 2 && is.null(numsample))
        stop("for large p and degree >=2, random sample is suggested\n")   
    oldx <- x 
    one <- rep(1,length(y))
    if(center){
        meanx <- drop(one %*% as.matrix(x))/length(y)
        x <- scale(x, meanx, FALSE) # centers x
    }
    else 
        meanx <- rep(0, length(y))
    if(ctrl$intercept && !ctrl$center) xbar <- matrix(apply(x, 2, mean), nrow=1)
    ens <- vector("list",mstop)
    Fboost <- offset <- 0
    if(family == "gaussian" || family=="tgaussianDC"){
        offset <- mean(y)
    }
    else if(family == "huber" || family == "thuberDC")
        offset <- median(y)
    else if(family %in% c("poisson", "tpoissonDC"))
        offset <- log(mean(y))
    else if(family %in% c("binom", "tbinomDC", "binomdDC", "expo", "texpoDC")){
        ly <- ifelse(y==1, 1-cost, cost)
        offset <- weighted.mean(y > 0, ly)
        offset <- log(offset/(1-offset))
        if(family %in% c("expo", "texpoDC"))
            offset <- 1/2*offset
    }
    Fboost <- rep(offset,length(y))
                                        #if(!twinboost && !is.null(f.init)) ### for nonconvex loss function, may need more work for twin boosting + nonconvex loss 
                                        #Fboost <- f.init
    ystar <- res <- matrix(NA,length(y),ncol(x))
    m <- 1
    coef0 <- sse <- minid <- rep(NA,ncol(x))
    sse <- minid <- rep(NA,ncol(x))
    risk <- coef <- rep(NA, mstop)
    if(ctrl$intercept) int <- rep(0, mstop)
                                        #if(ctrl$intercept && family!="gaussian") int <- rep(0, mstop)
    else int <- NULL
    if(learner=="tree"){
                                        #    maxdepth <- control.tree$maxdepth
        if(maxdepth==1){
            p1 <- ifelse(!twinboost, p,length(xselect.init))
                                        #      xselect <- rep(NA,mstop)
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
                                        #      xselect <- matrix(NA, ncol=2, nrow=mstop)
            xselect.new <- xselect.init[-(length(xselect.init))]
        }
    }
    if(learner=="tree")
        xselect <- vector("list", mstop)
    else xselect<- rep(NA,mstop)
    if(learner=="ls"){
### some calculations independent of mstop and memory allocation
### for each _column_ of the design matrix x, compute the corresponding
### Moore-Penrose inverse (which is a scalar in this case) for the raw
### and standardized input variables
        weights <- rep(1, dim(x)[1]) ###weights to be changed
        xw <- t(x * weights)
        xtx <- colSums(x^2 * weights)
        sxtx <- sqrt(xtx)
### MPinv <- (1 / xtx) * xw
        MPinvS <- (1 / sxtx) * xw
        if(twinboost) MPinvS1 <- abs(coefir) * MPinvS
                                        #see Bulhmann twin boosting paper equation (5)
        if (all(is.na(MPinvS)))
            warning("cannot compute column-wise inverses of design matrix")
    }
                                        #     if(maxdepth == 1) xselect <- rep(NA, mstop)
    while (m <= mstop){
        cor.w <- mse.w <- coef.w <- rep(NA,p)
        u <- ngradient(y, Fboost, cost = cost, family = family, s, sh, fk)
        if(!twinboost) xselect.init <- 1:p
### remove the following line
        cor.w1 <- mse.w1 <- mse.w2 <- mse.w3 <- rep(NA, p)  ### need to remove
        if(learner=="ls"){
### fit least squares to residuals _componentwise_, i.e.,
### compute regression coefficients for each _standardized_
### input variable and select the best variable       
                                        #        xselect <- which.max(abs(mu <- MPinvS %*% u))
### estimate regression coefficient (not standardized)
	    if(twintype == 1){
                mu <- MPinvS %*% u
                if(!twinboost)
                    xselect0 <- which.max(abs(mu))
                else xselect0 <- which.max(abs(MPinvS1 %*% u))
                coef1 <- mu[xselect0] / sxtx[xselect0]
                ind <- xselect0
            }
    else if(twintype==2){
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
        coef1 <- coef0[ind]
    }
	    ml.fit <- lm(u~x[, ind]-1)
	    coef[m] <- coef1 ### note: coef is not the final reported coefficient, see coef.bst
 	    xselect[m] <- ind
	    if(ctrl$intercept){
                if(center) int[m] <- mean(u)
                else {
                    int[m] <- (mean(u) - xbar[,ind] * coef[m]) * nu
                }
            }
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
                    xselect[[m]] <- which(colnames(x) %in% labs)
                }
                else{
                    tree.twin <- vector("list",nrow(inde))  ### check if correct for maxdepth=1
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
                    if(maxdepth==1) xselect[[m]] <- which.max(mse.w)
                    else {
                        tmp <- ml.fit$frame$var[ml.fit$frame$var%in%colnames(x)]
                        tmp <- unique(tmp)
                        if(length(tmp)!=0)
                            xselect[[m]] <- as.character(tmp)
                    }
                }
            }
### update prediction
        if(learner=="sm")
            Fboost <- Fboost + nu * fitted(ml.fit)
        else{
            Fboost <- Fboost + nu * predict(ml.fit)
            if(ctrl$intercept)
                Fboost <- Fboost + int[m]
        }
        if(family %in% c("huber", "tgaussianDC") && (is.null(s) || s < 0)) ### adaptive threshold for Huber loss, but the loss may not decrease along the path
            fk <- Fboost
        if(family %in% c("tgaussianDC", "thingeDC", "tbinomDC", "binomdDC", "texpoDC", "tpoissonDC", "thuberDC") && threshold=="adaptive")
            fk <- Fboost
        risk[m] <- mean(loss(y, Fboost, cost = cost, family = family, s=s, sh=sh, fk=fk))
	if(trace){
            if(m %% 10==0) cat("\nm=", m, "  risk = ", risk[m])
        } 
        ens[[m]] <- ml.fit
        if(m >= 2 && risk[m] > risk[m-1]){
                if(family %in% c("tgaussianDC", "thingeDC", "tbinomDC", "binomdDC", "texpoDC", "tpoissonDC", "thuberDC") && threshold=="standard"){ ###need to check
                    cat("loss value increases at m=", m, "\n") 
                ctrl$mstop <- m
                m <- mstop
            }
    }
        m <- m + 1
    }
    ensemble <- xselect
    xselect <- sort(unique(unlist(xselect)))
    RET <- list(y=y,x=oldx, cost=cost, family = family, learner=learner, yhat=Fboost,offset=offset, ens=ens, control.tree=control.tree, risk=risk, ctrl = ctrl, maxdepth=maxdepth, xselect=xselect, coef = coef, ensemble=ensemble, ml.fit=ml.fit, meanx=meanx, int=int)
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
    s <- object$ctrl$s
    sh <- object$ctrl$sh
    if(is.null(newdata) && is.null(newy))
        ynow <- y
    else ynow <- newy
    if(!missing(newdata)){
        if(object$ctrl$center){
            newdata <- scale(newdata, object$meanx, FALSE) # centers x because the linear predictor below is for centered x as well
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
                                        #    lp <- lp + attributes(coef(object, which=m))$intercept
                                        #if(!is.na(object$intercept))
                                        #if(family!="gaussian" && object$ctrl$intercept)
        if(object$ctrl$intercept)
            lp <- lp + object$int[m]
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
            else if(object$learner=="ls"){
                lp <- lp + nu * object$coef[m] * newdata[, object$ensemble[m]]
            }
        if(type=="all.res")
            res[,m] <- lp
        else if(type=="loss"){
            risk[m] <- mean(loss(ynow, lp, cost = cost, family = family, s=s, sh=sh, fk=lp))
        }
        else if(type == "error"){
            tmp <- sign(lp)
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
    function(x, y, K = 10, cost = 0.5, family = c("gaussian", "hinge", "hinge2", "binom", "expo", "poisson", "tgaussianDC", "thingeDC", "tbinomDC", "binomdDC", "texpoDC", "tpoissonDC"), 
             learner = c("ls", "sm", "tree"), ctrl = bst_control(), type = c("risk", "misc"), plot.it = TRUE, se = TRUE, n.cores=2, ...)
    {
        call <- match.call()
        family <- match.arg(family)
        type <- match.arg(type)
        learner <- match.arg(learner)
        family <- match.arg(family)
        mstop <- ctrl$mstop
        nu <- ctrl$nu
        df <- ctrl$df
        twinboost <- ctrl$twinboost
        twintype <- ctrl$twintype
        trace <- ctrl$trace
        s <- ctrl$s
        sh <- ctrl$sh
        fk <- ctrl$fk
        ctrl.cv <- ctrl
        if(family == "gaussian" && type =="misc") stop("Only risk option is implemented for gaussion family\n")
        all.folds <- cv.folds(length(y), K)
        fraction <- seq(from = 5, to = mstop, by=5)
                                        #fraction <- seq(from = 1, to = mstop, by=5)
        m1 <- 1:length(fraction)
        residmat <- matrix(0, length(fraction), K)
        registerDoParallel(cores=n.cores)
        i <- 1  ###needed to pass R CMD check with parallel code below
	residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
            omit <- all.folds[[i]]
            if(ctrl$twinboost){
                ctrl.cv$f.init <- ctrl$f.init[ - omit]
            }
            fit <- bst(x[ - omit,,drop=FALSE  ], y[ - omit], cost = cost, family = family, learner = learner, ctrl = ctrl.cv, ...)
            fit <- sapply(fraction, function(m) predict(fit, newdata = x[omit,  ,drop=FALSE], mstop = m))
            if(length(omit)==1)fit<-matrix(fit,nrow=1)
### hinge risk or misclassification error
            if(family %in% c("hinge", "hinge2")){
                if(type == "risk"){
                    tmp <- sapply(m1, function(m) loss(y[omit], fit[,m], cost = cost, family = family, s=s, sh=sh, fk=fk))
                    residmat[, i] <- apply(tmp, 2, mean)
                }
            else{
                tmp <- sapply(m1, function(m) {
                    tmp1 <- 0
                    for (ss in 1:length(omit)){
                        if(y[omit[ss]] == 1 && fit[ss,m] <= 0)
			tmp1 <- tmp1 + (1-cost)
                        else if(y[omit[ss]] == -1 && fit[ss,m] > 0)
			    tmp1 <- tmp1 + cost
                    }
                    tmp1 <- tmp1/length(omit)
                }
                )
                residmat[, i] <- 2*tmp  ### multiply 2 to compute misclassification rate if cost=0.5
            }
            }
            else{
                tmp <- sapply(m1, function(m) loss(y[omit], fit[,m], cost = cost, family = family, s=s, sh=sh, fk=fk))
                residmat[, i] <- apply(tmp, 2, mean)
            }
            residmat[, i]   ### return values for parallel computing only 12/9/2015
	}
	cv <- apply(residmat, 1, mean)
        cv.error <- sqrt(apply(residmat, 1, var)/K)
        object<-list(residmat = residmat, mstop = fraction, cv = cv, cv.error = cv.error, family=family)
        if(plot.it){
	         if(type=="risk") ylab <- "Cross-validation loss values"
		      else  if(type=="misc") ylab <- "Cross-validation misclassification errors"
		           plotCVbst(object,se=se, ylab=ylab)
		   }
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
                                        #if(!is.null(x$ctrl$twinboost))
    if(x$ctrl$twinboost)
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
coef.bst <- function(object, which=object$ctrl$mstop, ...) {
    if(object$learner != "ls")
        stop("Coefficients only implemented for linear least squares as base learner\n")
    ret <- numeric(NCOL(object$x))
    xselect <- object$ensemble[1:which]
    for (j in unique(xselect)){
        tmp <- which(xselect == j)
        ret[j] <- sum(object$coef[tmp])
    }
    names(ret) <- colnames(object$x)
    RET <- ret * object$ctrl$nu
    attr(RET, "offset") <- object$offset
    if(object$ctrl$intercept){
        if(object$family=="gaussian" && object$ctrl$center){
            int <- - object$meanx %*% RET
        }
        else{ 
            int <- cumsum(object$int)[which]
            if(object$ctrl$center)
                int <- int - object$meanx %*% RET
        }
        attr(RET, "intercept") <- int
        attr(RET, "offset2int") <- int+object$offset
    }
    RET
}

coefpath.bst <- function(object, ...) {
    if(object$learner != "ls")
        stop("Coefficients only implemented for linear least squares\n")
    vars <- colnames(object$x)
    if(is.null(vars)) stop("colnames of variables missing\n")
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


bst.sel <- function(x, y, q, type=c("firstq", "cv"), ...)
{
    type <- match.arg(type)
    if(is.null(colnames(x)))
        colnames(x) <- paste("x", 1:dim(x)[2], sep="")
    if(type=="cv"){ 
        if(!is.data.frame(x))
            x <- as.data.frame(x)
        fit.cv <- cv.bst(x, y, family="gaussian", plot.it=FALSE, ...)
        fit <- bst(x, y, family="gaussian", ctrl = bst_control(mstop=fit.cv$mstop[which.min(fit.cv$cv)]))
        sel <- coef(fit)
        sel <- which(abs(sel) > 0)
        return(sel)
    }

    fit <- bst(x, y, family="gaussian", ...)
    m <- NULL
    j <- 1
    mstop <- length(fit$ensemble)
    while(length(m) <= q && j <= mstop){
        m <- unique(c(m, fit$ensemble[j]))
        j <- j + 1
    }
    m
}

