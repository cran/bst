### R code from vignette source 'khan.Rnw'

###################################################
### code chunk number 1: khan.Rnw:45-46
###################################################
options(prompt = "R> ", continue = " ", width = 70, digits =4, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: khan.Rnw:48-94
###################################################
library("bst")
datafile <- system.file("extdata", "supplemental_data", package="bst")
dat0 <- read.delim(datafile, header=TRUE, skip=1)[,-(1:2)]
genename <- read.delim(datafile, header=TRUE, skip=1)[,(1:2)]
dat0 <- t(dat0)
dat1 <- dat0[rownames(dat0) %in% 
c("TEST.9", "TEST.13","TEST.5", "TEST.3", "TEST.11"),]
dat2 <- dat0[!rownames(dat0) %in% 
c("TEST.9", "TEST.13","TEST.5", "TEST.3", "TEST.11"),]
dat2 <- rbind(dat2, dat1)
train <- dat2[1:63,] ### training samples
test <- dat2[64:83,] ### test samples
train.classes <- substr(rownames(train), 1,2)
test.classes <- c("NB","RM","NB","EW","RM","BL","EW","RM","EW","EW","EW",
"RM","BL","RM","NB","NB","NB","NB","BL","EW")
train.classes <- as.numeric(factor(train.classes, levels=c("EW", "BL", "NB", "RM")))
test.classes <- as.numeric(factor(test.classes, levels=c("EW", "BL", "NB", "RM")))
### pre-processing training data: standardize predictors after log-transformation
train <- log10(train)
x <- train
meanx <- colMeans(x)
one <- rep(1,nrow(x))
normx <- sqrt(drop(one %*% (x^2)))
train <- scale(train, meanx, normx)
### compute a marginal relevance measure
tmp <- cbind(train, train.classes)
a0 <- b0 <- 0
  for(k in 1:length(table(train.classes))){
    tmp1 <- subset(tmp, tmp[,2309]==k)
    xc.bar <- colMeans(tmp1[,-2309])   ###average of gene j across class k 
    xa.bar <- colMeans(tmp[,-2309])    ###average of gene j across all samples
    a0 <- a0 + dim(tmp1)[1] * ((xc.bar - xa.bar)^2)
    b0 <- b0 + colSums((tmp[,-2309] - xc.bar)^2)
}
bw <- a0/b0
### select top 300 genes based on the ordered marginal relevance measure
npre <- 300
bw1 <- order(bw, decreasing=TRUE)[1:npre]
train <- train[,bw1]
### pre-processing test data: standardize predictors after log-transformation
### select the same 300 genes as in the training data
test <- log10(test)
test <- scale(test, meanx, normx)[, bw1]
test <- as.data.frame(test)
colnames(train) <- paste("x", 1:dim(train)[2], sep="")
colnames(test) <- paste("x", 1:dim(test)[2], sep="")


###################################################
### code chunk number 3: khan.Rnw:99-103
###################################################
m <- 30
set.seed(123)
dat.cvm <- cv.mhingebst(x=train, y=train.classes, balance=TRUE, K=5, 
ctrl = bst_control(mstop=m), family = "hinge", learner = "sm", type="misc", n.cores=2)


###################################################
### code chunk number 4: khan.Rnw:108-113
###################################################
m1 <- 20
dat.m1 <- mhingebst(x=train, y=train.classes, ctrl = bst_control(mstop=m1),
 family = "hinge", learner = "sm")
risk.te1 <- predict(dat.m1, newdata=test, newy=test.classes, type="error")
plot(risk.te1, type="l", xlab="Iteration", ylab="Test Error")


###################################################
### code chunk number 5: khan.Rnw:117-118
###################################################
plot(nsel(dat.m1, m1), ylab="No. Genes", xlab="Iteration", lty="solid", type="l")


###################################################
### code chunk number 6: khan.Rnw:123-130
###################################################
m2 <- 20
xinit <- unlist(dat.m1$xselect)
xinit <- subset(xinit, !is.na(xinit))
dat.m2 <- mhingebst(x=train, y=train.classes, family = "hinge", learner = "sm",
ctrl = bst_control(mstop=m2, twinboost=TRUE, f.init=dat.m1$yhat, xselect.init=xinit)) 
risk.te2 <- predict(dat.m2,newdata=test,newy=test.classes,type="error")
plot(risk.te2, type="l", xlab="Iteration", ylab="Test Error")


###################################################
### code chunk number 7: khan.Rnw:134-135
###################################################
plot(nsel(dat.m2, m2), ylab="No. Genes", xlab="Iteration", lty="solid", type="l")


