
library(MASS)
library(glmnet)
library(tidyverse)
# define a function to generate covariance matrix with dimension p
simu <- function(p,delta,seed=20191118){
  set.seed(seed)
  # create p^2 uniform random numbers
  rand <- runif(p*p)
  # with prob 0.1 to be 0.5, prob 0.9 to be 0
  rand <- ifelse(rand >= .1,0,.5)
  # create p*p matrix
  b <- matrix(rand,p,p)
  # transform b to be symmetric
  b[lower.tri(b)] <- t(b)[lower.tri(b)]
  # set diag to be delta
  diag(b) <- delta
  return(b)
}

# generate n random samples
simu2 <- function(n,p,delta,seed=20191118){
  a <- simu(p,delta,seed=seed)
  # standardise b to have unit diagonals
  theta <- cov2cor(a)
  # solve for the covariance matrix for Gaussian distribution
  sigma <- solve(theta)
  # simulate from multivariate normal distribution
  data <- mvrnorm(n,rep(0,p),sigma)
  return(data)
}

# node-wise lasso
# calc estimated edge matrix
node_wise<-function(n,p,delta,lambda=1,seed=20191118){
  dat <- simu2(n,p,delta,seed=seed)
  coef_<-rep(0,p)
  for (i in 1:p) {
    # the value of lambda to be continued...
    coefficients<-coef(glmnet(dat[,-i],dat[,i], lambda=lambda,intercept = F))
    coefficients<-as.vector(coefficients)
    swap<-coefficients[i]
    coefficients[i]<-coefficients[1]
    coefficients[1]<-swap
    coef_<-rbind(coef_,coefficients)
  }
  coef_ <- coef_[-1,]
  diag(coef_) <- 1
  return(coef_)
}

# transform
# edge<-function(n,p,delta,type,lambda=1,seed=20191118){
#   node<-node_wise(n,p,delta,lambda=lambda,seed=seed)
#   for (i in c(1:p-1)) {
#     for (j in seq(i+1,p)){
#       # specify node-wise lasso 1 or lasso 2
#       if (type=='1'){
#         # node-wise lasso 1: both are non-zero
#         if(isTRUE((node[i,j]!=0) && (node[j,i]!=0)))
#           node[i,j]<-1
#         else
#           node[i,j]<-0
#       }
#       else if (type=='2'){
#         # node-wise lasso 2: either is non-zero
#         if(isTRUE((node[i,j]!=0) || (node[j,i]!=0)))
#           node[i,j]<-1
#         else
#           node[i,j]<-0
#       }
#     }
#   }
#   node[lower.tri(node)] <- t(node)[lower.tri(node)]
#   return(node)
# }

# a more efficient way to transform edge
edge<-function(n,p,delta,type,lambda=1,seed=20191118){
  node<-node_wise(n,p,delta,lambda=lambda,seed=seed)
  upper<-node[upper.tri(node)]
  lower<-node[lower.tri(node)]
  node<-matrix(0,p,p)
  if(type=='1'){
    # node-wise lasso 1: both are non-zero
    upper[(upper!=0)&(lower!=0)]<-1
    node[upper.tri(node, diag=FALSE)] <- upper
    node[lower.tri(node)] <- t(node)[lower.tri(node)]
  }
  if(type=='2'){
    # node-wise lasso 2: either is non-zero
    upper[(upper!=0)|(lower!=0)]<-1
    node[upper.tri(node, diag=FALSE)] <- upper
    node[lower.tri(node)] <- t(node)[lower.tri(node)]
  }
  diag(node)<-1
  return(node)
}

res1=edge(1000,10,4,type='1')
res2=edge(1000,10,4,type='2')
original<-simu(10,4)

library(glasso)
# calculate edge of graphic lasso model
graphic <- function(n,p,delta,rho,seed=20191118){
  dat <- simu2(n,p,delta,seed=seed)
  cov <- cov(dat)
  glasso <- glasso(cov,rho)
  wi <- glasso$wi
  wi[!wi==0] <- 1
  return(wi)
}
res3 <- graphic(1000,10,4,0.1)

# compute TPR and FPR
# input for data and reference is matrix
tfpr <- function(data, reference){
  dim <- ncol(data)
  data[!data==0] <- 1
  reference[!reference==0] <- 1
  # flatten the matrix to vector
  data_flat <- as.vector(data)
  ref_flat <- as.vector(reference)
  # extract the diagnal elements
  diag_index <- (dim)*seq(0,dim-1)+seq(1,dim)
  data_flat <- data_flat[-diag_index]
  ref_flat <- ref_flat[-diag_index]
  # calculate TPR, TPR=TP/P
  tpr <- length(ref_flat[(data_flat==1)&(ref_flat==1)])/length(ref_flat[ref_flat==1])
  # calculate FPR, FPR=FP/N
  fpr <- length(ref_flat[(data_flat==1)&(ref_flat==0)])/length(ref_flat[ref_flat==0])
  # calculate overall prediction accuracy
  overall <- length(ref_flat[data_flat==ref_flat])/length(ref_flat)
  # calculate precision
  ppv <- length(ref_flat[(data_flat==1)&(ref_flat==1)])/length(data_flat[data_flat==1])
  # calculate F1 score
  f1 <- 2*ppv*tpr/(ppv+tpr)
  out <- c(tpr=tpr,fpr=fpr,overall=overall,ppv=ppv,f1=f1)
  return(out)
}

(out1 <- tfpr(res1,original))
(out2 <- tfpr(res2,original))
(out3 <- tfpr(res3,original))
par(mfrow=c(1,3))

# calc auc
calc_auc <- function(data){
  data <- data[order(data[,2]),]
  auc <- 0
  for(i in 2:nrow(data)){
    auc = auc + (data[i,2]-data[i-1,2]) * data[i,1]
  }
  return(auc)
}

# plot roc curve
# the range of lambda should be reconsidered
roc <- function(n,p,delta,type,shrink,k=100,seed=20191118,plot=T){
  roc_curve <- matrix(NA, k+2, 6)
  ref <- simu(p,delta,seed=seed)
  for (i in seq(1:k)){
    lambda <- i^3/shrink
    if(type=='g')
      res <- graphic(n,p,delta,rho = lambda,seed=seed)
    else
      res <- edge(n,p,delta,lambda = lambda,type=type,seed=seed)
    roc_curve[i,1:5] <- tfpr(res,ref)
    roc_curve[i,6] <- lambda
  }
  roc_curve[k+1,] <- c(1,1,NA,NA,NA,0)
  roc_curve[k+2,] <- c(0,0,NA,NA,NA,Inf)
  roc_curve <- as.data.frame(roc_curve) %>%
    arrange(V2,V1)
  colnames(roc_curve) <- c('tpr','fpr','overall','ppv','f1','lambda')
  auc <- calc_auc(roc_curve[,1:2])
  if (plot == T){
    plot(roc_curve[,2],roc_curve[,1],type='s',xlab = 'fpr',ylab = 'tpr',main = 'ROC curve',xlim = c(0,1),ylim = c(0,1))
    mtext(paste('auc =', round(auc,3)),3)
  }
  roc_curve <- roc_curve %>%
    arrange(lambda)
  return(list(value=roc_curve,auc=auc))
}

test11 <- roc(1000,10,4,type='g',10000,seed=20191118)
test12 <- roc(1000,10,4,type='1',10000,seed=20191120)
test13 <- roc(1000,10,4,type='2',10000,seed=20191120)


# try n=p=100
ref <- simu(50,4)
test21 <- roc(50,50,4,type = 'g',10000)
test22 <- roc(50,50,4,type = '1',10000)
test23 <- roc(50,50,4,type = '2',10000)

# try n=60, p=100
ref <- simu(40,4)
test31 <- roc(30,40,4,type = 'g',4000000,k=1500)
test32 <- roc(30,40,4,type = '1',100000,k=300)
test33 <- roc(30,40,4,type = '2',100000,k=300)

# find the best tuning parameter for each model
# select tuning parameter based on accuracy rate and F1 score
tunning <- function(list){
  value <- list$value
  value <- value[-1,]
  value <- value[-nrow(value),]
  lambda_overall <- value[which.max(value$overall),]
  lambda_f1 <- value[which.max(value$f1),]
  return(list(overall=lambda_overall,f1=lambda_f1))
}
(lambda11 <- tunning(test11))
(lambda12 <- tunning(test12))
(lambda13 <- tunning(test13))
(lambda21 <- tunning(test21))
(lambda22 <- tunning(test22))
(lambda23 <- tunning(test23))
(lambda31 <- tunning(test31))
(lambda32 <- tunning(test32))
(lambda33 <- tunning(test33))

# replicate for 50 times
rep_50 <- function(n,p,delta,type,shrink,k=100){
  overall <- data.frame(tpr=NA,fpr=NA,overall=NA,ppv=NA,f1=NA,lambda=NA)
  f1 <- data.frame(tpr=NA,fpr=NA,overall=NA,ppv=NA,f1=NA,lambda=NA)
  auc <- vector()
  for (i in seq(1:50)){
    original<-simu(p,delta,seed=20191117+i)
    value <- roc(n,p,delta,type,shrink,k=k,plot=F,seed=20191117+i)
    lambda <- tunning(value)
    overall[i,] <- lambda$overall
    f1[i,] <- lambda$f1
    auc[i] <- value$auc
  }
  return(list(overall=overall,f1=f1,auc=auc))
}

# try n=1000,p=10
rep11 <- rep_50(1000,10,4,type='g',10000)
rep12 <- rep_50(1000,10,4,type='1',10000)
rep13 <- rep_50(1000,10,4,type='2',10000)
# save(rep11,rep12,rep13,file='rep1.RData')
load('rep1.RData')

auc1=data.frame(g=rep11$auc,nw1=rep12$auc,nw2=rep13$auc)
boxplot(auc1,main='Boxplot of AUC values under three approaches')

fpr1=data.frame(g=rep11$overall['fpr'],nw1=rep12$overall['fpr'],nw2=rep13$overall['fpr'])
colnames(fpr1)<-c('g','nw1','nw2')
boxplot(fpr1,main='Boxplot of fpr values under three approaches')

lambda1<-data.frame(g=rep11$overall['lambda'],nw1=rep12$overall['lambda'],nw2=rep13$overall['lambda'])
summary(lambda1)

# try n=p=50
rep21 <- rep_50(50,50,4,type = 'g',10000)
rep22 <- rep_50(50,50,4,type = '1',10000)
rep23 <- rep_50(50,50,4,type = '2',10000)
# save(rep21,rep22,rep23,file='rep2.RData')
load('rep2.RData')

auc2=data.frame(g=rep21$auc,nw1=rep22$auc,nw2=rep23$auc)
boxplot(auc2,main='Boxplot of AUC values under three approaches')

fpr2=data.frame(g=rep21$overall['fpr'],nw1=rep22$overall['fpr'],nw2=rep23$overall['fpr'])
colnames(fpr2)<-c('g','nw1','nw2')
boxplot(fpr2,main='Boxplot of fpr values under three approaches')

lambda2<-data.frame(g=rep21$overall['lambda'],nw1=rep22$overall['lambda'],nw2=rep23$overall['lambda'])
summary(lambda2)


# try n=30, p=40
#rep31 <- rep_50(30,40,4,type = 'g',4000000,k=1500)
#rep32 <- rep_50(30,40,4,type = '1',100000,k=300)
#rep33 <- rep_50(30,40,4,type = '2',100000,k=300)
# save(rep31,rep32,rep33,file='rep3.RData')
# load('rep3.RData')

#try n=50, p=100
rep31 <- rep_50(50,100,4,type = 'g',4000000,k=1500)
rep32 <- rep_50(50,100,4,type = '1',100000,k=300)
rep33 <- rep_50(50,100,4,type = '2',100000,k=300)
save(rep31,rep32,rep33,file='rep31.RData')
load('rep31.RData')

auc3=data.frame(g=rep31$auc,nw1=rep32$auc,nw2=rep33$auc)
boxplot(auc3,main='Boxplot of AUC values under three approaches')

fpr3=data.frame(g=rep31$overall['fpr'],nw1=rep32$overall['fpr'],nw2=rep33$overall['fpr'])
colnames(fpr3)<-c('g','nw1','nw2')
boxplot(fpr3,main='Boxplot of fpr values under three approaches')

lambda3<-data.frame(g=rep31$overall['lambda'],nw1=rep32$overall['lambda'],nw2=rep33$overall['lambda'])
summary(lambda3)


# compare the time spend between models
library(microbenchmark)
# for n=1000, p=10
microbenchmark(edge(1000,10,4,type='1'))
microbenchmark(edge(1000,10,4,type='2'))
microbenchmark(graphic(1000,10,4,0.1))
# for n=50,p=50
microbenchmark(edge(50,50,4,type='1'))
microbenchmark(edge(50,50,4,type='2'))
microbenchmark(graphic(50,50,4,0.1))
# for n=30,p=40
microbenchmark(edge(30,40,4,type='1'))
microbenchmark(edge(30,40,4,type='2'))
microbenchmark(graphic(30,40,4,0.1))
# for n=50,p=100
microbenchmark(edge(50,100,4,type='1'))
microbenchmark(edge(50,100,4,type='2'))
microbenchmark(graphic(50,100,4,0.1))
# easily see that graphical lasso is far more efficient

# check for assymmetric of wi
b <- graphic(30,40,4,1/100000)
sum(abs((b$wi[lower.tri(b$wi)]-t(b$wi)[lower.tri(t(b$wi))])),na.rm=T)
# try different tunning parameter rho
dis <- vector()
rho <- vector()
for (i in seq(1,100)){
  rho[i] <- i/10000
  b <- graphic(30,40,4,i/10000)
  dis[i] <- sum(abs((b$wi[lower.tri(b$wi)]-t(b$wi)[lower.tri(t(b$wi))])),na.rm=T)
}
plot(rho,dis,type='l',main='degree of asymmetry',xlab='rho',ylab='distance')
############################################################################
