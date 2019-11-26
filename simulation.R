
library(MASS)
library(glmnet)

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
edge<-function(n,p,delta,type,lambda=1,seed=20191118){
  node<-node_wise(n,p,delta,lambda=lambda,seed=seed)
  for (i in c(1:p-1)) {
    for (j in seq(i+1,p)){
      # specify node-wise lasso 1 or lasso 2
      if (type=='1'){
      # node-wise lasso 1: both are non-zero
      if(isTRUE((node[i,j]!=0) && (node[j,i]!=0)))
        node[i,j]<-1
      else
        node[i,j]<-0
      }
      else if (type=='2'){
      # node-wise lasso 2: either is non-zero
      if(isTRUE((node[i,j]!=0) || (node[j,i]!=0)))
        node[i,j]<-1
      else
        node[i,j]<-0
      }
    }
  }
  node[lower.tri(node)] <- t(node)[lower.tri(node)]
  return(node)
}
res1=edge(1000,10,2,type='1')
res2=edge(1000,10,2,type='2')
original<-simu(10,2)

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
res3 <- graphic(1000,10,2,0.1)

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
  out <- c(tpr=tpr,fpr=fpr)
  return(out)
}
ref <- simu(10,2)
(out1 <- tfpr(res1,ref))
(out2 <- tfpr(res2,ref))
(out3 <- tfpr(res3,ref))
# plot roc curve for nodewise lasso1
roc1 <- matrix(NA, 100, 2)
for (i in seq(1:100)){
  lambda <- i/1000
  res1 <- edge(1000,10,2,lambda = lambda,type='1')
  roc1[i,] <- tfpr(res1,ref)
}
plot(roc1[,2],roc1[,1],type='l')
# for nodewise lasso2
roc2 <- matrix(NA, 100, 2)
for (i in seq(1:100)){
  lambda <- i/100
  res2 <- edge(1000,10,2,lambda = lambda,type='2')
  roc2[i,] <- tfpr(res2,ref)
}
plot(roc2[,2],roc2[,1],type='l')
# for graphic lasso
roc3 <- matrix(NA, 1000, 2)
for (i in seq(1:1000)){
  lambda <- i/1000
  res3 <- graphic(1000,10,2,rho = lambda)
  roc3[i,] <- tfpr(res3,ref)
}
plot(roc3[,2],roc3[,1],type='l')
