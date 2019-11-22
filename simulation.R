library(MASS)
library(glmnet)

# define a function to generate covariance matrix with dimension p
simu <- function(p,delta){
  set.seed(20191118)
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
simu2 <- function(n,p,delta){
  a <- simu(p,delta)
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
node_wise<-function(n,p,delta){
  dat <- simu2(n,p,delta)
  coef_<-rep(0,p)
  for (i in 1:p) {
    fit_node <-cv.glmnet(dat[,-i],dat[,i])
    
    # the value of lambda to be continued...
    coefficients<-coef(glmnet(dat[,-i],dat[,i], lambda=fit_node$lambda.min))
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
edge<-function(n,p,delta){
  node<-node_wise(n,p,delta)
  for (i in c(1:p-1)) {
    for (j in seq(i+1,p)){
      # node-wise lasso 1: both are non-zero
      if(isTRUE((node[i,j]!=0) && (node[j,i]!=0)))
      # node-wise lasso 2: either is non-zero
      # if(isTRUE((node[i,j]!=0) || (node[j,i]!=0)))
        node[i,j]<-1
      else
        node[i,j]<-0
    }
  }
  node[lower.tri(node)] <- t(node)[lower.tri(node)]
  return(node)
}
res=edge(1000,10,2)
original<-simu(10,2)
