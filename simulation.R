library(MASS)

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
  # standardise b to have unit diagonals
  theta <- cor(b)
  # solve for the covariance matrix for Gaussian distribution
  sigma <- solve(theta)
  return(sigma)
}

a=simu(5,1e-10)
# simulate from multivariate normal distribution
x <- mvrnorm(1000,rep(0,5),a)