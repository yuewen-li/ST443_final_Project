library(tidyverse)
library(lubridate)
library(stats)
library(glmnet)
library(FNN)
library(MASS)
library(caret)
library(dplyr)
library(tree)
library(randomForest)
library(gbm)

bike_raw <- as.tibble(read.csv('london_merged.csv'))

bike <- bike_raw %>%
  mutate(weather_code = as.factor(weather_code)) %>%
  mutate(is_holiday = as.factor(is_holiday)) %>%
  mutate(is_weekend = as.factor(is_weekend)) %>%
  mutate(timestamp = as.POSIXct(timestamp, format =  "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(date = as.Date(timestamp)) %>%
  mutate(time = as.factor(hour(timestamp))) %>%
  mutate(day = day(date)) %>%
  mutate(month = month(date)) %>%
  mutate(year = year(date)) %>%
  arrange(t1, t2) %>%
  filter(is.na(month) == F) %>%
  mutate(yday = yday(date)) %>%
  mutate(sin = sin(2*pi*yday/(365))) %>%
  mutate(cos = cos(2*pi*yday/(365))) %>%
  dplyr::select(-c(timestamp,yday,year,month,day,date))
bike
bike_y <- bike$cnt
## devide testing set and training set
set.seed(443)
index <- seq(1, 17413)
test_index <- sample(index, 3483)
train_index <- index[-test_index]
bike_y_test = bike_y[test_index]
bike_y = bike_y[-test_index]

## linear regression
summary(reg_bike <- lm(cnt~., data = bike[train_index,]))
yhat_lm <- predict(reg_bike, bike[test_index,])
summary((yhat_lm - bike$cnt[test_index])^2) #MSE 315967
par(mfrow=c(2,2))
plot(reg_bike)

## Ridge regression
bike_x <- model.matrix(cnt~., bike[train_index,])[,-1] #remove B0
fit.ridge <-glmnet(bike_x,bike_y , alpha=0)
plot(fit.ridge, xvar="lambda", label= TRUE)
plot(fit.ridge, xvar="dev", label= TRUE)
cv.ridge <-cv.glmnet(bike_x,bike_y, alpha=0)
plot(cv.ridge)
coef(cv.ridge)
coef(glmnet(bike_x,bike_y ,alpha=0, lambda=cv.ridge$lambda.min))
opt_lambda <- cv.ridge$lambda.min
opt_lambda
fit <- cv.ridge$glmnet.fit
summary(fit)
bike_x_test <- model.matrix(cnt~., bike[test_index,])[,-1]
yhat_ridge <- predict(cv.ridge, bike_x_test,type='response')
summary((yhat_ridge - bike_y_test)^2)



## lasso, first do the variables selection
bike_x <- model.matrix(cnt~., bike[train_index,])[,-1] #remove B0
bike_x
lasso_cv <- cv.glmnet(bike_x, bike_y,family=c('poisson'))
lasso_cv
coef(lasso_bike <- glmnet(bike_x, bike_y, lambda = lasso_cv$lambda.1se,family=c('poisson')))
bike_x_test <- model.matrix(cnt~., bike[test_index,])[,-1]
yhat_lasso <- predict(lasso_bike,bike_x_test,type='response')
summary((yhat_lasso - bike$cnt[test_index])^2)
coef(lasso_cv)
coef(glmnet(bike_x,bike_y, lambda=lasso_cv$lambda.min))

## Validation set approach to select best lambda in Lasso
set.seed(1)
train <-sample(nrow(bike_x), 3483, replace=FALSE)
lasso.train <-glmnet(bike_x[train,], bike_y[train])
pred.test <-predict(lasso.train, bike_x[-train,])
rmse <-sqrt(apply((bike_y[-train]-pred.test)^2,2,mean))
plot(log(lasso.train$lambda), rmse, type="b", xlab="Log(lambda)")
lambda.best <-lasso.train$lambda[order(rmse)[1]]
lambda.best

# poisson regression
count.fits <- function(model,newdata,y=NULL){
  fit.mu <- newdata
  fit.mu$muhat <- predict(model,newdata,type="response")
  if(!is.null(y)){
    ind.tmp <- seq(ncol(fit.mu)+1,ncol(fit.mu)+length(y))
    fit.mu <- cbind(fit.mu,matrix(NA,nrow(fit.mu),length(y)))
    for(j in seq(nrow(fit.mu))){
      if(class(model)[1]=="glm")
        fit.mu[j,ind.tmp] <- dpois(y,lambda=fit.mu$muhat[j])
      if(class(model)[1]=="negbin"){
        fit.mu[j,ind.tmp] <- dnbinom(y,size=model$theta,mu=fit.mu$muhat[j])  
      }  
    }
    colnames(fit.mu)[ind.tmp] <- paste("P(Y=",y,")",sep="")
  }
  fit.mu
}
summary(pois_bike <- glm(cnt~., data=bike[train_index,],family='poisson'))
yhat <- count.fits(pois_bike,bike[test_index,])
summary((yhat$muhat - bike$cnt[test_index])^2)
plot(pois_bike)

## Negative binomial regression
neg_bike <- glm.nb(cnt~., data=bike[train_index,])
summary(neg_bike)
Yhat_neg = predict(neg_bike, newdata=bike[test_index,],type="response")
mean((bike_y_test-Yhat_neg)^2)

## knn regression
# KNN Plot model accuracy vs different values of k=5
set.seed(123)
model <- train(
  cnt ~., data = bike[train_index,], method = "knn",
  trControl = trainControl("cv", number = 10),
  preProcess = c("center","scale"),
  tuneLength = 20
)
plot(model)
model$bestTune

yhat_knn <- knn.reg(bike_x, bike_x_test, bike_y, k = 5)
summary((yhat_knn$pred - bike$cnt[test_index])^2) #high mse



#######################
bike_month <- bike %>%
  group_by(month) %>%
  mutate(monthuse = mean(cnt)) %>%
  select(month, monthuse) %>%
  unique()
plot(bike_month$month, bike_month$monthuse)

bike_hour <- bike %>%
  group_by(time,is_holiday) %>%
  mutate(hourave = mean(cnt)) %>%
  select(hourave, time,is_holiday) %>%
  group_by() %>%
  mutate(time=as.numeric(time)) %>%
  unique()
plot(bike_hour$time, bike_hour$hourave, type = 'p')

bike_day <- bike %>%
  group_by(day) %>%
  mutate(dayavg = mean(cnt)) %>%
  select(day, dayavg)
plot(bike_day$day, bike_day$dayavg)  
ggplot(bike_hour,aes(time,hourave,color=is_holiday))+
  geom_line()

##Regressiontree
attach(bike)
tree.bike<-tree(cnt~. , data=bike, subset=train_index)
summary(tree.bike)
plot(tree.bike)
text(tree.bike, pretty=0)

cv.bike <- cv.tree(tree.bike)
plot(cv.bike$size, cv.bike$dev, type = "b")

#best=6

#Prune
prune.bike <- prune.tree(tree.bike, best = 5)
plot(prune.bike)
text(prune.bike, pretty = 0)

yhat <- predict(tree.bike, newdata = bike[-train_index,])
bike.test <- bike[-train_index, "cnt"]
plot(yhat, bike[-train_index,]$cnt)
abline(0, 1)

mean((yhat - bike[-train_index,]$cnt) ^ 2)



###########
# Bagging #
###########
bag.bike<-randomForest(cnt~. , data= bike, subset=train_index, mtry=11, importance=TRUE)
bag.bike

yhat.bag <-predict(bag.bike, newdata = bike[-train_index,])
plot(yhat.bag, bike[-train_index,]$cnt)
abline(0,1)

mean((yhat.bag-bike[-train_index,]$cnt) ^2 )

importance(bag.bike)
varImpPlot(bag.bike)

# time and is_weekend is most important

############
# Boosting #
############

boost.bike <- gbm( cnt ~ ., data = bike[train_index,], distribution = "gaussian", n.trees = 5000, interaction.depth = 4)

summary(boost.bike)
par(mfrow = c(1, 2))
plot(boost.bike, i="time")
plot(boost.bike, i='is_weekend')

yhat.boost = predict(boost.bike, newdata = bike[-train_index, ], n.trees = 5000)
mean((yhat.boost - bike[-train_index,]$cnt) ^ 2)
