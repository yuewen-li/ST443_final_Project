library(tidyverse)
library(lubridate)
library(stats)
library(glmnet)
library(FNN)

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
  mutate(sin = sin(2*pi*yday/(24*365))) %>%
  mutate(cos = cos(2*pi*yday/(24*365))) %>%
  select(-c(timestamp,yday,year,month,day,date))
bike
bike_y <- bike$cnt
## devide testing set and training set
set.seed(443)
index <- seq(1, 17413)
test_index <- sample(index, 3483)
train_index <- index[-test_index]

## linear regression
summary(reg_bike <- lm(cnt~., data = bike[train_index,]))
yhat_lm <- predict(reg_bike, bike[test_index,])
summary((yhat_lm - bike$cnt[test_index])^2) #MSE 315967
par(mfrow=c(2,2))
plot(reg_bike)

## lasso, first do the variables selection
bike_x <- model.matrix(cnt~., bike[train_index,])[,-1] #remove B0
bike_x
lasso_cv <- cv.glmnet(bike_x, bike_y[train_index],family=c('poisson'))
lasso_cv
coef(lasso_bike <- glmnet(bike_x, bike_y[train_index], lambda = lasso_cv$lambda.1se,family=c('poisson')))
bike_x_test <- model.matrix(cnt~., bike[test_index,])[,-1]
yhat_lasso <- predict(lasso_bike,bike_x_test,type='response')
summary((yhat_lasso - bike$cnt[test_index])^2)

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

## knn regression
yhat_knn <- knn.reg(bike_x, bike_x_test, bike_y, k = 10)
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
