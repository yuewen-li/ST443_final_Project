library(tidyverse)
library(lubridate)
library(stats)
library(glmnet)
library(FNN)
setwd('D:\\LSE\\ST443 Machine Learning and Data Mining\\pj')
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
  select(-c(timestamp, date)) %>%
  arrange(t1, t2) %>%
  filter(is.na(month) == F) %>%
  mutate(trend = seq(1, 17413)) %>%
  mutate(sin = sin(trend/(24*365))) %>%
  mutate(cos = cos(trend/(24*365))) 
bike
bike_y <- bike$cnt
## devide testing set and training set
set.seed(443)
index <- seq(1, 17413)
test_index <- sample(index, 3483)
train_index <- index[-test_index]
## linear regression
summary(reg_bike <- lm(cnt~.+I(trend^2), data = bike[train_index,]))
yhat_lm <- predict(reg_bike, bike[test_index,])
summary((yhat_lm - bike$cnt[test_index])^2)
## lasso, first do the variables selection
bike_x <- model.matrix(cnt~. +I(trend^2), bike[train_index,])[,-1]
lasso_cv <- cv.glmnet(bike_x, bike_y[train_index])
coef(lasso_bike <- glmnet(bike_x, bike_y[train_index], lambda = lasso_cv$lambda.1se))
bike_x_test <- model.matrix(cnt~. +I(trend^2), bike[test_index,])[,-1]
yhat_lasso <- predict(lasso_bike, bike_x_test)
summary((yhat_lasso - bike$cnt[test_index])^2)
## knn regression
yhat_knn <- knn.reg(bike_x, bike_x_test, bike_y, k = 10)
summary((yhat_knn$pred - bike$cnt[test_index])^2)
#######################
bike_month <- bike %>%
  group_by(month) %>%
  mutate(monthuse = mean(cnt)) %>%
  select(month, monthuse) %>%
  unique()
plot(bike_month$month, bike_month$monthuse)

bike_hour <- bike %>%
  group_by(time) %>%
  mutate(hourave = mean(cnt)) %>%
  select(hourave, time) %>%
  unique()
plot(bike_hour$time, bike_hour$hourave, type = 'p')

bike_day <- bike %>%
  group_by(day) %>%
  mutate(dayavg = mean(cnt)) %>%
  select(day, dayavg)
plot(bike_day$day, bike_day$dayavg)  
