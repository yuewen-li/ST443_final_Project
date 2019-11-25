library(tidyverse)
library(lubridate)
library(stats)
library(glmnet)
library(seasonal)
library(FNN)
setwd('C:/Users/mac/Desktop/london-bike-sharing-dataset')
bike_raw <- as.tibble(read.csv('london_merged.csv'))
bike <- bike_raw %>%
  mutate(weather_code = as.factor(weather_code)) %>%
  mutate(is_holiday = as.factor(is_holiday)) %>%
  mutate(is_weekend = as.factor(is_weekend)) %>%
  mutate(timestamp = as.POSIXct(timestamp, format =  "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(date = as.Date(timestamp)) %>%
  mutate(time = as.factor(hour(timestamp))) %>%
  #  mutate(hour = hour(timestamp)) %>%
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
bike$cnt<-log(bike$cnt)
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


Box.test(reg_bike$residuals, lag = 10, type = c( "Ljung-Box"), fitdf = 0)
Box.test(bike_y, lag = 10, type = c( "Ljung-Box"), fitdf = 0)


bikets<-ts(bike$cnt,frequency=12)
ts.plot(bikets)
bikedecomp<-decompose(bikets)
plot(bikedecomp)
bikeadj<-bikets-bikedecomp$seasonal
plot(bikeadj)
Box.test(bikeadj, lag = 10, type = c( "Ljung-Box"), fitdf = 0)

bikestl<-stl(bikets,s.window="period")
plot(bikestl)
bikestlts<-bikestl$time.series
Box.test(bikestlts[,1], lag = 10, type = c( "Ljung-Box"), fitdf = 0)

