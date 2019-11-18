library(tidyverse)
library(lubridate)
library(stats)
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
## linear regression
summary(reg_bike <- lm(cnt~.+I(trend^2), data = bike))
summary(fit <- fitted(reg_bike))
summary((fit-bike_y)^2)
## lasso
bike_x <- model.matrix(cnt~.-1 +I(trend^2), bike)
bike_y <- bike %>%
  select(cnt)
lasso_cv <- cv.glmnet(bike_x, bike_y)
coef(lasso <- glmnet(bike_x, bike_y, lambda = lasso_cv$lambda.1se))
fit_lasso <- predict(lasso, bike_x)
summary(fit_lasso-bike_y)
#######################
bike_month <- bike %>%
  mutate(month = month(date)) %>%
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
