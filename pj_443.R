library(tidyverse)
library(lubridate)
setwd('C:\\Users\\wangm\\Desktop\\')
bike_raw <- as.tibble(read.csv('london_merged.csv'))
bike <- bike_raw %>%
  mutate(weather_code = as.factor(weather_code)) %>%
  mutate(is_holiday = as.factor(is_holiday)) %>%
  mutate(is_weekend = as.factor(is_weekend)) %>%
  mutate(timestamp = as.POSIXct(timestamp, format =  "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(date = as.Date(timestamp)) %>%
  mutate(time = hour(timestamp)) %>%
  select(-timestamp)
bike <- bike[,c(10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9)]
df <- cbind(df[ncol(df)], df[1:(ncol(df) - 1)])
bike
summary(reg_bike <- lm(cnt~.-timestamp, data = bike1))
