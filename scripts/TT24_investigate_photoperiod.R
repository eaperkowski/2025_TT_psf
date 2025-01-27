# Libraries
library(tidyverse)
library(bigleaf)
library(lubridate)

# Read hourly weather .csv file from 2024 season
weather <- read.csv("../data/TT24_hourly_weather_data.csv") %>%
  dplyr::select(date, doy, solar_radiation_wm2) %>%
  separate(date, into = c("date_only", "time_only"), remove = FALSE, sep = " ") %>%
  mutate(date_only = ymd(date_only),
         solar_radiation_umol = Rg.to.PPFD(Rg = solar_radiation_wm2))
head(weather)

# Determine photoperiod for each day
light_data <- weather %>%
  filter(date_only > "2024-04-11" & date_only < "2024-09-15" & 
           solar_radiation_umol > 0) %>%
  group_by(date_only) %>%
  summarize(daylight_hours = length(solar_radiation_umol))
  
light_data_with_hourly_par <- light_data %>%
    left_join(weather) %>%
  dplyr::select(date, date_only, time_only, doy, solar_radiation_umol, daylight_hours)

# Some quick visualizations
## Photoperiod
ggplot(data = light_data, aes(x = date_only, y = daylight_hours)) +
  geom_line()

## PAR
ggplot(data = light_data_with_hourly_par, aes(x = as.Date(date), y = solar_radiation_umol)) +
  geom_line()

# What should photoperiod and max PAR be for biweekly growth chamber timesteps?

# April 16 - April 30, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-04-16" & date_only <= "2024-04-30") %>%
  summarize(photoperiod = mean(daylight_hours))
## 14 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-04-16" & date_only <= "2024-04-30") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 653 umol/m2/s


# May 1 - May 15, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-05-01" & date_only <= "2024-05-15") %>%
  summarize(photoperiod = mean(daylight_hours))
## 13.2 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-05-01" & date_only <= "2024-05-15") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 653 umol/m2/s

# May 16 - May 31, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-05-16" & date_only <= "2024-05-31") %>%
  summarize(photoperiod = mean(daylight_hours))
## 14 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-05-16" & date_only <= "2024-05-31") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 653 umol/m2/s

# June 1 - June 15, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-06-01" & date_only <= "2024-06-15") %>%
  summarize(photoperiod = mean(daylight_hours))
## 13.1 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-06-01" & date_only <= "2024-06-15") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 653 umol/m2/s

# June 16 - June 30, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-06-16" & date_only <= "2024-06-30") %>%
  summarize(photoperiod = mean(daylight_hours))
## 12.9 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-06-16" & date_only <= "2024-06-30") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 653 umol/m2/s

# July 1 - July 15, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-07-01" & date_only <= "2024-07-15") %>%
  summarize(photoperiod = mean(daylight_hours))
## 14 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-07-01" & date_only <= "2024-07-15") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 250 umol/m2/s

# July 16 - July 31, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-07-16" & date_only <= "2024-07-31") %>%
  summarize(photoperiod = mean(daylight_hours))
## 13.3 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-07-16" & date_only <= "2024-07-31") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 265 umol/m2/s

# August 1 - August 15, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-08-01" & date_only <= "2024-08-15") %>%
  summarize(photoperiod = mean(daylight_hours))
## 12.6 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-08-01" & date_only <= "2024-08-15") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## 280 umol/m2/s

# August 16 - August 31, 2024
## Photoperiod
light_data %>%
  filter(date_only >= "2024-08-16" & date_only <= "2024-08-31") %>%
  summarize(photoperiod = mean(daylight_hours))
## 12 hour photoperiod

## PAR
light_data_with_hourly_par %>%
  filter(date_only >= "2024-08-16" & date_only <= "2024-08-31") %>%
  group_by(date_only) %>%
  summarize(max_par = max(solar_radiation_umol)) %>%
  ungroup() %>%
  summarize(mean_max_par = mean(max_par))
## PAR = 174 umol/m2/s






