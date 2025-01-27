## Trillium Trail Integrative Bio
## TOMST sensor data (air & soil temperatures, soil moisture)
## Code written and annotated by JM 01-2025

# Set working directory to folder where TOMST data are saved.
setwd("../2024_tomst_probe_data/")

# I'm gonna plot some stuff using ggplot.
require(ggplot2)

## =============================================================================
## (1) Read in and format data.
## =============================================================================

# List files to read in (.csv with file names including "data" and "09-03" 
#  [data read in on 3 Sept 2024]).
filenames <- list.files(pattern = "09_03")[grepl("data", list.files(pattern = "09_03"))]
# Read in all TOMST .csv files listed in 'filenames' object.
ldf <- lapply(filenames, function (x) { read.csv(x, header = FALSE, sep = ";") } )

# Name columns and apply to each list element
cols <- c('index', 'date.time', 'time.zone', 'soil.temp', 'surface.temp', 'air.temp',
          'soil.moisture', 'shake', 'errFlag')
ldf <- lapply(ldf, setNames, cols)

# Name each list element with its corresponding .csv file name
names(ldf) <- gsub(filenames, pattern="\\..*", replacement="")

# Combine all list elements into one dataframe
temp <- do.call(rbind, ldf)
# Assign rownames (previously names of list elements) to 'sensor.id' column
temp$sensor.id <- substr(rownames(temp), start = 6, stop = 13)
head(temp) #looks good!
rownames(temp) <- NULL #remove rownames

## Reformatting date and time information
temp$r_date <- strptime(temp$date.time, '%Y.%m.%d %H:%M') #convert to POSIX
temp$r_date_only <- format(temp$r_date, '%F') #date (YY-MM-DD) without time
temp$r_year <- as.numeric(format(temp$r_date, '%Y')) #year
temp$r_month <- as.numeric(format(temp$r_date, '%m')) #month
temp$r_day <- as.numeric(format(temp$r_date, '%d')) #day
temp$r_doy <- as.numeric(format(temp$r_date, '%j')) #DOY
temp$r_hour <- as.numeric(format(temp$r_date, '%H')) #hour
# if hour is greater than 8 and less than 20 --> day, else night
# (not 100% sure about whether this is correct -- check time zone, sensors might
#  be on GMT rather than local time zone)
temp$day.night <- with(temp, ifelse(r_hour > 8 & r_hour < 20, 'day', 'night'))

temp$r_date_ct <- as.POSIXct(temp$r_date) #convert to POSIXct

# Subset to observations starting after 125 DOY of 2023
temp2 <- temp[which(temp$r_year == 2024 |
                      temp$r_year == 2023 & temp$r_doy > 125),]


## =============================================================================
## (2) Summarize temperature and moisture data.
## =============================================================================

# Summarize (mean, standard deviation, number of observations) of soil 
# temperature, soil surface temperature, air temperature, and
#  soil moisture by day/night for each day
temp_summary <- doBy::summaryBy(soil.temp + surface.temp + air.temp + soil.moisture ~
                                  #sensor.id + #can add or omit summarizing by sensor ID 
                                  r_date_only + day.night,
                                id =~ r_doy + r_year + r_month + r_day,
                                data = temp2, FUN = c(mean, sd, length))

# Summarize (mean, standard deviation, number of observations) of soil 
#  temperature, soil surface temperature, air temperature, and
#  soil moisture by day/night for each day
temp_min_max <- doBy::summaryBy(soil.temp + surface.temp + air.temp + soil.moisture ~
                                  #sensor.id + #can add or omit summarizing by sensor ID
                                  r_date_only,
                                id =~ r_doy + r_year + r_month + r_day,
                                data = temp2, FUN = c(min, max, median, sd, length))

# Summarize (mean, median, and quantiles) of soil temperature, soil surface 
#  temperature, air temperature, and soil moisture for each day
temp_quan <- doBy::summaryBy(soil.temp + surface.temp + air.temp + soil.moisture ~
                                  #sensor.id + #can add or omit summarizing by sensor ID 
                                  r_date_only,
                                id =~ r_doy + r_year + r_month + r_day,
                                data = temp2, FUN = c(mean, median,
                                                      function(x) quantile(x, probs = 0.25), 
                                                      function(x) quantile(x, probs = 0.75),
                                                      sd, length),
                             fun.names = c('mean', 'median', 'quan25', 'quan75',
                                           'sd', 'length'))


## =============================================================================
## (3) Initial visualizations.
## =============================================================================

# Note: These vary in whether they are plotting air temperature vs soil temperature
# as the response. Play around as needed!

## Using full data set:

# Plot soil temperature through time with GAM smoothed fit line
ggplot(data = subset(temp2, r_year == 2024 & r_month > 3), 
       aes(x = as.POSIXct(r_date_only), y = soil.temp)) +
  geom_line() +
  stat_smooth() +
  theme_bw() +
  labs(x = "", y = "Soil temperature")

# Plot soil temperature through time (separated by day or night) with GAM smoothed
#  fit line
ggplot(data = temp2, aes(x = as.POSIXct(r_date_only), y = air.temp,
                         color = day.night)) +
  geom_line() +
  stat_smooth()+
  theme_bw() +
  labs(x = "", y = "Soil temperature")

## Using summary data sets:

# Plot mean soil surface temperature per day, through time
ggplot(data = temp_summary,
       aes(x = as.POSIXct(r_date_only), y = surface.temp.mean,
           color = day.night)) +
  geom_line() +
  stat_smooth()+
  theme_bw() +
  labs(x = "", y = "Air temperature")

# Plot max and min air temperature through time
ggplot(data = subset(temp_min_max, r_year == 2024 & r_month > 3)) +
  geom_line(aes(x = as.POSIXct(r_date_only), y = air.temp.min),
            color = "blue") +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = air.temp.min),
              color = "blue") +
  geom_line(aes(x = as.POSIXct(r_date_only), y = air.temp.max),
            color = "red") +
stat_smooth(aes(x = as.POSIXct(r_date_only), y = air.temp.max),
            color = "red") +
  geom_line(aes(x = as.POSIXct(r_date_only), y = air.temp.median)) +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = air.temp.median), color = "black")+
  theme_bw() +
  labs(x = "", y = "Air temperature")

# Plot quantiles of air temperature through time
ggplot(data = temp_quan) +
  geom_line(aes(x = as.POSIXct(r_date_only), y = surface.temp.quan25),
            color = "blue") +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = surface.temp.quan25),
              color = "blue") +
  geom_line(aes(x = as.POSIXct(r_date_only), y = surface.temp.quan75),
            color = "red") +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = surface.temp.quan75),
              color = "red") +
  theme_bw() +
  labs(x = "", y = "Air temperature")

# Plot max and min soil temperature through time
ggplot(data = subset(temp_summary, r_year == 2024 & r_month > 3)) +
  geom_line(aes(x = as.POSIXct(r_date_only), y = soil.temp.mean, 
                color = day.night)) +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = soil.temp.mean, 
                  color = day.night)) +
  theme_bw() +
  labs(x = "", y = "Soil temperature")

# Plot max and min soil temperature through time
ggplot(data = subset(temp_min_max, r_year == 2024 & r_month > 3)) +
  geom_line(aes(x = as.POSIXct(r_date_only), y = soil.temp.min),
            color = "blue") +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = soil.temp.min),
              color = "blue") +
  geom_line(aes(x = as.POSIXct(r_date_only), y = soil.temp.max),
            color = "red") +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = soil.temp.max),
              color = "red") +
  geom_line(aes(x = as.POSIXct(r_date_only), y = soil.temp.median)) +
  stat_smooth(aes(x = as.POSIXct(r_date_only), y = soil.temp.median), color = "black")+
  theme_bw() +
  labs(x = "", y = "Soil temperature")


## Looking at variation among sensors (field plots).

# During growing season months
ggplot(data = temp2[which((temp2$r_month > 3 & temp2$r_month < 8) & 
                            temp2$r_year == 2024),]) +
  stat_smooth(aes(x = r_hour, y = soil.temp, 
                  group = interaction(r_doy, as.factor(sensor.id)), 
                          color = as.factor(sensor.id))) +
  facet_wrap(~r_month)

# Daily for the month of May
ggplot(data = temp2[which((temp2$r_month  == 5) & 
                            temp2$r_year == 2024),]) +
  stat_smooth(aes(x = r_hour, y = soil.temp, 
                  group = interaction(r_doy, as.factor(sensor.id)), 
                  color = as.factor(sensor.id))) +
  facet_wrap(~r_doy)

## =============================================================================
## (3) Crude daily soil temp min/max values for two-week periods starting Apr 1,
##     2024
## =============================================================================

# April 1 - April 15, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-04-01" & r_date_only <= "2024-04-15") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))
## mean_min_soil_temp mean_max_soil_temp mean_median_soil_temp
##                7.7            11.1125              8.929167

# April 16 - April 30, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-04-16" & r_date_only <= "2024-04-30") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))
## mean_min_soil_temp mean_max_soil_temp mean_median_soil_temp
##              9.475           12.86667              10.98333

# May 1 - May 15, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-05-01" & r_date_only <= "2024-05-15") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))
## mean_min_soil_temp mean_max_soil_temp mean_median_soil_temp
##           12.46667           14.58333              13.27083


# May 16 - May 31, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-05-16" & r_date_only <= "2024-05-31") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# June 1 - June 15, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-06-01" & r_date_only <= "2024-06-15") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# June 16 - June 30, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-06-15" & r_date_only <= "2024-06-30") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# July 1 - July 15, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-07-01" & r_date_only <= "2024-07-15") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# July 16 - July 31, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-07-16" & r_date_only <= "2024-07-31") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# August 1 - August 15, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-08-01" & r_date_only <= "2024-08-15") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))

# August 16 - August 31, 2024
temp_min_max %>%
  filter(r_date_only >= "2024-08-16" & r_date_only <= "2024-08-31") %>%
  summarize(mean_min_soil_temp = mean(soil.temp.min, na.rm = TRUE),
            mean_max_soil_temp = mean(soil.temp.max, na.rm = TRUE),
            mean_median_soil_temp = mean(soil.temp.median, na.rm = TRUE))





