# Libraries
library(tidyverse)
library(lubridate)

# Read emergence file
emergence <- read.csv("../data_sheets/TT25_plant_emergence.csv") %>%
  mutate(#emerged = mdy(emerged),
         is_emerged = ifelse(is.na(emerged), "no", "yes"),
         #unwhorled = mdy(unwhorled),
         is_unwhorled = ifelse(is.na(unwhorled), "no", "yes"),
         #senesced = mdy(senesced),
         is_senesced = ifelse(is.na(senesced), "no", "yes")) %>%
  data.frame()
head(emergence)


is.na(emergence$emerged)

# Some summary statistics
emerge_summary <- emergence %>%
  group_by(spp, FullExpTrt) %>%
    
    
    
    
    
    
    
    n = length(color),
    n_emerged = n(is_emerged),
    percent_emerged = n_emerged / n * 100)
    
    ,
    
    n = length(color),
    n_unwhorled = length(!is.na(unwhorled)),
    percent_unwhorled = n_unwhorled / n,
    
    n = length(color),
    n_senesced = length(!is.na(senesced)),
    percent_senesced = n_senesced / n
    
    
  )
