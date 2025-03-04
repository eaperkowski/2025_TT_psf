# Libraries
library(tidyverse)
library(lubridate)

# Read emergence file
emergence <- read.csv("../data_sheets/TT25_plant_emergence.csv") %>%
  dplyr::select(-notes) %>%
  mutate(emerged = mdy(emerged),
         is_emerged = ifelse(is.na(emerged), "no", "yes"),
         unwhorled = mdy(unwhorled),
         is_unwhorled = ifelse(is.na(unwhorled), "no", "yes"),
         senesced = mdy(senesced),
         is_senesced = ifelse(is.na(senesced), "no", "yes")) %>%
  data.frame()
head(emergence)

# Some summary statistics (workaround bc for some reason dplyr
# is not allowing me to summarize counts by trt)
emerge_summary <- emergence %>%
  group_by(spp, FullExpTrt) %>%
  summarize(total_n = length(FullExpTrt))

emerged <- emergence %>%
  group_by(spp, FullExpTrt) %>%
  count(emerged_bool = is_emerged == "yes") %>%
  filter(emerged_bool == TRUE) %>%
  select(spp, FullExpTrt, emerged_n = n)

unwhorled <- emergence %>%
  group_by(spp, FullExpTrt) %>%
  count(unwhorled_bool = is_unwhorled == "yes") %>%
  filter(unwhorled_bool == TRUE) %>%
  select(spp, FullExpTrt, unwhorled_n = n)

emerge_summary_final <- emerge_summary %>%
  full_join(emerged) %>%
  full_join(unwhorled) %>%
  mutate(emerged_n = ifelse(is.na(emerged_n), 0, emerged_n),
         unwhorled_n = ifelse(is.na(unwhorled_n), 0, unwhorled_n),
         percent_emerged = emerged_n / total_n * 100,
         percent_unwhorled = unwhorled_n / total_n * 100) %>%
  separate(FullExpTrt, c("PlantGMTrt", "ExpSoilSource", "ExpFungSource"), 
           remove = FALSE) %>%
  mutate(PlantGMTrt = factor(str_replace(PlantGMTrt ,"Plant", ""),
                             levels = c("W", "NW")),
         ExpSoilSource = factor(str_replace(ExpSoilSource ,"Soil", ""),
                                levels = c("W", "NW")),
         ExpFungSource = factor(str_replace(ExpFungSource ,"Fung", ""),
                                levels = c("W", "NW", "Wsterile", "NWsterile")))

# Emergence plots
ggplot(data = subset(emerge_summary_final, spp == "Mai"),
       aes(x = PlantGMTrt, y = percent_emerged)) +
  geom_bar(stat = "identity") +
  labs(x = "Plant GM source", y = "% emerged") +
  facet_grid(ExpSoilSource ~ ExpFungSource) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggplot(data = subset(emerge_summary_final, spp == "Tri"),
       aes(x = PlantGMTrt, y = percent_emerged)) +
  geom_bar(stat = "identity") +
  labs(x = "Plant GM source", y = "% emerged") +
  scale_y_continuous(limits = c(0,100)) +
  facet_grid(ExpSoilSource ~ ExpFungSource) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Unwhorled plots
ggplot(data = subset(emerge_summary_final, spp == "Mai"),
       aes(x = PlantGMTrt, y = percent_unwhorled)) +
  geom_bar(stat = "identity") +
  labs(x = "Plant GM source", y = "% unwhorled") +
  facet_grid(ExpSoilSource ~ ExpFungSource) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = subset(emerge_summary_final, spp == "Tri"),
       aes(x = PlantGMTrt, y = percent_unwhorled)) +
  geom_bar(stat = "identity") +
  labs(x = "Plant GM source", y = "% unwhorled") +
  scale_y_continuous(limits = c(0,100)) +
  facet_grid(ExpSoilSource ~ ExpFungSource) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




