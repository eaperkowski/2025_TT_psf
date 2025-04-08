# TT25_create_phys_traits.R
# R script that merges LI-6800 files from "../cleaned_li6800/" and
# uses merged file to fit CO2 response curves
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries and custom functions
#####################################################################
# Libraries
library(tidyverse)

# Helper function for calculating Maianthemum total leaf area
# based on stem length
calc_leafarea_mai <- function(x) exp(log(x)*1.6) # where x = stem length in cm

#####################################################################
# Read data files and merge into single cleaned file
#####################################################################

# Read treatment summary
trt_summary <- read.csv("../data/TT25_treatment_summary.csv")

# Read photosynthesis and allometry datasets
photosynthesis <- read.csv("../data/TT25_gasExchange.csv") 
allometry <- read.csv("../data/TT25_stem_lengths.csv") %>%
  mutate(total_leaf_area = calc_leafarea_mai(stem_length_cm)) %>%
  group_by(id, trt) %>%
  summarize(n_stems = length(id),
            total_leaf_area = sum(total_leaf_area)) %>%
  mutate(spp = "Mai")

# Read fluorescence dataset, with some light cleaning for tidy merge
fluorescence <- read.csv("../data/li600/TT25_psf_2025_03_31_25.csv") %>%
  dplyr::select(Obs.:E_apparent, VPDleaf, Fs:PhiPS2, ETR) %>%
  separate(id, into = c("id", "spp", "trt"), sep = "_") %>%
  mutate(id = gsub(pattern = "\\.", replacement = "_", id),
         spp = "Mai") %>%
  dplyr::select(-Obs., -Time, -Date, -gsw)
head(fluorescence)

# Merge photosynthesis and allometry data frames
full_df <- photosynthesis %>%
  full_join(allometry) %>%
  full_join(fluorescence, by = c("id", "spp", "trt")) %>%
  full_join(trt_summary) %>%
  mutate(id = ifelse(id == "235", "TT24_235", id),
         full_trt = str_c(plantGMtrt, "_", expSoilSource, "_", ExpFungSource)) %>%
  dplyr::select(id, spp, full_trt, plantGMtrt:ExpFungSource, anet:ETR)

#####################################################################
# Save cleaned dataset
#####################################################################
write.csv(full_df, "../data/TT25_full_physiology.csv", row.names = F)
