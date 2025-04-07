# TT25_fit_response_curves.R
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
library(plantecophys)

# Load custom functions for cleaning LI-6800 files,
# standardizing Vcmax/Jmax/Rd to single temperature,
# and calculating 
R.utils::sourceDirectory("../functions/")

#####################################################################
# Clean raw licor files and put in cleaned subfolder
#####################################################################
# clean_licor_files(directory_path = "../data/li6800_aci/raw/",
#                   write_directory = "../data/li6800_aci/clean/")

#####################################################################
# Merge cleaned LI-6800 files into single file
#####################################################################
# List and set file names within cleaned_li6800 subfolder
files <- list.files(path = "../data/li6800_aci/clean/",
                    recursive = T,
                    pattern = "\\.csv$",
                    full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))

# Read all files and merge into central data frame
li6800_merged <- plyr::rbind.fill(lapply(files, read.csv)) %>%
  mutate(date = lubridate::ymd_hms(date),
         date_only = stringr::word(date, 1)) %>%
  mutate(doy = yday(date_only),
         Qin_cuvette = 300) %>%
  dplyr::select(obs, time, elapsed, date, date_only, doy, hhmmss, id, machine,
                A:Ci, gsw, Tleaf, VPDleaf, CO2_r, Asty, Adyn, Qin_cuvette, Flow_r) %>%
  arrange(machine, date, obs)

#####################################################################
# 04-01-2025 ALBERT
#####################################################################

#########################################
# 5015_mai_orange (leaf width = 1.9 cm)
#########################################

# Leaf area correction
mai_5015_orange_LA <- pi * (1.9/2)^2
mai_5015_orange_multiplier <- 6 / mai_5015_orange_LA
mai_5015_orange <- subset(li6800_merged, id == "5015_mai_orange") %>%
  mutate(A_corrected = Asty * mai_5015_orange_multiplier,
         gsw_corrected = gsw * mai_5015_orange_multiplier)

# Snapshot measurement
mai_5015_orange_snap <- mai_5015_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 5015,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5015_orange_fit <- mai_5015_orange %>% fitaci(
  varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                  Ci = "Ci", PPFD = "Qin_cuvette"),
  fitTPU = TRUE, Tcorrect = FALSE, citransition = 200)
plot(mai_5015_orange_fit)
summary(mai_5015_orange_fit)

# Write to data frame
aci_coefs <- data.frame(mai_5015_orange_snap, t(coef(mai_5015_orange_fit)))

#########################################
# 1941_mai_orange 
#########################################
mai_1941_orange <- subset(li6800_merged, id == "1941_mai_orange")

# Snapshot measurement
mai_1941_orange_snap <- mai_1941_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 1941,
         spp = "Mai",
         trt = "orange",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_1941_orange_fit <-mai_1941_orange %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1941_orange_fit)
summary(mai_1941_orange_fit)

# Write to data frame
aci_coefs[2,] <- c(mai_1941_orange_snap, t(coef(mai_1941_orange_fit)))

#########################################
# 3993_mai_orange (Leaf width = 2.2 cm) 
#########################################

# Leaf area correction
mai_3993_orange_LA <- pi * (2.2/2)^2
mai_3993_orange_multiplier <- 6 / mai_3993_orange_LA 
mai_3993_orange <- subset(li6800_merged, id == "3993_mai_orange") %>%
  mutate(A_corrected = Asty * mai_3993_orange_multiplier,
         gsw_corrected = gsw * mai_3993_orange_multiplier)

# Snapshot measurement
mai_3993_orange_snap <- mai_3993_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 3993,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_3993_orange_fit <- mai_3993_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3993_orange_fit)
summary(mai_3993_orange_fit)

# Write to data frame
aci_coefs[3,] <- c(mai_3993_orange_snap, t(coef(mai_3993_orange_fit)))

#####################################################################
# 04-01-2025 GIBSON
#####################################################################

#########################################
# 5009_mai_pink (leaf width = 1.6 cm)
#########################################

# Leaf area correction
mai_5009_pink_LA <- pi * (1.6/2)^2
mai_5009_pink_multiplier <- 6 / mai_5009_pink_LA 
mai_5009_pink <- subset(li6800_merged, id == "5009_mai_pink") %>%
  mutate(A_corrected = Asty * mai_5009_pink_multiplier,
         gsw_corrected = gsw * mai_5009_pink_multiplier)

# Snapshot measurement
mai_5009_pink_snap <- mai_5009_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = 5009,
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5009_pink_fit <- mai_5009_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5009_pink_fit)
summary(mai_5009_pink_fit)
## Note: really bad A/Ci curve, not including in dataset

# Write to data frame
aci_coefs[4,] <- c(mai_5009_pink_snap, NA, NA, NA, NA)

#########################################
# TT24_135_mai_pink
#########################################
mai_TT24_135_pink <- subset(li6800_merged, id == "TT24_135_mai_pink") 

# Snapshot measurement
mai_TT24_135_pink_snap <- mai_TT24_135_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_135",
         spp = "Mai",
         trt = "pink",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_135_pink_fit <- mai_TT24_135_pink %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_135_pink_fit)
summary(mai_TT24_135_pink_fit)

# Write to data frame
aci_coefs[5,] <- c(mai_TT24_135_pink_snap, t(coef(mai_TT24_135_pink_fit)))

#########################################
# flag2_mai_pink
#########################################
mai_flag2_pink <- subset(li6800_merged, id == "flag2_mai_pink") 

# Snapshot measurement
mai_flag2_pink_snap <- mai_flag2_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = "flag2",
         spp = "Mai",
         trt = "pink",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_flag2_pink_fit <- mai_flag2_pink %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_flag2_pink_fit)
summary(mai_flag2_pink_fit)

# Write to data frame
aci_coefs[6,] <- c(mai_flag2_pink_snap, t(coef(mai_flag2_pink_fit)))

#########################################
# 2608_mai_pink (Leaf width = 1.8 cm)
#########################################

# Leaf area correction
mai_2608_pink_LA <- pi * (1.8/2)^2
mai_2608_pink_multiplier <- 6 / mai_2608_pink_LA 
mai_2608_pink <- subset(li6800_merged, id == "2608_mai_pink") %>%
  mutate(A_corrected = Asty * mai_2608_pink_multiplier,
         gsw_corrected = gsw * mai_2608_pink_multiplier)

# Snapshot measurement
mai_2608_pink_snap <- mai_2608_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = 2608,
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2608_pink_fit <- mai_2608_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2608_pink_fit)
summary(mai_2608_pink_fit)
## Note: really bad A/Ci curve, not including in dataset

# Write to data frame
aci_coefs[7,] <- c(mai_2608_pink_snap, NA, NA, NA, NA)

#####################################################################
# 04-01-2025 OZZIE
#####################################################################

#########################################
# TT24_mai_green (Leaf width = 2.0 cm)
#########################################
# Note: TT24_230 but missing number

# Leaf area correction
mai_TT24_230_green_LA <- pi * (2.0/2)^2
mai_TT24_230_green_multiplier <- 6 / mai_TT24_230_green_LA 
mai_TT24_230_green <- subset(li6800_merged, id == "TT24_mai_green") %>%
  mutate(A_corrected = Asty * mai_TT24_230_green_multiplier,
         gsw_corrected = gsw * mai_TT24_230_green_multiplier)

# Snapshot measurement
mai_TT24_230_green_snap <- mai_TT24_230_green %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_230",
         spp = "Mai",
         trt = "green",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_230_green_fit <- mai_TT24_230_green %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 600)
plot(mai_TT24_230_green_fit)
summary(mai_TT24_230_green_fit)

# Write to data frame
aci_coefs[8,] <- c(mai_TT24_230_green_snap, t(coef(mai_TT24_230_green_fit)))

#########################################
# 4781_mai_green (note duplicate - using 
# curve on 4-3-2025 as curve looks cleaner)
#########################################

#########################################
# 2285_mai_green (note duplicate - using 
# curve on 4-1-2025 because Ci extends 
# beyond curve on 4-3-2025 and fluxes 
# are very similar between measurement 
# days at lower CO2 values)
#########################################
mai_2285_green <- subset(li6800_merged, 
                         id == "2285_mai_green" & 
                           date_only == "2025-04-01") 

# Snapshot measurement
mai_2285_green_snap <- mai_2285_green %>%
  filter(row_number() == 1) %>%
  mutate(id = 2285,
         spp = "Mai",
         trt = "green",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2285_green_fit <- mai_2285_green %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2285_green_fit)
summary(mai_2285_green_fit)

# Write to data frame
aci_coefs[9,] <- c(mai_2285_green_snap, t(coef(mai_2285_green_fit)))

#####################################################################
# 04-01-2025 STAN
#####################################################################

#########################################
# 5567_mai_blue
#########################################
mai_5567_blue <- subset(li6800_merged, id == "5567_mai_blue") 

# Snapshot measurement
mai_5567_blue_snap <- mai_5567_blue %>%
  filter(row_number() == 1) %>%
  mutate(id = 5567,
         spp = "Mai",
         trt = "blue",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5567_blue_fit <- mai_5567_blue %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE,  citransition = 500)
plot(mai_5567_blue_fit)
summary(mai_5567_blue_fit)

# Write to data frame
aci_coefs[10,] <- c(mai_5567_blue_snap, t(coef(mai_5567_blue_fit)))

#########################################
# 4216_mai_blue
#########################################
mai_4216_blue <- subset(li6800_merged, id == "4216_mai_blue") 

# Snapshot measurement
mai_4216_blue_snap <- mai_4216_blue %>%
  filter(row_number() == 1) %>%
  mutate(id = 4216,
         spp = "Mai",
         trt = "blue",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4216_blue_fit <- mai_4216_blue %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE,  citransition = 500)
plot(mai_4216_blue_fit)
summary(mai_4216_blue_fit)

# Write to data frame
aci_coefs[11,] <- c(mai_4216_blue_snap, t(coef(mai_4216_blue_fit)))

#########################################
# 4883_mai_blue
#########################################
mai_4883_blue <- subset(li6800_merged, id == "4883_mai_blue") 

# Snapshot measurement
mai_4883_blue_snap <- mai_4883_blue %>%
  filter(row_number() == 1) %>%
  mutate(id = 4883,
         spp = "Mai",
         trt = "blue",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4883_blue_fit <- mai_4883_blue %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4883_blue_fit)
summary(mai_4883_blue_fit)

# Write to data frame
aci_coefs[12,] <- c(mai_4883_blue_snap, t(coef(mai_4883_blue_fit)))

#########################################
# 4883_mai_blue
#########################################
mai_4883_blue <- subset(li6800_merged, id == "4883_mai_blue") 

# Snapshot measurement
mai_4883_blue_snap <- mai_4883_blue %>%
  filter(row_number() == 1) %>%
  mutate(id = 4883,
         spp = "Mai",
         trt = "blue",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4883_blue_fit <- mai_4883_blue %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4883_blue_fit)
summary(mai_4883_blue_fit)

# Write to data frame
aci_coefs[12,] <- c(mai_4883_blue_snap, t(coef(mai_4883_blue_fit)))

#########################################
# 4859_mai_blue
#########################################
mai_4859_blue <- subset(li6800_merged, id == "4859_mai_blue") 

# Snapshot measurement
mai_4859_blue_snap <- mai_4859_blue %>%
  filter(row_number() == 1) %>%
  mutate(id = 4859,
         spp = "Mai",
         trt = "blue",
         anet = Asty,
         gsw = gsw,
         iwue = Asty / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4859_blue_fit <- mai_4859_blue %>%
  fitaci(varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4859_blue_fit)
summary(mai_4859_blue_fit)
## Bad curve fit, omitting from dataset

# Write to data frame
aci_coefs[13,] <- c(mai_4859_blue_snap, NA, NA, NA, NA)

#####################################################################
# 04-02-2025 OZZIE
#####################################################################

#########################################
# 54_mai_white
#########################################
mai_54_white <- subset(li6800_merged, id == "54_mai_white") 

# Snapshot measurement
mai_54_white_snap <- mai_54_white %>%
  filter(row_number() == 1) %>%
  mutate(id = 54,
         spp = "Mai",
         trt = "white",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_54_white_fit <- mai_54_white %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_54_white_fit)
summary(mai_54_white_fit)

# Write to data frame
aci_coefs[14,] <- c(mai_54_white_snap, t(coef(mai_54_white_fit)))

#########################################
# 1468_mai_white
#########################################
mai_1468_white <- subset(li6800_merged, id == "1468_mai_white") 

# Snapshot measurement
mai_1468_white_snap <- mai_1468_white %>%
  filter(row_number() == 1) %>%
  mutate(id = 1468,
         spp = "Mai",
         trt = "white",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_1468_white_fit <- mai_1468_white %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1468_white_fit)
summary(mai_1468_white_fit)

# Write to data frame
aci_coefs[15,] <- c(mai_1468_white_snap, t(coef(mai_1468_white_fit)))

#########################################
# 3438_mai_white
#########################################
mai_3438_white <- subset(li6800_merged, id == "3438_mai_white") 

# Snapshot measurement
mai_3438_white_snap <- mai_3438_white %>%
  filter(row_number() == 1) %>%
  mutate(id = 3438,
         spp = "Mai",
         trt = "white",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_3438_white_fit <- mai_3438_white %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3438_white_fit)
summary(mai_3438_white_fit)

# Write to data frame
aci_coefs[16,] <- c(mai_3438_white_snap, t(coef(mai_3438_white_fit)))

#########################################
# 2265_mai_white
#########################################
mai_2265_white <- subset(li6800_merged, id == "2265_mai_white") 

# Snapshot measurement
mai_2265_white_snap <- mai_2265_white %>%
  filter(row_number() == 1) %>%
  mutate(id = 2265,
         spp = "Mai",
         trt = "white",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2265_white_fit <- mai_2265_white %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2265_white_fit)
summary(mai_2265_white_fit)

# Write to data frame
aci_coefs[17,] <- c(mai_2265_white_snap, t(coef(mai_2265_white_fit)))

#########################################
# 2265_mai_red
#########################################
mai_1264_red <- subset(li6800_merged, id == "1264_mai_red") 

# Snapshot measurement
mai_1264_red_snap <- mai_1264_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 1264,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_1264_red_fit <- mai_1264_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1264_red_fit)
summary(mai_1264_red_fit)

# Write to data frame
aci_coefs[18,] <- c(mai_1264_red_snap, t(coef(mai_1264_red_fit)))

#########################################
# 4740_mai_red
#########################################
mai_4740_red <- subset(li6800_merged, id == "4740_mai_red") 

# Snapshot measurement
mai_4740_red_snap <- mai_4740_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 4740,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4740_red_fit <- mai_4740_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4740_red_fit)
summary(mai_4740_red_fit)

# Write to data frame
aci_coefs[19,] <- c(mai_4740_red_snap, t(coef(mai_4740_red_fit)))

#####################################################################
# 04-02-2025 STAN
#####################################################################

#########################################
# 5074_mai_red
#########################################
mai_5074_red <- subset(li6800_merged, id == "5074_mai_red") 

# Snapshot measurement
mai_5074_red_snap <- mai_5074_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 5074,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5074_red_fit <- mai_5074_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 200)
plot(mai_5074_red_fit)
summary(mai_5074_red_fit)

# Write to data frame
aci_coefs[20,] <- c(mai_5074_red_snap, t(coef(mai_5074_red_fit)))

#########################################
# 5596_mai_red
#########################################
mai_5596_red <- subset(li6800_merged, id == "5596_mai_red") 

# Snapshot measurement
mai_5596_red_snap <- mai_5596_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 5596,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5596_red_fit <- mai_5596_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 400)
plot(mai_5596_red_fit)
summary(mai_5596_red_fit)

# Write to data frame
aci_coefs[21,] <- c(mai_5596_red_snap, t(coef(mai_5596_red_fit)))

#########################################
# 2639_mai_red
#########################################
mai_2639_red <- subset(li6800_merged, id == "2639_mai_red") 

# Snapshot measurement
mai_2639_red_snap <- mai_2639_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 2639,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2639_red_fit <- mai_2639_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2639_red_fit)
summary(mai_2639_red_fit)

# Write to data frame
aci_coefs[22,] <- c(mai_2639_red_snap, t(coef(mai_2639_red_fit)))

#########################################
# TT24_106_mai_red
#########################################
mai_TT24_106_red <- subset(li6800_merged, id == "TT24_106_mai_red") 

# Snapshot measurement
mai_TT24_106_red_snap <- mai_TT24_106_red %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_106",
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_106_red_fit <- mai_TT24_106_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_106_red_fit)
summary(mai_TT24_106_red_fit)

# Write to data frame
aci_coefs[23,] <- c(mai_TT24_106_red_snap, t(coef(mai_TT24_106_red_fit)))

#########################################
# TT24_111_mai_red
#########################################
mai_TT24_111_red <- subset(li6800_merged, id == "TT24_111_mai_red") 

# Snapshot measurement
mai_TT24_111_red_snap <- mai_TT24_111_red %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_111",
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_111_red_fit <- mai_TT24_111_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_111_red_fit)
summary(mai_TT24_111_red_fit)

# Write to data frame
aci_coefs[24,] <- c(mai_TT24_111_red_snap, t(coef(mai_TT24_111_red_fit)))

#########################################
# 4844_mai_red
#########################################
mai_4844_red <- subset(li6800_merged, id == "4844_mai_red" & A > 0) 

# Snapshot measurement
mai_4844_red_snap <- mai_4844_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 4844,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4844_red_fit <- mai_4844_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4844_red_fit)
summary(mai_4844_red_fit)

# Write to data frame
aci_coefs[25,] <- c(mai_4844_red_snap, t(coef(mai_4844_red_fit)))

#########################################
# 5600_mai_red
#########################################
mai_5600_red <- subset(li6800_merged, id == "5600_mai_red") 

# Snapshot measurement
mai_5600_red_snap <- mai_5600_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 5600,
         spp = "Mai",
         trt = "red",
         anet = A,
         gsw = gsw,
         iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5600_red_fit <- mai_5600_red %>%
  fitaci(varnames = list(ALEAF = "A", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5600_red_fit)
summary(mai_5600_red_fit)

# Write to data frame
aci_coefs[26,] <- c(mai_5600_red_snap, t(coef(mai_5600_red_fit)))

#####################################################################
# 04-02-2025 GIBSON
#####################################################################

#########################################
# 4193_mai_pink
#########################################
# Leaf area correction
mai_4193_pink_LA <- pi * (1.9/2)^2
mai_4193_pink_multiplier <- 6 / mai_4193_pink_LA 
mai_4193_pink <- subset(li6800_merged, id == "4193_mai_pink") %>%
  mutate(A_corrected = A * mai_4193_pink_multiplier,
         gsw_corrected = gsw * mai_4193_pink_multiplier)

# Snapshot measurement
mai_4193_pink_snap <- mai_4193_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = 4193,
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4193_pink_fit <- mai_4193_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4193_pink_fit)
summary(mai_4193_pink_fit)

# Write to data frame
aci_coefs[27,] <- c(mai_4193_pink_snap, t(coef(mai_4193_pink_fit)))

#########################################
# 4061_mai_pink
#########################################
# Leaf area correction
mai_4061_pink_LA <- pi * (1.5/2)^2
mai_4061_pink_multiplier <- 6 / mai_4061_pink_LA 
mai_4061_pink <- subset(li6800_merged, id == "4061_mai_pink") %>%
  mutate(A_corrected = A * mai_4061_pink_multiplier,
         gsw_corrected = gsw * mai_4061_pink_multiplier)

# Snapshot measurement
mai_4061_pink_snap <- mai_4061_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = 4061,
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4061_pink_fit <- mai_4061_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_4061_pink_fit)
summary(mai_4061_pink_fit)

# Write to data frame
aci_coefs[28,] <- c(mai_4061_pink_snap, t(coef(mai_4061_pink_fit)))

#########################################
# TT24_104_mai_pink
#########################################
# Leaf area correction
mai_TT24_104_pink_LA <- pi * (2.7/2)^2
mai_TT24_104_pink_multiplier <- 6 / mai_TT24_104_pink_LA 
mai_TT24_104_pink <- subset(li6800_merged, id == "TT24_104_mai_pink") %>%
  mutate(A_corrected = A * mai_TT24_104_pink_multiplier,
         gsw_corrected = gsw * mai_TT24_104_pink_multiplier)

# Snapshot measurement
mai_TT24_104_pink_snap <- mai_TT24_104_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_104",
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_104_pink_fit <- mai_TT24_104_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_104_pink_fit)
summary(mai_TT24_104_pink_fit)

# Write to data frame
aci_coefs[29,] <- c(mai_TT24_104_pink_snap, t(coef(mai_TT24_104_pink_fit)))

#########################################
# 4237_mai_pink
#########################################
# Leaf area correction
mai_4237_pink_LA <- pi * (1.5/2)^2
mai_4237_pink_multiplier <- 6 / mai_4237_pink_LA 
mai_4237_pink <- subset(li6800_merged, id == "4237_mai_pink") %>%
  mutate(A_corrected = A * mai_4237_pink_multiplier,
         gsw_corrected = gsw * mai_4237_pink_multiplier)

# Snapshot measurement
mai_4237_pink_snap <- mai_4237_pink %>%
  filter(row_number() == 1) %>%
  mutate(id = 4237,
         spp = "Mai",
         trt = "pink",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_4237_pink_fit <- mai_4237_pink %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 200)
plot(mai_4237_pink_fit)
summary(mai_4237_pink_fit)

# Write to data frame
aci_coefs[30,] <- c(mai_4237_pink_snap, t(coef(mai_4237_pink_fit)))

#########################################
# 5178_mai_red
#########################################
# Leaf area correction
mai_5178_red_LA <- pi * (2.0/2)^2
mai_5178_red_multiplier <- 6 / mai_5178_red_LA 
mai_5178_red <- subset(li6800_merged, id == "5178_mai_red") %>%
  mutate(A_corrected = A * mai_5178_red_multiplier,
         gsw_corrected = gsw * mai_5178_red_multiplier)

# Snapshot measurement
mai_5178_red_snap <- mai_5178_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 5178,
         spp = "Mai",
         trt = "red",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5178_red_fit <- mai_5178_red %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 200)
plot(mai_5178_red_fit)
summary(mai_5178_red_fit)

# Write to data frame
aci_coefs[31,] <- c(mai_5178_red_snap, t(coef(mai_5178_red_fit)))

#########################################
# 5122_mai_red
#########################################
# Leaf area correction
mai_5122_red_LA <- pi * (2.1/2)^2
mai_5122_red_multiplier <- 6 / mai_5122_red_LA 
mai_5122_red <- subset(li6800_merged, id == "5122_mai_red") %>%
  mutate(A_corrected = A * mai_5122_red_multiplier,
         gsw_corrected = gsw * mai_5122_red_multiplier)

# Snapshot measurement
mai_5122_red_snap <- mai_5122_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 5122,
         spp = "Mai",
         trt = "red",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_5122_red_fit <- mai_5122_red %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_5122_red_fit)
summary(mai_5122_red_fit)

# Write to data frame
aci_coefs[32,] <- c(mai_5122_red_snap, t(coef(mai_5122_red_fit)))

#########################################
# 545_mai_red
#########################################
# Leaf area correction
mai_545_red_LA <- pi * (1.9/2)^2
mai_545_red_multiplier <- 6 / mai_545_red_LA 
mai_545_red <- subset(li6800_merged, id == "545_mai_red") %>%
  mutate(A_corrected = A * mai_545_red_multiplier,
         gsw_corrected = gsw * mai_545_red_multiplier)

# Snapshot measurement
mai_545_red_snap <- mai_545_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 545,
         spp = "Mai",
         trt = "red",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_545_red_fit <- mai_545_red %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 400)
plot(mai_545_red_fit)
summary(mai_545_red_fit)

# Write to data frame
aci_coefs[33,] <- c(mai_545_red_snap, t(coef(mai_545_red_fit)))

#####################################################################
# 04-02-2025 ALBERT
#####################################################################

#########################################
# 2799_mai_orange
#########################################
# Leaf area correction
mai_2799_orange_LA <- pi * (2.4/2)^2
mai_2799_orange_multiplier <- 6 / mai_2799_orange_LA 
mai_2799_orange <- subset(li6800_merged, id == "2799_mai_orange" & A > 0.1) %>%
  mutate(A_corrected = A * mai_2799_orange_multiplier,
         gsw_corrected = gsw * mai_2799_orange_multiplier)

# Snapshot measurement
mai_2799_orange_snap <- mai_2799_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 2799,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2799_orange_fit <- mai_2799_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE, citransition = 600)
plot(mai_2799_orange_fit)
summary(mai_2799_orange_fit)

# Write to data frame
aci_coefs[34,] <- c(mai_2799_orange_snap, t(coef(mai_2799_orange_fit)))

#########################################
# 3136_mai_orange
#########################################
# Leaf area correction
mai_3136_orange_LA <- pi * (1.8/2)^2
mai_3136_orange_multiplier <- 6 / mai_3136_orange_LA 
mai_3136_orange <- subset(li6800_merged, id == "3136_mai_orange") %>%
  mutate(A_corrected = A * mai_3136_orange_multiplier,
         gsw_corrected = gsw * mai_3136_orange_multiplier)

# Snapshot measurement
mai_3136_orange_snap <- mai_3136_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 3136,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_3136_orange_fit <- mai_3136_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_3136_orange_fit)
summary(mai_3136_orange_fit)

# Write to data frame
aci_coefs[35,] <- c(mai_3136_orange_snap, t(coef(mai_3136_orange_fit)))

#########################################
# TT24_130_mai_orange
#########################################
# Leaf area correction
mai_TT24_130_orange_LA <- pi * (2.3/2)^2
mai_TT24_130_orange_multiplier <- 6 / mai_TT24_130_orange_LA 
mai_TT24_130_orange <- subset(li6800_merged, id == "TT24_130_mai_orange") %>%
  mutate(A_corrected = A * mai_TT24_130_orange_multiplier,
         gsw_corrected = gsw * mai_TT24_130_orange_multiplier)

# Snapshot measurement
mai_TT24_130_orange_snap <- mai_TT24_130_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = "TT24_130",
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_TT24_130_orange_fit <- mai_TT24_130_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_TT24_130_orange_fit)
summary(mai_TT24_130_orange_fit)

# Write to data frame
aci_coefs[36,] <- c(mai_TT24_130_orange_snap, t(coef(mai_TT24_130_orange_fit)))

#########################################
# 2546_mai_orange
#########################################
# Leaf area correction
mai_2546_orange_LA <- pi * (1.8/2)^2
mai_2546_orange_multiplier <- 6 / mai_2546_orange_LA 
mai_2546_orange <- subset(li6800_merged, id == "2546_mai_orange") %>%
  mutate(A_corrected = A * mai_2546_orange_multiplier,
         gsw_corrected = gsw * mai_2546_orange_multiplier)

# Snapshot measurement
mai_2546_orange_snap <- mai_2546_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 2546,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2546_orange_fit <- mai_2546_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2546_orange_fit)
summary(mai_2546_orange_fit)

# Write to data frame
aci_coefs[37,] <- c(mai_2546_orange_snap, t(coef(mai_2546_orange_fit)))

#########################################
# 2546_mai_orange
#########################################
# Leaf area correction
mai_2546_orange_LA <- pi * (1.8/2)^2
mai_2546_orange_multiplier <- 6 / mai_2546_orange_LA 
mai_2546_orange <- subset(li6800_merged, id == "2546_mai_orange") %>%
  mutate(A_corrected = A * mai_2546_orange_multiplier,
         gsw_corrected = gsw * mai_2546_orange_multiplier)

# Snapshot measurement
mai_2546_orange_snap <- mai_2546_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 2546,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2546_orange_fit <- mai_2546_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2546_orange_fit)
summary(mai_2546_orange_fit)

# Write to data frame
aci_coefs[37,] <- c(mai_2546_orange_snap, t(coef(mai_2546_orange_fit)))

#########################################
# 1689_mai_orange
#########################################
# Leaf area correction
mai_1689_orange_LA <- pi * (2.0/2)^2
mai_1689_orange_multiplier <- 6 / mai_1689_orange_LA 
mai_1689_orange <- subset(li6800_merged, id == "1689_mai_orange") %>%
  mutate(A_corrected = A * mai_1689_orange_multiplier,
         gsw_corrected = gsw * mai_1689_orange_multiplier)

# Snapshot measurement
mai_1689_orange_snap <- mai_1689_orange %>%
  filter(row_number() == 1) %>%
  mutate(id = 1689,
         spp = "Mai",
         trt = "orange",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_1689_orange_fit <- mai_1689_orange %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_1689_orange_fit)
summary(mai_1689_orange_fit)

# Write to data frame
aci_coefs[38,] <- c(mai_1689_orange_snap, t(coef(mai_1689_orange_fit)))

#########################################
# 2689_mai_red
#########################################
# Leaf area correction
mai_2689_red_LA <- pi * (1.9/2)^2
mai_2689_red_multiplier <- 6 / mai_2689_red_LA 
mai_2689_red <- subset(li6800_merged, id == "2689_mai_red") %>%
  mutate(A_corrected = A * mai_2689_red_multiplier,
         gsw_corrected = gsw * mai_2689_red_multiplier)

# Snapshot measurement
mai_2689_red_snap <- mai_2689_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 2689,
         spp = "Mai",
         trt = "red",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2689_red_fit <- mai_2689_red %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2689_red_fit)
summary(mai_2689_red_fit)

# Write to data frame
aci_coefs[39,] <- c(mai_2689_red_snap, t(coef(mai_2689_red_fit)))

#####################################################################
# 04-03-2025 OZZIE
#####################################################################

#########################################
# 2285_mai_green
#########################################
# Used curve from 04-01-2025

#########################################
# 4698_mai_green
#########################################

mai_2689_red <- subset(li6800_merged, id == "4698_mai_orange")

# Snapshot measurement
mai_2689_red_snap <- mai_2689_red %>%
  filter(row_number() == 1) %>%
  mutate(id = 2689,
         spp = "Mai",
         trt = "red",
         anet = A_corrected,
         gsw = gsw_corrected,
         iwue = A_corrected / gsw_corrected,
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, spp:ci.ca)

# A/Ci curve fit
mai_2689_red_fit <- mai_2689_red %>%
  fitaci(varnames = list(ALEAF = "A_corrected", Tleaf = "Tleaf", 
                         Ci = "Ci", PPFD = "Qin_cuvette"),
         fitTPU = TRUE, Tcorrect = FALSE)
plot(mai_2689_red_fit)
summary(mai_2689_red_fit)

# Write to data frame
aci_coefs[39,] <- c(mai_2689_red_snap, t(coef(mai_2689_red_fit)))

