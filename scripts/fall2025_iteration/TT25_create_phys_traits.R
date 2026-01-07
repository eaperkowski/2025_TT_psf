# Script cleans and fits A/Ci curve data from Trillium.
# Curves collected on Dec. 9 through Dec. 11. Conditions: 500 umol/s flow
# rate, 10 000 rpm mixing fan, 400 umol/mol CO2 starting point, 50% relative 
# humidity, 800 umol/m2/s light, 25degC leaf temperature

# Libraries
library(tidyverse)
library(plantecophys)
library(lme4)
library(car)
library(emmeans)

# Load function for cleaning raw LI-6800 files
source("../../functions/clean_licor_files.R")

# Load function for temperature standardizing (probably not needed)
source("../../functions/temp_standardize.R")

# Read treatment summary
trt_summary <- read.csv("../../data/2025_2026/TT25_tri_trt_list.csv")
head(trt_summary)

# Clean raw LI-6800 files
## clean_licor_files(directory_path = "../../data/2025_2026/li6800/raw/",
##                   write_directory = "../../data/2025_2026/li6800/clean/")

# Create dataframe that merges all cleaned LI-6800 files together
files <- list.files(path = "../../data/2025_2026/li6800/clean/",
                      recursive = T,
                      pattern = "\\.csv$",
                      full.names = T)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))

# Read all files and merge into central data frame
aci_merged_tri <- plyr::rbind.fill(lapply(files, read.csv)) %>%
  mutate(date = lubridate::ymd_hms(date),
         date_only = stringr::word(date, 1)) %>%
  mutate(doy = yday(date_only),
         Qin_cuvette = 800,
         keep_row = "yes") %>%
  dplyr::select(obs, time, elapsed, date, date_only, doy, hhmmss, id, machine,
                A:Ci, gsw, Tleaf, VPDleaf, CO2_r, Asty, Qin_cuvette, Flow_r, keep_row) %>%
  arrange(machine, date, obs)

# Get first row of each measurement for snapshop Anet, gsw measurements
aci_snapshot_tri <- aci_merged_tri %>%
  group_by(id, machine) %>%
  filter(row_number() == 1) %>%
  mutate(anet = Asty,
         iwue = ifelse((anet / gsw) < 0, NA, anet/gsw),
         ci.ca = Ci / Ca) %>%
  dplyr::select(id, machine, anet, gsw, iwue, ci.ca) %>%
  left_join(trt_summary)


##############################################################################
# Fit A/Ci curves (Trillium)
##############################################################################

# Remove measurements to improve curve fit
aci_merged_tri$keep_row[c(39, 52, 78, 92, 105, 118, 131, 144, 157, 170, 183, 196, 
                          222, 223, 224, 225, 228, 234, 235, 274, 300, 325, 326, 366, 368, 371, 376, 378, 391, 417, 
                          430, 432, 433, 443, 454, 457, 458, 469, 471, 472, 482, 
                          534, 547, 549, 560, 573, 586, 605, 612, 703, 716, 722, 729, 
                          730, 768, 784,
                          794, 820, 847, 898, 900, 904, 911, 924, 930, 937, 950, 963,
                          976, 996, 997, 998, 999, 1000, 1001, 1258,
                          1015, 1021, 1024, 1028, 1034, 1054, 1067, 1073, 1099, 
                          1275, 1289, 1290, 1291, 1296, 1317, 1331,
                          1345, 1359, 1373, 1429, 1443, 1457, 1471, 1485,
                          1513, 1534, 1569)] <- "no"

########
# 5436
########
tri_5436 <- aci_merged_tri %>%
  filter(id == 5436 & keep_row == "yes" & Ci > 380 & Asty > 0 & Ci < 1000) %>% 
  fitaci(varnames = list(ALEAF = "Asty", 
                         Tleaf = "Tleaf", 
                         Ci = "Ci",
                         PPFD = "Qin_cuvette"),
  Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5436)
summary(tri_5436)

# Write to data frame
aci_coefs <- data.frame(id = 5436, spp = "Tri", t(coef(tri_5436)), leaf_length = 5.7)

########
# 9145
########
tri_9145 <- aci_merged_tri %>%
  filter(id == 9145 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, fitmethod = "bilinear")
plot(tri_9145)
summary(tri_9145)

aci_coefs <- aci_coefs %>%
  add_row(id = 9145, spp = "Tri", 
          Vcmax = coef(tri_9145)[1], Jmax = coef(tri_9145)[2],
          Rd = coef(tri_9145)[3], TPU = coef(tri_9145)[4],
          leaf_length = 6.4)


########
# 2927
########
tri_2927 <- aci_merged_tri %>%
  filter(id == 2927 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_2927)
summary(tri_2927)

aci_coefs <- aci_coefs %>%
  add_row(id = 2927, spp = "Tri", 
          Vcmax = coef(tri_2927)[1], Jmax = coef(tri_2927)[2],
          Rd = coef(tri_2927)[3], TPU = coef(tri_2927)[4],
          leaf_length = 6.0)

########
# 2912
########
tri_2912 <- aci_merged_tri %>%
  filter(id == "2912_2", keep_row == "yes" & Ci > 100) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_2912)
summary(tri_2912)

aci_coefs <- aci_coefs %>%
  add_row(id = 2912, spp = "Tri", 
          Vcmax = coef(tri_2912)[1], Jmax = coef(tri_2912)[2],
          Rd = coef(tri_2912)[3], TPU = coef(tri_2912)[4],
          leaf_length = 5.0)

########
# 5381
########
tri_5381 <- aci_merged_tri %>%
  filter(id == 5381  & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_5381)
summary(tri_5381)

aci_coefs <- aci_coefs %>%
  add_row(id = 5381, spp = "Tri", 
          Vcmax = coef(tri_5381)[1], Jmax = coef(tri_5381)[2],
          Rd = coef(tri_5381)[3], TPU = coef(tri_5381)[4],
          leaf_length = 7.3)

########
# 5488 (no plant)
########
aci_coefs <- aci_coefs %>%
  add_row(id = 5488, spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 1374
########
tri_1374 <- aci_merged_tri %>%
  filter(id == 1374 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_1374)
summary(tri_1374)

aci_coefs <- aci_coefs %>%
  add_row(id = 1374, spp = "Tri", 
          Vcmax = coef(tri_1374)[1], Jmax = coef(tri_1374)[2],
          Rd = coef(tri_1374)[3], TPU = coef(tri_1374)[4],
          leaf_length = 6.5)

########
# 6875
########
tri_6875 <- aci_merged_tri %>%
  filter(id == 6875 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_6875)
summary(tri_6875)

aci_coefs <- aci_coefs %>%
  add_row(id = 6875, spp = "Tri", 
          Vcmax = coef(tri_6875)[1], Jmax = coef(tri_6875)[2],
          Rd = coef(tri_6875)[3], TPU = coef(tri_6875)[4],
          leaf_length = 7.1)

########
# 3031
########
tri_3031 <- aci_merged_tri %>%
  filter(id == "3031_2" & keep_row == "yes" & Ci < 1000) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3031)
summary(tri_3031)

aci_coefs <- aci_coefs %>%
  add_row(id = 3031, spp = "Tri", 
          Vcmax = coef(tri_3031)[1], Jmax = coef(tri_3031)[2],
          Rd = coef(tri_3031)[3], TPU = coef(tri_3031)[4],
          leaf_length = 6.5)

########
# 3004
########
tri_3004 <- aci_merged_tri %>%
  filter(id == 3004 & keep_row == "yes" & Ci > 0) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition = 400)
plot(tri_3004)
summary(tri_3004)

aci_coefs <- aci_coefs %>%
  add_row(id = 3004, spp = "Tri", 
          Vcmax = coef(tri_3004)[1], Jmax = coef(tri_3004)[2],
          Rd = coef(tri_3004)[3], TPU = coef(tri_3004)[4],
          leaf_length = 6.5)


########
# 552
########
tri_552 <- aci_merged_tri %>%
  filter(id == 552  & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_552)
summary(tri_552)

aci_coefs <- aci_coefs %>%
  add_row(id = 552, spp = "Tri", 
          Vcmax = coef(tri_552)[1], Jmax = coef(tri_552)[2],
          Rd = coef(tri_552)[3], TPU = coef(tri_552)[4],
          leaf_length = 5.9)

########
# 5904
########
tri_5904 <- aci_merged_tri %>%
  filter(id == 5904 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5904)
summary(tri_5904)

aci_coefs <- aci_coefs %>%
  add_row(id = 5904, spp = "Tri", 
          Vcmax = coef(tri_5904)[1], Jmax = coef(tri_5904)[2],
          Rd = coef(tri_5904)[3], TPU = coef(tri_5904)[4],
          leaf_length = 4.5)

########
# 5714
########
tri_5714 <- aci_merged_tri %>%
  filter(id == 5714 & keep_row == "yes" & Ci > 200 & Ci < 650) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5714)
summary(tri_5714)

aci_coefs <- aci_coefs %>%
  add_row(id = 5714, spp = "Tri", 
          Vcmax = coef(tri_5714)[1], Jmax = coef(tri_5714)[2],
          Rd = coef(tri_5714)[3], TPU = coef(tri_5714)[4],
          leaf_length = 5.3)

########
# 3379
########
aci_coefs <- aci_coefs %>%
  add_row(id = 3379, spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 4777
########
aci_coefs <- aci_coefs %>%
  add_row(id = 4777, spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 5722
########
tri_5722 <- aci_merged_tri %>%
  filter(id == 5722 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5722)
summary(tri_5722)

aci_coefs <- aci_coefs %>%
  add_row(id = 5722, spp = "Tri", 
          Vcmax = coef(tri_5722)[1], Jmax = coef(tri_5722)[2],
          Rd = coef(tri_5722)[3], TPU = coef(tri_5722)[4],
          leaf_length = 6.1)

########
# 2416 (update leaf area)
########
tri_2416 <- aci_merged_tri %>%
  filter(id == 2416 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2416)
summary(tri_2416)

aci_coefs <- aci_coefs %>%
  add_row(id = 2416, spp = "Tri", 
          Vcmax = coef(tri_2416)[1], Jmax = coef(tri_2416)[2],
          Rd = coef(tri_2416)[3], TPU = coef(tri_2416)[4],
          leaf_length = 4.4)

########
# 2276 (update leaf area)
########
tri_2276 <- aci_merged_tri %>%
  filter(id == 2276 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2276)
summary(tri_2276)

aci_coefs <- aci_coefs %>%
  add_row(id = 2276, spp = "Tri", 
          Vcmax = coef(tri_2276)[1], Jmax = coef(tri_2276)[2],
          Rd = coef(tri_2276)[3], TPU = coef(tri_2276)[4],
          leaf_length = 4.8)

########
# 1647 (update leaf area)
########
tri_1647 <- aci_merged_tri %>%
  filter(id == 1647 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1647)
summary(tri_1647)

aci_coefs <- aci_coefs %>%
  add_row(id = 1647, spp = "Tri", 
          Vcmax = coef(tri_1647)[1], Jmax = coef(tri_1647)[2],
          Rd = coef(tri_1647)[3], TPU = coef(tri_1647)[4],
          leaf_length = 5.5)

########
# 2563 (update leaf area)
########
tri_2563 <- aci_merged_tri %>%
  filter(id == 2563 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2563)
summary(tri_2563)

aci_coefs <- aci_coefs %>%
  add_row(id = 2563, spp = "Tri", 
          Vcmax = coef(tri_2563)[1], Jmax = coef(tri_2563)[2],
          Rd = coef(tri_2563)[3], TPU = coef(tri_2563)[4],
          leaf_length = 5.6)

########
# 4990 (update leaf area)
########
tri_4990 <- aci_merged_tri %>%
  filter(id == 4990 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4990)
summary(tri_4990)

aci_coefs <- aci_coefs %>%
  add_row(id = 4990, spp = "Tri", 
          Vcmax = coef(tri_4990)[1], Jmax = coef(tri_4990)[2],
          Rd = coef(tri_4990)[3], TPU = coef(tri_4990)[4],
          leaf_length = 7.4)

########
# 4177
########
tri_4177 <- aci_merged_tri %>%
  filter(id == 4177 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4177)
summary(tri_4177)

aci_coefs <- aci_coefs %>%
  add_row(id = 4177, spp = "Tri", 
          Vcmax = coef(tri_4177)[1], Jmax = coef(tri_4177)[2],
          Rd = coef(tri_4177)[3], TPU = coef(tri_4177)[4],
          leaf_length = 7.1)

########
# 989
########
tri_989 <- aci_merged_tri %>%
  filter(id == "989_2" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_989)
summary(tri_989)

aci_coefs <- aci_coefs %>%
  add_row(id = 989, spp = "Tri", 
          Vcmax = coef(tri_989)[1], Jmax = coef(tri_989)[2],
          Rd = coef(tri_989)[3], TPU = coef(tri_989)[4],
          leaf_length = 5.5)

########
# flag_1
########
tri_flag1 <- aci_merged_tri %>%
  filter(id == "flag1" & keep_row == "yes" & Ci < 950) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_flag1)
summary(tri_flag1)

aci_coefs <- aci_coefs %>%
  mutate(id = as.character(id)) %>%
  add_row(id = "flag1", spp = "Tri", 
          Vcmax = coef(tri_flag1)[1], Jmax = coef(tri_flag1)[2],
          Rd = coef(tri_flag1)[3], TPU = coef(tri_flag1)[4],
          leaf_length = 5.1)

########
# 5619
########
tri_5619 <- aci_merged_tri %>%
  filter(id == "5619" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5619)
summary(tri_5619)

aci_coefs <- aci_coefs %>%
  add_row(id = "5619", spp = "Tri", 
          Vcmax = coef(tri_5619)[1], Jmax = coef(tri_5619)[2],
          Rd = coef(tri_5619)[3], TPU = coef(tri_5619)[4],
          leaf_length = 6.1)

########
# 583
########
tri_583 <- aci_merged_tri %>%
  filter(id == 583 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_583)
summary(tri_583)

aci_coefs <- aci_coefs %>%
  add_row(id = "583", spp = "Tri", 
          Vcmax = coef(tri_583)[1], Jmax = coef(tri_583)[2],
          Rd = coef(tri_583)[3], TPU = coef(tri_583)[4],
          leaf_length = 6.0)

########
# 6507
########
tri_6507 <- aci_merged_tri %>%
  filter(id == 6507 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6507)
summary(tri_6507)

aci_coefs <- aci_coefs %>%
  add_row(id = "6507", spp = "Tri", 
          Vcmax = coef(tri_6507)[1], Jmax = coef(tri_6507)[2],
          Rd = coef(tri_6507)[3], TPU = coef(tri_6507)[4],
          leaf_length = 5.9)

########
# 5877
########
tri_5877 <- aci_merged_tri %>%
  filter(id == 5877 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE, citransition =  400)
plot(tri_5877)
summary(tri_5877)

aci_coefs <- aci_coefs %>%
  add_row(id = "5877", spp = "Tri", 
          Vcmax = coef(tri_5877)[1], Jmax = coef(tri_5877)[2],
          Rd = coef(tri_5877)[3], TPU = coef(tri_5877)[4],
          leaf_length = 7.2)

########
# 5865
########
tri_5865 <- aci_merged_tri %>%
  filter(id == 5865 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5865)
summary(tri_5865)

aci_coefs <- aci_coefs %>%
  add_row(id = "5865", spp = "Tri", 
          Vcmax = coef(tri_5865)[1], Jmax = coef(tri_5865)[2],
          Rd = coef(tri_5865)[3], TPU = coef(tri_5865)[4],
          leaf_length = 4.8)

########
# 1851
########
tri_1851 <- aci_merged_tri %>%
  filter(id == "1851_2"  & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1851)
summary(tri_1851)

aci_coefs <- aci_coefs %>%
  add_row(id = "1851", spp = "Tri", 
          Vcmax = coef(tri_1851)[1], Jmax = coef(tri_1851)[2],
          Rd = coef(tri_1851)[3], TPU = coef(tri_1851)[4],
          leaf_length = 5.3)

########
# 43
########
tri_43 <- aci_merged_tri %>%
  filter(id == 43 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_43)
summary(tri_43)

aci_coefs <- aci_coefs %>%
  add_row(id = "43", spp = "Tri", 
          Vcmax = coef(tri_43)[1], Jmax = coef(tri_43)[2],
          Rd = coef(tri_43)[3], TPU = coef(tri_43)[4],
          leaf_length = 7.1)

########
# 4760
########
tri_4760 <- aci_merged_tri %>%
  filter(id == 4760 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4760)
summary(tri_4760)

aci_coefs <- aci_coefs %>%
  add_row(id = "4760", spp = "Tri", 
          Vcmax = coef(tri_4760)[1], Jmax = coef(tri_4760)[2],
          Rd = coef(tri_4760)[3], TPU = coef(tri_4760)[4],
          leaf_length = 4.9)

########
# 7120
########
aci_coefs <- aci_coefs %>%
  add_row(id = "7120", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 6879
########
tri_6879 <- aci_merged_tri %>%
  filter(id == 6879 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6879)
summary(tri_6879)

aci_coefs <- aci_coefs %>%
  add_row(id = "6879", spp = "Tri", 
          Vcmax = coef(tri_6879)[1], Jmax = coef(tri_6879)[2],
          Rd = coef(tri_6879)[3], TPU = coef(tri_6879)[4],
          leaf_length = 4.7)

########
# 5403
########
aci_coefs <- aci_coefs %>%
  add_row(id = "5403", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 3305
########
tri_3305 <- aci_merged_tri %>%
  filter(id == 3305 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3305)
summary(tri_3305)

aci_coefs <- aci_coefs %>%
  add_row(id = "3305", spp = "Tri", 
          Vcmax = coef(tri_3305)[1], Jmax = coef(tri_3305)[2],
          Rd = coef(tri_3305)[3], TPU = coef(tri_3305)[4],
          leaf_length = 4.4)

########
# 5040
########
tri_5040 <- aci_merged_tri %>%
  filter(id == 5040 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5040)
summary(tri_5040)

aci_coefs <- aci_coefs %>%
  add_row(id = "5040", spp = "Tri", 
          Vcmax = coef(tri_5040)[1], Jmax = coef(tri_5040)[2],
          Rd = coef(tri_5040)[3], TPU = coef(tri_5040)[4],
          leaf_length = 7.0)

########
# 902
########
tri_902 <- aci_merged_tri %>%
  filter(id == 902 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_902)
summary(tri_902)

aci_coefs <- aci_coefs %>%
  add_row(id = "902", spp = "Tri", 
          Vcmax = coef(tri_902)[1], Jmax = coef(tri_902)[2],
          Rd = coef(tri_902)[3], TPU = coef(tri_902)[4],
          leaf_length = 7.2)

########
# 5852
########
tri_5852 <- aci_merged_tri %>%
  filter(id == 5852 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5852)
summary(tri_5852)

aci_coefs <- aci_coefs %>%
  add_row(id = "5852", spp = "Tri", 
          Vcmax = coef(tri_5852)[1], Jmax = coef(tri_5852)[2],
          Rd = coef(tri_5852)[3], TPU = coef(tri_5852)[4],
          leaf_length = 5.4)

########
# 4444(2988)
########
tri_4444 <- aci_merged_tri %>%
  filter(id == "4444(2988)" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4444)
summary(tri_4444)

aci_coefs <- aci_coefs %>%
  add_row(id = "4444(2988)", spp = "Tri", 
          Vcmax = coef(tri_4444)[1], Jmax = coef(tri_4444)[2],
          Rd = coef(tri_4444)[3], TPU = coef(tri_4444)[4],
          leaf_length = 7.5)

########
# 5060
########
tri_5060 <- aci_merged_tri %>%
  filter(id ==  "5060" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5060)
summary(tri_5060)

aci_coefs <- aci_coefs %>%
  add_row(id = "5060", spp = "Tri", 
          Vcmax = coef(tri_5060)[1], Jmax = coef(tri_5060)[2],
          Rd = coef(tri_5060)[3], TPU = coef(tri_5060)[4],
          leaf_length = 5.8)

########
# 916
########
tri_916 <- aci_merged_tri %>%
  filter(id == 916 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_916)
summary(tri_916)

aci_coefs <- aci_coefs %>%
  add_row(id = "916", spp = "Tri", 
          Vcmax = coef(tri_916)[1], Jmax = coef(tri_916)[2],
          Rd = coef(tri_916)[3], TPU = coef(tri_916)[4],
          leaf_length = 4.6)

########
# 5031
########
aci_coefs <- aci_coefs %>%
  add_row(id = "5031", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 2980
########
tri_2980 <- aci_merged_tri %>%
  filter(id == 2980 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2980)
summary(tri_2980)

aci_coefs <- aci_coefs %>%
  add_row(id = "2980", spp = "Tri", 
          Vcmax = coef(tri_2980)[1], Jmax = coef(tri_2980)[2],
          Rd = coef(tri_2980)[3], TPU = coef(tri_2980)[4],
          leaf_length = 5.6)

########
# 774
########
tri_774 <- aci_merged_tri %>%
  filter(id == 774 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_774)
summary(tri_774)

aci_coefs <- aci_coefs %>%
  add_row(id = "774", spp = "Tri", 
          Vcmax = coef(tri_774)[1], Jmax = coef(tri_774)[2],
          Rd = coef(tri_774)[3], TPU = coef(tri_774)[4],
          leaf_length = 6.5)

########
# 4777c
########
tri_4777c <- aci_merged_tri %>%
  filter(id == "4774C" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4777c)
summary(tri_4777c)

aci_coefs <- aci_coefs %>%
  add_row(id = "4777C", spp = "Tri", 
          Vcmax = coef(tri_4777c)[1], Jmax = coef(tri_4777c)[2],
          Rd = coef(tri_4777c)[3], TPU = coef(tri_4777c)[4],
          leaf_length = 5.9)

########
# 1827
########
tri_1827 <- aci_merged_tri %>%
  filter(id == 1827 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1827)
summary(tri_1827)

aci_coefs <- aci_coefs %>%
  add_row(id = "1827", spp = "Tri", 
          Vcmax = coef(tri_1827)[1], Jmax = coef(tri_1827)[2],
          Rd = coef(tri_1827)[3], TPU = coef(tri_1827)[4],
          leaf_length = 5.5)

########
# 5771
########
tri_5771 <- aci_merged_tri %>%
  filter(id == 5771 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5771)
summary(tri_5771)

aci_coefs <- aci_coefs %>%
  add_row(id = "5771", spp = "Tri", 
          Vcmax = coef(tri_5771)[1], Jmax = coef(tri_5771)[2],
          Rd = coef(tri_5771)[3], TPU = coef(tri_5771)[4],
          leaf_length = 5.5)

########
# 3371
########
tri_3371 <- aci_merged_tri %>%
  filter(id == 3371 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3371)
summary(tri_3371)

aci_coefs <- aci_coefs %>%
  add_row(id = "3371", spp = "Tri", 
          Vcmax = coef(tri_3371)[1], Jmax = coef(tri_3371)[2],
          Rd = coef(tri_3371)[3], TPU = coef(tri_3371)[4],
          leaf_length = 6)

########
# 3708
########
tri_3708 <- aci_merged_tri %>%
  filter(id == 3708 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3708)
summary(tri_3708)

aci_coefs <- aci_coefs %>%
  add_row(id = "3708", spp = "Tri", 
          Vcmax = coef(tri_3708)[1], Jmax = coef(tri_3708)[2],
          Rd = coef(tri_3708)[3], TPU = coef(tri_3708)[4],
          leaf_length = 6.6)

########
# 5229
########
tri_5229red <- aci_merged_tri %>%
  filter(id == 5229 & keep_row == "yes" & machine == "gibson") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5229red)
summary(tri_5229red)

aci_coefs <- aci_coefs %>%
  add_row(id = "5229_red", spp = "Tri", 
          Vcmax = coef(tri_5229red)[1], Jmax = coef(tri_5229red)[2],
          Rd = coef(tri_5229red)[3], TPU = coef(tri_5229red)[4],
          leaf_length = 6.2)

########
# 4149
########
tri_4149 <- aci_merged_tri %>%
  filter(id == 4149 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4149)
summary(tri_4149)

aci_coefs <- aci_coefs %>%
  add_row(id = "4149", spp = "Tri", 
          Vcmax = coef(tri_4149)[1], Jmax = coef(tri_4149)[2],
          Rd = coef(tri_4149)[3], TPU = coef(tri_4149)[4],
          leaf_length = 7.1)

########
# 4000
########
tri_4000 <- aci_merged_tri %>%
  filter(id ==  4000 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4000)
summary(tri_4000)

aci_coefs <- aci_coefs %>%
  add_row(id = "4000", spp = "Tri", 
          Vcmax = coef(tri_4000)[1], Jmax = coef(tri_4000)[2],
          Rd = coef(tri_4000)[3], TPU = coef(tri_4000)[4],
          leaf_length = 7.5)

########
# 2988
########
tri_2988 <- aci_merged_tri %>%
  filter(id == 2988 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2988)
summary(tri_2988)

aci_coefs <- aci_coefs %>%
  add_row(id = "2988", spp = "Tri", 
          Vcmax = coef(tri_2988)[1], Jmax = coef(tri_2988)[2],
          Rd = coef(tri_2988)[3], TPU = coef(tri_2988)[4],
          leaf_length = 6.1)

########
# 452
########
tri_452 <- aci_merged_tri %>%
  filter(id == 452 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_452)
summary(tri_452)

aci_coefs <- aci_coefs %>%
  add_row(id = "452", spp = "Tri", 
          Vcmax = coef(tri_452)[1], Jmax = coef(tri_452)[2],
          Rd = coef(tri_452)[3], TPU = coef(tri_452)[4],
          leaf_length = 5.5)

########
# 6558
########
tri_6558 <- aci_merged_tri %>%
  filter(id == 6558 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6558)
summary(tri_6558)

aci_coefs <- aci_coefs %>%
  add_row(id = "6558", spp = "Tri", 
          Vcmax = coef(tri_6558)[1], Jmax = coef(tri_6558)[2],
          Rd = coef(tri_6558)[3], TPU = coef(tri_6558)[4],
          leaf_length = 6.5)

########
# 4505
########
tri_4505 <- aci_merged_tri %>%
  filter(id == 4505 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4505)
summary(tri_4505)

aci_coefs <- aci_coefs %>%
  add_row(id = "4505", spp = "Tri", 
          Vcmax = coef(tri_4505)[1], Jmax = coef(tri_4505)[2],
          Rd = coef(tri_4505)[3], TPU = coef(tri_4505)[4],
          leaf_length = 5.9)

########
# 2329
########
tri_2329 <- aci_merged_tri %>%
  filter(id == 2329 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2329)
summary(tri_2329)

aci_coefs <- aci_coefs %>%
  add_row(id = "2329", spp = "Tri", 
          Vcmax = coef(tri_2329)[1], Jmax = coef(tri_2329)[2],
          Rd = coef(tri_2329)[3], TPU = coef(tri_2329)[4],
          leaf_length = 7.8)

########
# 2996
########
tri_2996a2 <- aci_merged_tri %>%
  filter(id == "2996A2" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2996a2)
summary(tri_2996a2)

aci_coefs <- aci_coefs %>%
  add_row(id = "2996A2", spp = "Tri", 
          Vcmax = coef(tri_2996a2)[1], Jmax = coef(tri_2996a2)[2],
          Rd = coef(tri_2996a2)[3], TPU = coef(tri_2996a2)[4],
          leaf_length = 5)

########
# 632
########
tri_632 <- aci_merged_tri %>%
  filter(id == 632 & keep_row == "yes" & Ci < 1100) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_632)
summary(tri_632)

aci_coefs <- aci_coefs %>%
  add_row(id = "632", spp = "Tri", 
          Vcmax = coef(tri_632)[1], Jmax = coef(tri_632)[2],
          Rd = coef(tri_632)[3], TPU = coef(tri_632)[4],
          leaf_length = 5)

########
# 5228
########
tri_5228 <- aci_merged_tri %>%
  filter(id == 5228 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5228)
summary(tri_5228)

aci_coefs <- aci_coefs %>%
  add_row(id = "5228", spp = "Tri", 
          Vcmax = coef(tri_5228)[1], Jmax = coef(tri_5228)[2],
          Rd = coef(tri_5228)[3], TPU = coef(tri_5228)[4],
          leaf_length = 7.1)

########
# 4712
########
tri_4712 <- aci_merged_tri %>%
  filter(id == 4712 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4712)
summary(tri_4712)

aci_coefs <- aci_coefs %>%
  add_row(id = "4712", spp = "Tri", 
          Vcmax = coef(tri_4712)[1], Jmax = coef(tri_4712)[2],
          Rd = coef(tri_4712)[3], TPU = coef(tri_4712)[4],
          leaf_length = 5.9)

########
# 6894
########
tri_6894 <- aci_merged_tri %>%
  filter(id == 6894 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6894)
summary(tri_6894)

aci_coefs <- aci_coefs %>%
  add_row(id = "6894", spp = "Tri", 
          Vcmax = coef(tri_6894)[1], Jmax = coef(tri_6894)[2],
          Rd = coef(tri_6894)[3], TPU = coef(tri_6894)[4],
          leaf_length = 2.1)

########
# 1669
########
tri_1669 <- aci_merged_tri %>%
  filter(id == 1669 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1669)
summary(tri_1669)

aci_coefs <- aci_coefs %>%
  add_row(id = "1669", spp = "Tri", 
          Vcmax = coef(tri_1669)[1], Jmax = coef(tri_1669)[2],
          Rd = coef(tri_1669)[3], TPU = coef(tri_1669)[4],
          leaf_length = 6.1)

########
# 978
########
tri_978 <- aci_merged_tri %>%
  filter(id == 978 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_978)
summary(tri_978)

aci_coefs <- aci_coefs %>%
  add_row(id = "978", spp = "Tri", 
          Vcmax = coef(tri_978)[1], Jmax = coef(tri_978)[2],
          Rd = coef(tri_978)[3], TPU = coef(tri_978)[4],
          leaf_length = 4.8)

########
# 2382
########
tri_2382 <- aci_merged_tri %>%
  filter(id == 2382 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2382)
summary(tri_2382)

aci_coefs <- aci_coefs %>%
  add_row(id = "2382", spp = "Tri", 
          Vcmax = coef(tri_2382)[1], Jmax = coef(tri_2382)[2],
          Rd = coef(tri_2382)[3], TPU = coef(tri_2382)[4],
          leaf_length = 4.9)

########
# flag3
########
tri_flag3 <- aci_merged_tri %>%
  filter(id == "flag3" & keep_row == "yes" & Asty >-0.2) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_flag3)
summary(tri_flag3)

aci_coefs <- aci_coefs %>%
  add_row(id = "flag3", spp = "Tri", 
          Vcmax = coef(tri_flag3)[1], Jmax = coef(tri_flag3)[2],
          Rd = coef(tri_flag3)[3], TPU = coef(tri_flag3)[4],
          leaf_length = 5.7)

########
# 6885
########
tri_6885 <- aci_merged_tri %>%
  filter(id == "6885" & keep_row == "yes" & Ci < 1000) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6885)
summary(tri_6885)

aci_coefs <- aci_coefs %>%
  add_row(id = "6885", spp = "Tri", 
          Vcmax = coef(tri_6885)[1], Jmax = coef(tri_6885)[2],
          Rd = coef(tri_6885)[3], TPU = coef(tri_6885)[4],
          leaf_length = 4.9)

########
# 4109
########
tri_4109 <- aci_merged_tri %>%
  filter(id == 4109 & keep_row == "yes" & Ci > 300 & Ci < 1100) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4109)
summary(tri_4109)

aci_coefs <- aci_coefs %>%
  add_row(id = "4109", spp = "Tri", 
          Vcmax = coef(tri_4109)[1], Jmax = coef(tri_4109)[2],
          Rd = coef(tri_4109)[3], TPU = coef(tri_4109)[4],
          leaf_length = 4.8)

########
# 4714
########
tri_4714 <- aci_merged_tri %>%
  filter(id == 4714 & keep_row == "yes" & Ci > 100) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4714)
summary(tri_4714)

aci_coefs <- aci_coefs %>%
  add_row(id = "4714", spp = "Tri", 
          Vcmax = coef(tri_4714)[1], Jmax = coef(tri_4714)[2],
          Rd = coef(tri_4714)[3], TPU = coef(tri_4714)[4],
          leaf_length = 6.1)

########
# 3829
########
tri_3829 <- aci_merged_tri %>%
  filter(id == 3829 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3829)
summary(tri_3829)

aci_coefs <- aci_coefs %>%
  add_row(id = "3829", spp = "Tri", 
          Vcmax = coef(tri_3829)[1], Jmax = coef(tri_3829)[2],
          Rd = coef(tri_3829)[3], TPU = coef(tri_3829)[4],
          leaf_length = 4.6)

########
# 482
########
tri_482 <- aci_merged_tri %>%
  filter(id == 482 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_482)
summary(tri_482)

aci_coefs <- aci_coefs %>%
  add_row(id = "482", spp = "Tri", 
          Vcmax = coef(tri_482)[1], Jmax = coef(tri_482)[2],
          Rd = coef(tri_482)[3], TPU = coef(tri_482)[4],
          leaf_length = 6.4)

########
# 2074
########
aci_coefs <- aci_coefs %>%
  add_row(id = "2074", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 4511
########
tri_4511 <- aci_merged_tri %>%
  filter(id == 4511 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4511)
summary(tri_4511)

aci_coefs <- aci_coefs %>%
  add_row(id = "4511", spp = "Tri", 
          Vcmax = coef(tri_4511)[1], Jmax = coef(tri_4511)[2],
          Rd = coef(tri_4511)[3], TPU = coef(tri_4511)[4],
          leaf_length = 8.1)

########
# 5229
########
tri_5229blue <- aci_merged_tri %>%
  filter(id == 5229 & keep_row == "yes" & is.na(machine)) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5229blue)
summary(tri_5229blue)

aci_coefs <- aci_coefs %>%
  add_row(id = "5229_blue", spp = "Tri", 
          Vcmax = coef(tri_5229blue)[1], Jmax = coef(tri_5229blue)[2],
          Rd = coef(tri_5229blue)[3], TPU = coef(tri_5229blue)[4],
          leaf_length = 6.4)

########
# 392
########
tri_392 <- aci_merged_tri %>%
  filter(id == 392 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_392)
summary(tri_392)

aci_coefs <- aci_coefs %>%
  add_row(id = "392", spp = "Tri", 
          Vcmax = coef(tri_392)[1], Jmax = coef(tri_392)[2],
          Rd = coef(tri_392)[3], TPU = coef(tri_392)[4],
          leaf_length = 5.2)

########
# 2996A3
########
tri_2996a3 <- aci_merged_tri %>%
  filter(id == "2996A3" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2996a3)
summary(tri_2996a3)

aci_coefs <- aci_coefs %>%
  add_row(id = "2996A3", spp = "Tri", 
          Vcmax = coef(tri_2996a3)[1], Jmax = coef(tri_2996a3)[2],
          Rd = coef(tri_2996a3)[3], TPU = coef(tri_2996a3)[4],
          leaf_length = 7.0)

########
# 1750
########
tri_1750 <- aci_merged_tri %>%
  filter(id == 1750 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1750)
summary(tri_1750)

aci_coefs <- aci_coefs %>%
  add_row(id = "1750", spp = "Tri", 
          Vcmax = coef(tri_1750)[1], Jmax = coef(tri_1750)[2],
          Rd = coef(tri_1750)[3], TPU = coef(tri_1750)[4],
          leaf_length = 6.4)

########
# 4265
########
tri_4265 <- aci_merged_tri %>%
  filter(id == 4265 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4265)
summary(tri_4265)

aci_coefs <- aci_coefs %>%
  add_row(id = "4265", spp = "Tri", 
          Vcmax = coef(tri_4265)[1], Jmax = coef(tri_4265)[2],
          Rd = coef(tri_4265)[3], TPU = coef(tri_4265)[4],
          leaf_length = 6.4)

########
# 4576
########
tri_4576 <- aci_merged_tri %>%
  filter(id == 4576 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4576)
summary(tri_4576)

aci_coefs <- aci_coefs %>%
  add_row(id = "4576", spp = "Tri", 
          Vcmax = coef(tri_4576)[1], Jmax = coef(tri_4576)[2],
          Rd = coef(tri_4576)[3], TPU = coef(tri_4576)[4],
          leaf_length = 6.6)

########
# 5742
########
aci_coefs <- aci_coefs %>%
  add_row(id = "5742", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = 5.6)

########
# 1184
########
aci_coefs <- aci_coefs %>%
  add_row(id = "1184", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 2996A1
########
tri_2996a1 <- aci_merged_tri %>%
  filter(id == "2996A1" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2996a1)
summary(tri_2996a1)

aci_coefs <- aci_coefs %>%
  add_row(id = "2996A1", spp = "Tri", 
          Vcmax = coef(tri_2996a1)[1], Jmax = coef(tri_2996a1)[2],
          Rd = coef(tri_2996a1)[3], TPU = coef(tri_2996a1)[4],
          leaf_length = NA)

########
# 1666
########
tri_1666 <- aci_merged_tri %>%
  filter(id == 1666 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1666)
summary(tri_1666)

aci_coefs <- aci_coefs %>%
  add_row(id = "1666", spp = "Tri", 
          Vcmax = coef(tri_1666)[1], Jmax = coef(tri_1666)[2],
          Rd = coef(tri_1666)[3], TPU = coef(tri_1666)[4],
          leaf_length = 5.1)

########
# 4942
########
tri_4942 <- aci_merged_tri %>%
  filter(id == 4942 & keep_row == "yes" & Asty > -0.6) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4942)
summary(tri_4942)

aci_coefs <- aci_coefs %>%
  add_row(id = "4942", spp = "Tri", 
          Vcmax = coef(tri_4942)[1], Jmax = coef(tri_4942)[2],
          Rd = coef(tri_4942)[3], TPU = coef(tri_4942)[4],
          leaf_length = 5.3)

########
# flag2
########
tri_flag2 <- aci_merged_tri %>%
  filter(id == "flag2" & keep_row == "yes" & Asty > -0.4) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_flag2)
summary(tri_flag2)

aci_coefs <- aci_coefs %>%
  add_row(id = "flag2", spp = "Tri", 
          Vcmax = coef(tri_flag2)[1], Jmax = coef(tri_flag2)[2],
          Rd = coef(tri_flag2)[3], TPU = coef(tri_flag2)[4],
          leaf_length = 4)

########
# 972
########
tri_972 <- aci_merged_tri %>%
  filter(id == 972 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_972)
summary(tri_972)

aci_coefs <- aci_coefs %>%
  add_row(id = "972", spp = "Tri", 
          Vcmax = coef(tri_972)[1], Jmax = coef(tri_972)[2],
          Rd = coef(tri_972)[3], TPU = coef(tri_972)[4],
          leaf_length = 6.5)

########
# 3043
########
tri_3043 <- aci_merged_tri %>%
  filter(id == 3043 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3043)
summary(tri_3043)

aci_coefs <- aci_coefs %>%
  add_row(id = "3043", spp = "Tri", 
          Vcmax = coef(tri_3043)[1], Jmax = coef(tri_3043)[2],
          Rd = coef(tri_3043)[3], TPU = coef(tri_3043)[4],
          leaf_length = 5.6)

########
# 4414
########
tri_4414 <- aci_merged_tri %>%
  filter(id == 4414 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4414)
summary(tri_4414)

aci_coefs <- aci_coefs %>%
  add_row(id = "4414", spp = "Tri", 
          Vcmax = coef(tri_4414)[1], Jmax = coef(tri_4414)[2],
          Rd = coef(tri_4414)[3], TPU = coef(tri_4414)[4],
          leaf_length = 5.5)

########
# 4105
########
tri_4105 <- aci_merged_tri %>%
  filter(id == 4105 & keep_row == "yes" & Asty > -0.4) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4105)
summary(tri_4105)

aci_coefs <- aci_coefs %>%
  add_row(id = "4105", spp = "Tri", 
          Vcmax = coef(tri_4105)[1], Jmax = coef(tri_4105)[2],
          Rd = coef(tri_4105)[3], TPU = coef(tri_4105)[4],
          leaf_length = 5.4)

########
# 3563
########
tri_3563 <- aci_merged_tri %>%
  filter(id == 3563 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3563)
summary(tri_3563)

aci_coefs <- aci_coefs %>%
  add_row(id = "3563", spp = "Tri", 
          Vcmax = coef(tri_3563)[1], Jmax = coef(tri_3563)[2],
          Rd = coef(tri_3563)[3], TPU = coef(tri_3563)[4],
          leaf_length = 4.6)

########
# 2886
########
aci_coefs <- aci_coefs %>%
  add_row(id = "2886", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 4431
########
aci_coefs <- aci_coefs %>%
  add_row(id = "4431", spp = "Tri", 
          Vcmax = NA, Jmax = NA,
          Rd = NA, TPU = NA,
          leaf_length = NA)

########
# 7574
########
tri_7574 <- aci_merged_tri %>%
  filter(id == 7574 & keep_row == "yes" & Asty < 4) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_7574)
summary(tri_7574)

aci_coefs <- aci_coefs %>%
  add_row(id = "7574", spp = "Tri", 
          Vcmax = coef(tri_7574)[1], Jmax = coef(tri_7574)[2],
          Rd = coef(tri_7574)[3], TPU = coef(tri_7574)[4],
          leaf_length = 4.7)

########
# 4959
########
tri_4959 <- aci_merged_tri %>%
  filter(id == 4959 & keep_row == "yes" & Asty > -0.8) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4959)
summary(tri_4959)

aci_coefs <- aci_coefs %>%
  add_row(id = "4959", spp = "Tri", 
          Vcmax = coef(tri_4959)[1], Jmax = coef(tri_4959)[2],
          Rd = coef(tri_4959)[3], TPU = coef(tri_4959)[4],
          leaf_length = 6.6)

########
# 4081
########
tri_4081 <- aci_merged_tri %>%
  filter(id == 4081 & keep_row == "yes" & Asty < 1 & Asty > -0.5) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4081)
summary(tri_4081)

aci_coefs <- aci_coefs %>%
  add_row(id = "4081", spp = "Tri", 
          Vcmax = coef(tri_4081)[1], Jmax = coef(tri_4081)[2],
          Rd = coef(tri_4081)[3], TPU = coef(tri_4081)[4],
          leaf_length = 4.8)

########
# 955
########
tri_955 <- aci_merged_tri %>%
  filter(id == 955 & keep_row == "yes" & elapsed > 6000) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_955)
summary(tri_955)

aci_coefs <- aci_coefs %>%
  add_row(id = "955", spp = "Tri", 
          Vcmax = coef(tri_955)[1], Jmax = coef(tri_955)[2],
          Rd = coef(tri_955)[3], TPU = coef(tri_955)[4],
          leaf_length = 6.5)

########
# 4452
########
tri_4452 <- aci_merged_tri %>%
  filter(id == 4452 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4452)
summary(tri_4452)

aci_coefs <- aci_coefs %>%
  add_row(id = "4452", spp = "Tri", 
          Vcmax = coef(tri_4452)[1], Jmax = coef(tri_4452)[2],
          Rd = coef(tri_4452)[3], TPU = coef(tri_4452)[4],
          leaf_length = 6.9)

########
# 3077
########
tri_3077 <- aci_merged_tri %>%
  filter(id == 3077 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3077)
summary(tri_3077)

aci_coefs <- aci_coefs %>%
  add_row(id = "3077", spp = "Tri", 
          Vcmax = coef(tri_3077)[1], Jmax = coef(tri_3077)[2],
          Rd = coef(tri_3077)[3], TPU = coef(tri_3077)[4],
          leaf_length = 5.5)

########
# 1975
########
tri_1975 <- aci_merged_tri %>%
  filter(id == 1975 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1975)
summary(tri_1975)

aci_coefs <- aci_coefs %>%
  add_row(id = "1975", spp = "Tri", 
          Vcmax = coef(tri_1975)[1], Jmax = coef(tri_1975)[2],
          Rd = coef(tri_1975)[3], TPU = coef(tri_1975)[4],
          leaf_length = 6.5)

########
# 2941
########
tri_2941 <- aci_merged_tri %>%
  filter(id == 2941 & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_2941)
summary(tri_2941)

aci_coefs <- aci_coefs %>%
  add_row(id = "2941", spp = "Tri", 
          Vcmax = coef(tri_2941)[1], Jmax = coef(tri_2941)[2],
          Rd = coef(tri_2941)[3], TPU = coef(tri_2941)[4],
          leaf_length = 6.3)

########
# 1926
########
tri_1926 <- aci_merged_tri %>%
  filter(id == 1926 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_1926)
summary(tri_1926)

aci_coefs <- aci_coefs %>%
  add_row(id = "1926", spp = "Tri", 
          Vcmax = coef(tri_1926)[1], Jmax = coef(tri_1926)[2],
          Rd = coef(tri_1926)[3], TPU = coef(tri_1926)[4],
          leaf_length = 6.1)

########
# 4745
########
tri_4745 <- aci_merged_tri %>%
  filter(id == 4745 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_4745)
summary(tri_4745)

aci_coefs <- aci_coefs %>%
  add_row(id = "4745", spp = "Tri", 
          Vcmax = coef(tri_4745)[1], Jmax = coef(tri_4745)[2],
          Rd = coef(tri_4745)[3], TPU = coef(tri_4745)[4],
          leaf_length = 5.2)

########
# 3526
########
tri_3526 <- aci_merged_tri %>%
  filter(id == 3526 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_3526)
summary(tri_3526)

aci_coefs <- aci_coefs %>%
  add_row(id = "3526", spp = "Tri", 
          Vcmax = coef(tri_3526)[1], Jmax = coef(tri_3526)[2],
          Rd = coef(tri_3526)[3], TPU = coef(tri_3526)[4],
          leaf_length = 5.2)

########
# 5626
########
tri_5626 <- aci_merged_tri %>%
  filter(id == 5626 & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5626)
summary(tri_5626)

aci_coefs <- aci_coefs %>%
  add_row(id = "5626", spp = "Tri", 
          Vcmax = coef(tri_5626)[1], Jmax = coef(tri_5626)[2],
          Rd = coef(tri_5626)[3], TPU = coef(tri_5626)[4],
          leaf_length = 5.9)

########
# 6881
########
tri_6881 <- aci_merged_tri %>%
  filter(id == 6881 & keep_row == "yes" & Ci > 200 & Ci < 1200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_6881)
summary(tri_6881)

aci_coefs <- aci_coefs %>%
  add_row(id = "6881", spp = "Tri", 
          Vcmax = coef(tri_6881)[1], Jmax = coef(tri_6881)[2],
          Rd = coef(tri_6881)[3], TPU = coef(tri_6881)[4],
          leaf_length = 4.3)

########
# 7147
########
tri_7147 <- aci_merged_tri %>%
  filter(id == "7147_2" & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_7147)
summary(tri_7147)

aci_coefs <- aci_coefs %>%
  add_row(id = "7147", spp = "Tri", 
          Vcmax = coef(tri_7147)[1], Jmax = coef(tri_7147)[2],
          Rd = coef(tri_7147)[3], TPU = coef(tri_7147)[4],
          leaf_length = 6.6)

########
# 5267
########
tri_5267 <- aci_merged_tri %>%
  filter(id == "5267_2" & keep_row == "yes" & Ci > 200) %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot(tri_5267)
summary(tri_5267)

aci_coefs <- aci_coefs %>%
  add_row(id = "5267", spp = "Tri", 
          Vcmax = coef(tri_5267)[1], Jmax = coef(tri_5267)[2],
          Rd = coef(tri_5267)[3], TPU = coef(tri_5267)[4],
          leaf_length = 5.1)

########
# 
########
tri_ <- aci_merged_tri %>%
  filter(id ==  & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot()
summary()

aci_coefs <- aci_coefs %>%
  add_row(id = "", spp = "Tri", 
          Vcmax = coef()[1], Jmax = coef()[2],
          Rd = coef()[3], TPU = coef()[4],
          leaf_length = )

########
# 
########
tri_ <- aci_merged_tri %>%
  filter(id ==  & keep_row == "yes") %>% fitaci(
    varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                    Ci = "Ci", PPFD = "Qin_cuvette"),
    Tcorrect = FALSE, fitTPU = TRUE)
plot()
summary()

aci_coefs <- aci_coefs %>%
  add_row(id = "", spp = "Tri", 
          Vcmax = coef()[1], Jmax = coef()[2],
          Rd = coef()[3], TPU = coef()[4],
          leaf_length = )





