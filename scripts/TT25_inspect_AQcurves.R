# Libraries
library(tidyverse)
library(ggpubr)

# Load helper function for cleaning licor files
source("../functions/clean_licor_files.R")

# Clean A/Q curve files
clean_licor_files(directory_path = "../data/li6800_lightResponse/raw/",
                  write_directory = "../data/li6800_lightResponse/clean/")

# Locate folder where cleaned A/Q curve files are located
files <- list.files(path = "../data/li6800_lightResponse/clean", 
                    pattern = "\\.csv$",
                    recursive = TRUE,
                    full.names = TRUE)
files <- setNames(files, stringr::str_extract(basename(files),
                                              ".*(?=\\.csv)"))

# Merge A/Q curve files into data frame
aq_merged <- plyr::rbind.fill(lapply(files, read.csv)) 

# Inspect merged data frame
head(aq_merged)

# Plot A/Q curves
five105_mai_plot <- ggplot(data = subset(aq_merged, id == "5105_mai_orange" & Qin < 602), 
                           aes(x = Qin, y = Asty)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 300, color = "red", linewidth = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
  scale_y_continuous(limits = c(0.6, 1.4), breaks = seq(0.6, 1.4, 0.2)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")), 
       y = expression(bold("A"["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  theme_classic(base_size = 18)

six483_mai_plot <- ggplot(data = subset(aq_merged, id == "6483_mai_yellow" & Qin < 602), 
                          aes(x = Qin, y = A)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 300, color = "red", linewidth = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
  scale_y_continuous(limits = c(1.5, 3.5), breaks = seq(1.5, 3.5, 0.5)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")), 
       y = expression(bold("A"["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  theme_classic(base_size = 18)

TT24_101_mai_plot <- ggplot(data = subset(aq_merged, 
                                          id == "TT24_101_mai_pink" & Qin < 602), 
                            aes(x = Qin, y = A)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 300, color = "red", linewidth = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
  scale_y_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.15)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")), 
       y = expression(bold("A"["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  theme_classic(base_size = 18)

fifty600_mai_plot <- ggplot(data = subset(aq_merged, id == "5600_mai_red" & Qin < 602), 
                            aes(x = Qin, y = A)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 300, color = "red", linewidth = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
  scale_y_continuous(limits = c(1.0, 2.0), breaks = seq(1.0, 2.0, 0.25)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")), 
       y = expression(bold("A"["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  theme_classic(base_size = 18)

four781_mai_plot <- ggplot(data = subset(aq_merged, id == "4781_mai_green" & Qin < 602), 
                           aes(x = Qin, y = A)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 300, color = "red", linewidth = 2, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
  scale_y_continuous(limits = c(0.6, 1.2), breaks = seq(0.6, 1.2, 0.15)) +
  labs(x = expression(bold("PAR ("*mu*"mol m"^"-2"*" s"^"-1"*")")), 
       y = expression(bold("A"["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  theme_classic(base_size = 18)

png("../drafts/figs/TT25_AQ_response.png", height = 16, width = 12, units = "in",
    res = 600)
ggarrange(five105_mai_plot, six483_mai_plot, TT24_101_mai_plot,
          fifty600_mai_plot, four781_mai_plot, ncol = 2, nrow = 3)
dev.off()


