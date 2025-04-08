# TT25_phys_analyses.R
# R script that analyzes cleaned dataset for 2025 TT plant-soil
# feedback experiment
# 
# Note: all paths assume that the folder containing this R script
# is the working root directory

#####################################################################
# Libraries and data file
#####################################################################
# Libraries
library(tidyverse)
library(car)
library(emmeans)
library(multcomp)

# Read data file
df <- read.csv("../data/TT25_full_physiology.csv")

# Remove sterile treatments from analyses (too little power)
df_model <- filter(df, ExpFungSource != "NWsterile" & ExpFungSource != "Wsterile")

## Add code for facet labels
facet.labs <- c("Treatment legacy: non-weeded", 
                "Treatment legacy: weeded")
names(facet.labs) <- c("NW", "W")

## Color palettes
gm.colors <- c("#F1B700", "#00B2BE")

#####################################################################
#####################################################################
# Fit models!
#####################################################################
#####################################################################

#####################################################################
# Anet (without sterile treatments)
#####################################################################
# Fit model
anet_lm <- lm(anet ~ plantGMtrt * expSoilSource * ExpFungSource,
              data = subset(df_model, anet > 0))

# Check model assumptions
qqnorm(residuals(anet_lm))
qqline(residuals(anet_lm))
hist(residuals(anet_lm))
densityPlot(residuals(anet_lm))
shapiro.test(residuals(anet_lm))
outlierTest(anet_lm)

# Model results
summary(anet_lm)
Anova(anet_lm)

# Post hoc comparisons: none, no significant effect of any trt combo

#####################################################################
# gsw (without sterile treatments)
#####################################################################
# Remove a strong outlier (gsw = 0.002)
df_model$gsw[34] <- NA

# Fit model
gsw_lm <- lm(log(gsw) ~ plantGMtrt * expSoilSource * ExpFungSource,
              data = df_model)

# Check model assumptions
qqnorm(residuals(gsw_lm))
qqline(residuals(gsw_lm))
hist(residuals(gsw_lm))
densityPlot(residuals(gsw_lm))
shapiro.test(residuals(gsw_lm))
outlierTest(gsw_lm)

# Model results
summary(gsw_lm)
Anova(gsw_lm)

# Post hoc comparisons

## Main plant GM trt effect
emmeans(gsw_lm, pairwise~plantGMtrt, type = "response") 
### nonweeded treatment (GM presence) legacy = reduced stomatal conductance

## Main fungal source effect
emmeans(gsw_lm, pairwise~ExpFungSource, type = "response")
### nonweeded AM fungal community = increased stomatal conductance

## Weak interaction between plant GM trt and fungal source
cld(emmeans(gsw_lm, pairwise~plantGMtrt * ExpFungSource, type = "response"))
### increased stomatal conductance from nonweeded AM fungal community
### only observed in plants with nonweeded treatment legacy

# Plot prep
gsw_plot_prep <- cld(
  emmeans(gsw_lm, ~plantGMtrt * expSoilSource*ExpFungSource, type = "response"), 
  Letters = LETTERS) %>%
  mutate(full_trt = str_c("plant", plantGMtrt, "_soil", expSoilSource, "_fung", ExpFungSource),
         full_trt = factor(full_trt, levels = c("plantNW_soilNW_fungNW",
                                                "plantNW_soilW_fungNW",
                                                "plantNW_soilNW_fungW",
                                                "plantNW_soilW_fungW",
                                                "plantW_soilW_fungW",
                                                "plantW_soilNW_fungW",
                                                "plantW_soilW_fungNW",
                                                "plantW_soilNW_fungNW")),
         plot_trt = str_c("soil", expSoilSource, "_fung", ExpFungSource),
         plot_trt = factor(plot_trt, levels = c("soilW_fungW",
                                                "soilW_fungNW",
                                                "soilNW_fungW",
                                                "soilNW_fungNW")))


gsw_plot <- ggplot() +
  geom_errorbar(data = gsw_plot_prep, 
                aes(x = plot_trt, y = response, ymin = lower.CL, ymax = upper.CL),
                width = 0.3) +
  geom_point(data = gsw_plot_prep, 
             aes(x = plot_trt, 
                 y = response, 
                 color = plantGMtrt), size = 4) +
  geom_text(data = gsw_plot_prep,
            aes(x = plot_trt, y = 0.03, label = .group), 
            fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01)) +
  scale_color_manual(values = gm.colors) +
  labs(x = NULL, 
       y = expression("Stomatal conductance (mol m"^"-2"*"s"^"-1"*")"),
       color = NULL) +
  facet_grid(~plantGMtrt, labeller = labeller(plantGMtrt = facet.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("../plots/TT25_gsw_plot.png", 
    width = 10, height = 6, units = "in", res = 600)
gsw_plot
dev.off()

#####################################################################
# iWUE (without sterile treatments)
#####################################################################
# Remove a strong outlier (gsw = 0.002)
df_model$iwue[64] <- NA

# Fit model
iwue_lm <- lm(iwue ~ plantGMtrt * expSoilSource * ExpFungSource,
             data = df_model)

# Check model assumptions
qqnorm(residuals(iwue_lm))
qqline(residuals(iwue_lm))
hist(residuals(iwue_lm))
densityPlot(residuals(iwue_lm))
shapiro.test(residuals(iwue_lm))
outlierTest(iwue_lm)

# Model results
summary(iwue_lm)
Anova(iwue_lm)

# Pairwise comparisons
emmeans(iwue_lm, pairwise~plantGMtrt)




#####################################################################
# Vcmax (without sterile treatments)
#####################################################################

# Fit model
vcmax_lm <- lm(Vcmax ~ plantGMtrt * expSoilSource * ExpFungSource,
             data = df_model)

# Check model assumptions
qqnorm(residuals(vcmax_lm))
qqline(residuals(vcmax_lm))
hist(residuals(vcmax_lm))
densityPlot(residuals(vcmax_lm))
shapiro.test(residuals(vcmax_lm))
outlierTest(vcmax_lm)

# Model results
summary(vcmax_lm)
Anova(vcmax_lm)

# Post hoc comparisons: none, no significant effect of any trt combo

#####################################################################
# Jmax (without sterile treatments)
#####################################################################

# Fit model
jmax_lm <- lm(log(Jmax) ~ plantGMtrt * expSoilSource * ExpFungSource,
               data = df_model)

# Check model assumptions
qqnorm(residuals(jmax_lm))
qqline(residuals(jmax_lm))
hist(residuals(jmax_lm))
densityPlot(residuals(jmax_lm))
shapiro.test(residuals(jmax_lm))
outlierTest(jmax_lm)

# Model results
summary(jmax_lm)
Anova(jmax_lm)

# Post hoc comparisons: none, no significant effect of any trt combo

#####################################################################
# Quantum yield (without sterile treatments)
#####################################################################
# Remove strong outliers
df_model$PhiPS2[c(7, 79)] <- NA
df_model$PhiPS2[c(72, 86)] <- NA
df_model$PhiPS2[c(80, 84)] <- NA

# Fit model
phips2_lm <- lm(PhiPS2 ~ plantGMtrt * expSoilSource * ExpFungSource,
              data = subset(df_model, PhiPS2 > 0))

# Check model assumptions
qqnorm(residuals(phips2_lm))
qqline(residuals(phips2_lm))
hist(residuals(phips2_lm))
densityPlot(residuals(phips2_lm))
shapiro.test(residuals(phips2_lm))
outlierTest(phips2_lm)

# Model results
summary(phips2_lm)
Anova(phips2_lm)

# Post hoc comparisons

## Main plant GM trt effect
emmeans(phips2_lm, pairwise~plantGMtrt, type = "response") 
### nonweeded treatment (GM presence) legacy = reduced phips2

## Weak interaction between plant GM trt and fungal source
cld(emmeans(phips2_lm, pairwise~plantGMtrt * ExpFungSource, type = "response"))
### no pairwise effects... unsure why this is pinging as an interaction

## Three-way interaction between treatment legacy, soil source, fungal source
phiPS2_plot_prep <- cld(
  emmeans(phips2_lm, ~plantGMtrt * expSoilSource*ExpFungSource, type = "response"), 
  Letters = LETTERS) %>%
  mutate(full_trt = str_c("plant", plantGMtrt, "_soil", expSoilSource, "_fung", ExpFungSource),
         full_trt = factor(full_trt, levels = c("plantNW_soilNW_fungNW",
                                                "plantNW_soilW_fungNW",
                                                "plantNW_soilNW_fungW",
                                                "plantNW_soilW_fungW",
                                                "plantW_soilW_fungW",
                                                "plantW_soilNW_fungW",
                                                "plantW_soilW_fungNW",
                                                "plantW_soilNW_fungNW")),
         plot_trt = str_c("soil", expSoilSource, "_fung", ExpFungSource),
         plot_trt = factor(plot_trt, levels = c("soilW_fungW",
                                                "soilW_fungNW",
                                                "soilNW_fungW",
                                                "soilNW_fungNW")))
### PhiPS2 is significantly lower in nonweeded plants grown in weeded soil and
### with the weeded AM fungal community than weeded plants grown in non-weeded soil
### and nonweeded AM fungal community AND  in weeded plants grown in weeded soil
### and weeded AM fungal community

phiPSII_plot <- ggplot() +
  geom_errorbar(data = phiPS2_plot_prep, 
                aes(x = plot_trt, y = emmean, ymin = lower.CL, ymax = upper.CL),
                width = 0.3) +
  geom_point(data = phiPS2_plot_prep, 
             aes(x = plot_trt, 
                 y = emmean, 
                 color = plantGMtrt), size = 4) +
  geom_text(data = phiPS2_plot_prep,
            aes(x = plot_trt, y = 0.8, label = .group), 
            fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0.4, 0.8), breaks = seq(0.4, 0.8, 0.1)) +
  scale_color_manual(values = gm.colors) +
  labs(x = NULL, 
       y = expression(phi["PSII"]),
       color = NULL) +
  facet_grid(~plantGMtrt, labeller = labeller(plantGMtrt = facet.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("../plots/TT25_phipsii_plot.png", 
    width = 10, height = 6, units = "in", res = 600)
phiPSII_plot
dev.off()

#####################################################################
# Total leaf area (without sterile treatments)
#####################################################################

# Fit model
tla_lm <- lm(sqrt(total_leaf_area) ~ plantGMtrt * expSoilSource * ExpFungSource,
              data = df_model)

# Check model assumptions
qqnorm(residuals(tla_lm))
qqline(residuals(tla_lm))
hist(residuals(tla_lm))
densityPlot(residuals(tla_lm))
shapiro.test(residuals(tla_lm))
outlierTest(tla_lm)

# Model results
summary(tla_lm)
Anova(tla_lm)

# Post hoc comparisons

## Main plant GM trt effect
emmeans(tla_lm, pairwise~plantGMtrt, type = "response") 
### nonweeded treatment (GM presence) legacy = reduced total leaf area

## Main plant GM trt effect
emmeans(tla_lm, pairwise~ExpFungSource, type = "response") 
### nonweeded treatment (GM presence) legacy = reduced total leaf area

## Three-way interaction between treatment legacy, soil source, fungal source
tla_plot_prep <- cld(
  emmeans(tla_lm, ~plantGMtrt * expSoilSource*ExpFungSource, type = "response"),
  Letters = LETTERS) %>%
  mutate(full_trt = str_c("plant", plantGMtrt, "_soil", expSoilSource, "_fung", ExpFungSource),
         full_trt = factor(full_trt, levels = c("plantNW_soilNW_fungNW",
                                                "plantNW_soilW_fungNW",
                                                "plantNW_soilNW_fungW",
                                                "plantNW_soilW_fungW",
                                                "plantW_soilW_fungW",
                                                "plantW_soilNW_fungW",
                                                "plantW_soilW_fungNW",
                                                "plantW_soilNW_fungNW")),
         plot_trt = str_c("soil", expSoilSource, "_fung", ExpFungSource),
         plot_trt = factor(plot_trt, levels = c("soilW_fungW",
                                                "soilW_fungNW",
                                                "soilNW_fungW",
                                                "soilNW_fungNW")))
### iiiinteresting. Really no effect of treatment combinations on total leaf area
### except for one. Nonweeded plants grown in nonweeded soil but inoculated with
### weeded AMF exhibited significantly reduced total leaf area compared to nonweeded
### plants grown in nonweeded soil and nonweeded AMF

tla_plot <- ggplot() +
  geom_errorbar(data = tla_plot_prep, 
                aes(x = plot_trt, y = response, ymin = lower.CL, ymax = upper.CL),
                width = 0.3) +
  geom_point(data = tla_plot_prep, 
             aes(x = plot_trt, 
                 y = response, 
                 color = plantGMtrt), size = 4) +
  geom_text(data = tla_plot_prep,
            aes(x = plot_trt, y = 150, label = .group), 
            fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_color_manual(values = gm.colors) +
  labs(x = NULL, 
       y = expression("Total leaf area (cm"^"2"*")"),
       color = NULL) +
  facet_grid(~plantGMtrt, labeller = labeller(plantGMtrt = facet.labs)) +
  guides(color = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tla_plot

png("../plots/TT25_tla_plot.png", 
    width = 10, height = 6, units = "in", res = 600)
tla_plot
dev.off()

