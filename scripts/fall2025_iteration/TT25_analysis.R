##############################################################################
# Prep
##############################################################################

# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)

# Read in compiled dataset
df <- read.csv("../../data/2025_2026/TT25_tri_photo_traits.csv") %>%
  mutate(jmax_vcmax = Vcmax / Jmax)
head(df)

# Remove sterile plant treatment from analyses
df_noSterile <- df %>%
  filter(ExpFungSource != "Wsterile" & ExpFungSource != "NWsterile")


df_noSterile %>%
  filter(!is.na(anet)) %>%
  group_by(PlantGMTrt, ExpSoilSource, ExpFungSource) %>%
  summarize(length(anet))

###############
# Anet
###############
df_noSterile$anet[75] <- NA

anet_tri <- lmer(log(anet + 1) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | machine), data = df_noSterile)

# Check normality assumptions
plot(anet_tri)
qqnorm(residuals(anet_tri))
qqline(residuals(anet_tri))
hist(residuals(anet_tri))
shapiro.test(residuals(anet_tri))
outlierTest(anet_tri)

# Model output
summary(anet_tri)
Anova(anet_tri)
r.squaredGLMM(anet_tri)

# Pairwise comparisons
cld(emmeans(anet_tri, pairwise~PlantGMTrt*ExpFungSource, type = "response"), alpha = 0.1)
## Plants inoculated with AM fungi from weeded treatment have greater net
## photosynthesis than plants inoculated with AM fungi from non-weeded treatment
## However, this response is only observed in plants that historically grew in
## non-weeded treatment

cld(emmeans(anet_tri, pairwise~ExpSoilSource*ExpFungSource), alpha = 0.1)

###############
# gsw
###############
gsw_tri <- lmer(log(gsw) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | machine), data = filter(df_noSterile, gsw > 0))

# Check normality assumptions
plot(gsw_tri)
 qqnorm(residuals(gsw_tri))
qqline(residuals(gsw_tri))
hist(residuals(gsw_tri))
shapiro.test(residuals(gsw_tri))
outlierTest(gsw_tri)

# Model output
summary(gsw_tri)
Anova(gsw_tri)
r.squaredGLMM(gsw_tri)

# Pairwise comparisons
cld(emmeans(gsw_tri, pairwise~ExpSoilSource*ExpFungSource, type = "response"), alpha = 0.1)
## Plants growing in non-weeded soil have greater stomatal conductance than
## plants growing in weeded soil; however, this response is only observed
## when plants are inoculated with nonweeded AMF community


###############
# Vcmax25
###############
df_noSterile$Vcmax[11] <- NA

vcmax_tri <- lmer(log(Vcmax) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    (1 | machine), data = df_noSterile)

# Check normality assumptions
plot(vcmax_tri)
qqnorm(residuals(vcmax_tri))
qqline(residuals(vcmax_tri))
hist(residuals(vcmax_tri))
shapiro.test(residuals(vcmax_tri))
outlierTest(vcmax_tri)

# Model output
summary(vcmax_tri)
Anova(vcmax_tri)
r.squaredGLMM(vcmax_tri)

# Pairwise comparisons
emmeans(vcmax_tri, pairwise~ExpFungSource)
## Weeded fungal source has marginally greater Vcmax

cld(emmeans(vcmax_tri, pairwise~ExpSoilSource*ExpFungSource))
## Exp fungal source response only observed in weeded soils

cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpSoilSource))
## Plants with weeded treatment history have greater Vcmax
## in non-weeded treatment soils

###############
# Jmax
###############
jmax_tri <- lmer(log(Jmax) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    (1 | machine), data = df_noSterile)

# Check normality assumptions
plot(jmax_tri)
qqnorm(residuals(jmax_tri))
qqline(residuals(jmax_tri))
hist(residuals(jmax_tri))
shapiro.test(residuals(jmax_tri))
outlierTest(jmax_tri)

# Model output
summary(jmax_tri)
Anova(jmax_tri)
r.squaredGLMM(jmax_tri)

# Pairwise comparisons
emmeans(jmax_tri, pairwise~ExpFungSource)
## Weeded fungal source has marginally greater Jmax than non-weeded fungal source

cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpFungSource))
## Weeded fungal source response is only observed in plants with
## non-weeded treatment history

cld(emmeans(jmax_tri, pairwise~ExpSoilSource*ExpFungSource))
## Weeded fungal source response is only observed in weeded soils

cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpSoilSource), alpha = 0.1)
## Non-weeded soil source marginally increases Jmax compared to weeded soil 
## source, but only in plants with weeded treatment history

###############
# Jmax:Vcmax
###############
jmaxvcmax_tri <- lmer(log(jmax_vcmax) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | machine), data = df_noSterile)

# Check normality assumptions
plot(jmaxvcmax_tri)
qqnorm(residuals(jmaxvcmax_tri))
qqline(residuals(jmaxvcmax_tri))
hist(residuals(jmaxvcmax_tri))
shapiro.test(residuals(jmaxvcmax_tri))
outlierTest(jmaxvcmax_tri)

# Model output
summary(jmaxvcmax_tri)
Anova(jmaxvcmax_tri)
r.squaredGLMM(jmaxvcmax_tri)

# Pairwise comparisons
cld(emmeans(jmaxvcmax_tri, pairwise~PlantGMTrt*ExpSoilSource, type = "response"))
## Weeded plants have greater Jmax:Vcmax (approaching 1), but this 
## response is only observed in non-weeded soils

###############
# Rd
###############
rd_tri <- lmer(log(Rd) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                 (1 | machine), data = df_noSterile)

# Check normality assumptions
plot(rd_tri)
qqnorm(residuals(rd_tri))
qqline(residuals(rd_tri))
hist(residuals(rd_tri))
shapiro.test(residuals(rd_tri))
outlierTest(rd_tri)

# Model output
summary(rd_tri)
Anova(rd_tri)
r.squaredGLMM(rd_tri)

# Pairwise comparisons
emmeans(rd_tri, pairwise~ExpFungSource, type = "response")
## 

###############
# iWUE
###############
iwue_tri <- lmer(iwue ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                        (1 | machine), data = filter(df_noSterile, iwue > 0 & iwue < 200))

# Check normality assumptions
plot(iwue_tri)
qqnorm(residuals(iwue_tri))
qqline(residuals(iwue_tri))
hist(residuals(iwue_tri))
shapiro.test(residuals(iwue_tri))
outlierTest(iwue_tri)

# Model output
summary(iwue_tri)
Anova(iwue_tri) # No effect of any treatment combination
r.squaredGLMM(iwue_tri)


###############
# TLA
###############
tla_tri <- lm(tla ~ PlantGMTrt * ExpSoilSource * ExpFungSource, 
              data = df_noSterile)

# Check normality assumptions
plot(tla_tri)
qqnorm(residuals(tla_tri))
qqline(residuals(tla_tri))
hist(residuals(tla_tri))
shapiro.test(residuals(tla_tri))
outlierTest(tla_tri)

# Model output
summary(tla_tri)
Anova(tla_tri) # No treatment effect
r.squaredGLMM(tla_tri)




