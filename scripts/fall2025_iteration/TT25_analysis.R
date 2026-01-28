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
library(ggpubr)

# Read photosynthetic dataset
## df_photo <- read.csv("../../data/2025_2026/TT25_tri_photo_traits.csv") %>%
##   mutate(jmax_vcmax = Vcmax / Jmax)
## head(df)

# Read harvest dataset
## df_harvest <- read.csv("../../data/2025_2026/TT25_tri_harvest.csv")

# Merge photosynthetic and harvest dataset
## df <- df_photo %>%
##   full_join(df_harvest, by = "id") %>%
##   write.csv("../../data/2025_2026/TT25_tri_compiled.csv", row.names = F)

# Read compiled dataset
df <- read.csv("../../data/2025_2026/TT25_tri_compiled.csv")

# Remove sterile plant treatment from analyses
df_noSterile <- df %>%
  filter(ExpFungSource != "Wsterile" & ExpFungSource != "NWsterile") %>%
  mutate(FullExpTrt = gsub("-", "_", FullExpTrt))

# Some plot aesthetics
gm.colors <- c("#F1B700", "#00B2BE")
facet_lab <- c("Plant history: ambient", "Plant history: weeded")
names(facet_lab) <- c("NW", "W")

###############
# Anet
###############
df_noSterile$anet[c(75)] <- NA

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

# Plot prep (full)
anet_tri_results <- cld(emmeans(anet_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS, reversed = TRUE, alpha = 0.1) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

# Plot prep (ExpSoilSource * ExpFungSource interaction)
anet_soilFung_int <- cld(emmeans(anet_tri, ~ExpSoilSource*ExpFungSource, 
                                 type = "response"), 
                                 Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../../drafts/figs/TT25_tri_anet_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = anet_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 3, 
                group = full_trt, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("A"["net"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_anet_soilAMFint_plot <-
  ggplot(data = anet_soilFung_int,
         x = ExpSoilSource, 
         y = anet) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 2, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("A"["net"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

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

# Plot prep (full)
gsw_tri_results <- cld(emmeans(gsw_tri, 
                               ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                               type = "response"), 
                       Letters = LETTERS, reversed = TRUE) %>%
  data.frame() %>% 
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
gsw_soilFung_int <- cld(emmeans(gsw_tri, ~ExpSoilSource*ExpFungSource, 
                                 type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../../drafts/figs/TT25_tri_gsw_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = gsw_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, 
                width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, 
                y = 0.03, 
                group = full_trt, 
                label = .group), 
            size = 6, position = position_dodge(width = 0.75), 
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("g"["sw"]* " (mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_gsw_soilAMFint_plot <-
  ggplot(data = gsw_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 0.02, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.005)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("g"["sw"]* " (mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_gsw_soilAMFint_plot

###############
# Vcmax25
###############
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
## Weeded fungal source has greater Vcmax

cld(emmeans(vcmax_tri, pairwise~ExpSoilSource*ExpFungSource))
## Exp fungal source response only observed in weeded soils

cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpSoilSource))
## Plants with weeded treatment history have greater Vcmax
## in non-weeded treatment soils

# Plot prep (full)
vcmax_tri_results <- cld(emmeans(vcmax_tri, 
                                 ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                 type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>%
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
vcmax_soilFung_int <- cld(emmeans(vcmax_tri, ~ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

# Plot prep (ExpFungSource)
vcmax_soilFung_ind <- cld(emmeans(vcmax_tri, ~ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_vcmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = vcmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 30, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_vcmax_soilAMFint_plot <-
  ggplot(data = vcmax_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, 
                    y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 20, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_vcmax_soilAMFint_plot

## Individual AMF response
vcmax_AMF_ind_plot <- ggplot(data = vcmax_soilFung_ind) +
  geom_errorbar(aes(x = ExpFungSource,
                    y = response,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5) +
  geom_point(aes(x = ExpFungSource,
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21) +
  geom_text(aes(x = ExpFungSource, y = 20, label = .group), 
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Fungal Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

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

# Plot prep (full)
jmax_tri_results <- cld(emmeans(jmax_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>%
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
jmax_soilFung_int <- cld(emmeans(jmax_tri, ~ExpSoilSource*ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

# Plot prep (ExpFungSource)
jmax_soilFung_ind <- cld(emmeans(jmax_tri, ~ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


# Full plot
png("../../drafts/figs/TT25_tri_jmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = jmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 60, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_jmax_soilAMFint_plot <-
  ggplot(data = jmax_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, 
                    y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 40, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## AM fungal source (individual) plot
jmax_AMF_ind_plot <- ggplot(data = jmax_soilFung_ind) +
  geom_errorbar(aes(x = ExpFungSource,
                    y = response,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5) +
  geom_point(aes(x = ExpFungSource,
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21) +
  geom_text(aes(x = ExpFungSource, y = 40, label = .group), 
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Fungal Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

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

# Plot prep (full)
jmaxvcmax_tri_results <- cld(emmeans(jmaxvcmax_tri, 
                                     ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                             Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_jmaxvcmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = jmaxvcmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 0.6, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, 0.1)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* ": V"["cmax25"]*" (unitless)")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

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


emmeans(rd_tri, pairwise~ExpSoilSource*ExpFungSource, type = "response")
## Plants with weeded fungal source have greater Rd25 than
## plants with ambient fungal source

cld(emmeans(rd_tri, pairwise~ExpFungSource*PlantGMTrt, type = "response"))
## Weeded fungal source response is driven by plants that have historically
## grown in the weeded treatment

# Plot prep (full)
rd_tri_results <- cld(emmeans(rd_tri, 
                              ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                      Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_rd_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = rd_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 3, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("Rd"["25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

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

# Plot prep (full)
iwue_tri_results <- cld(emmeans(iwue_tri, 
                              ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                      Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_iwue_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = iwue_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = emmean,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = emmean,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 150, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("iWUE ("*mu*"mol mol"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

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

# Plot prep (full)
tla_tri_results <- cld(emmeans(tla_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                        Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_tla_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = tla_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = emmean,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = emmean,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 80, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(20, 80), breaks = seq(20, 80, 20)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# total biomass
###############
tbio_tri <- lm(log(total_biomass) ~ PlantGMTrt * ExpSoilSource * ExpFungSource, 
               data = df_noSterile)

# Check normality assumptions
plot(tbio_tri)
qqnorm(residuals(tbio_tri))
qqline(residuals(tbio_tri))
hist(residuals(tbio_tri))
shapiro.test(residuals(tbio_tri))
outlierTest(tbio_tri)

# Model output
summary(tbio_tri)
Anova(tbio_tri)

###############
# Root:shoot
###############
df_noSterile$root_shoot[c(25, 107, 109, 111, 113, 114, 115)] <- NA

rootshoot_tri <- lm(log(root_shoot) ~ PlantGMTrt * ExpSoilSource * ExpFungSource, 
               data = df_noSterile)

# Check normality assumptions
plot(rootshoot_tri)
qqnorm(residuals(rootshoot_tri))
qqline(residuals(rootshoot_tri))
hist(residuals(rootshoot_tri))
shapiro.test(residuals(rootshoot_tri))
outlierTest(rootshoot_tri)

# Model output
summary(rootshoot_tri)
Anova(rootshoot_tri)

###############
# Rhizome mass
###############
rhizome_tri <- lm(log(rhizome_mass_g) ~ PlantGMTrt * ExpSoilSource * ExpFungSource, 
                    data = df_noSterile)

# Check normality assumptions
plot(rhizome_tri)
qqnorm(residuals(rhizome_tri))
qqline(residuals(rhizome_tri))
hist(residuals(rhizome_tri))
shapiro.test(residuals(rhizome_tri))
outlierTest(rhizome_tri)

# Model output
summary(rhizome_tri)
Anova(rhizome_tri)

###############
# Plot arrangements
###############
png("../../drafts/figs/TT25_soilAMint_plot.png", height = 8, 
    width = 8, units = "in", res = 600)
ggarrange(tri_anet_soilAMFint_plot, tri_gsw_soilAMFint_plot,
          tri_vcmax_soilAMFint_plot, tri_jmax_soilAMFint_plot,
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)"), align = "hv")
dev.off()

png("../../drafts/figs/TT25_AMind_plot.png", height = 4, 
    width = 10, units = "in", res = 600)
ggarrange(vcmax_AMF_ind_plot, jmax_AMF_ind_plot,
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)"), align = "hv",
          font.label = list(size = 18))
dev.off()


