setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# complexes
dat1 = readRDS("data/analysis/min_pairs/complexes.rds") %>%
  filter(quant_mode == 'iBAQ')
# filter to those where min_fractions is not NA
dat1 %<>% filter(!is.na(min_fractions))

# calculate mean AUC per min_fractions/min_pairs
means1 = dat1 %>%
  group_by(min_pairs, min_fractions) %>%
  summarise(n = n(), auroc = mean(auroc, na.rm = T)) %>%
  ungroup()

# add in mean AUC from min_fractions experiment, where min_pairs == 0
dat2 = readRDS("data/analysis/min_fractions/complexes.rds") %>%
  filter(quant_mode == 'iBAQ')
  
# calculate mean over # fractions
means2 = dat2 %>%
  mutate(min_pairs = 0) %>%
  dplyr::rename(min_fractions = n_fractions) %>%
  group_by(min_pairs, min_fractions) %>%
  summarise(n = n(), auroc = mean(auroc, na.rm = T)) %>%
  ungroup()

# plot
p1 = bind_rows(means1, means2) %>%
  filter(min_fractions <= 10) %>%
  filter(min_fractions >= min_pairs) %>%
  ggplot(aes(x = factor(min_fractions), y = factor(min_pairs), fill = auroc)) +
  geom_tile(color = 'white') +
  scale_fill_paletteer_c("pals::coolwarm", name = 'AUC ', 
                         breaks = seq(0.665, 0.675, 0.01)) +
  scale_x_discrete('Min. # of fractions', expand = c(0, 0)) +
  scale_y_discrete('Min. # of pairs', expand = c(0, 0)) +
  guides(fill = guide_colorbar(ticks = F, frame.colour = 'black')) +
  coord_fixed() +
  boxed_theme() +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.25, 'lines'),
        legend.key.width = unit(0.25, 'lines'))
p1
ggsave("fig/analysis/min_pairs/min_pairs-heatmap.pdf", p1,
       width = 4.6, height = 4.6, units = "cm", useDingbats = F)
