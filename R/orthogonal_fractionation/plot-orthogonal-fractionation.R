setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read AUC
auc = readRDS("data/analysis/orthogonal_fractionation/AUC.rds") %>% 
  filter(complex_set == 'held-out')

# plot everything
pal = brewer.pal(7, 'Blues')[-1]
p1 = auc %>%
  mutate(pct_sec = n_sec / n_datasets,
         n_datasets = ifelse(n_datasets == 2, 2, 
                             paste(n_datasets, 'datasets'))) %>%
  ggplot(aes(x = factor(n_sec), y = auc, fill = pct_sec,
             color = pct_sec)) +
  facet_grid(~ n_datasets, scales = 'free', space = 'free') +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
               outlier.size = 0.3, size = 0.4) + 
  scale_fill_distiller(palette = 'Blues', direction = 1,
                       name = '% SEC datasets ', limits = c(-0.2, 1),
                       breaks = c(-0.2, 1), labels = c(0, 100)) +
  scale_color_distiller(palette = 'Blues', direction = 1,
                        name = '% SEC datasets ', limits = c(-0.2, 1),
                        breaks = c(-0.2, 1), labels = c(0, 100)) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  scale_y_continuous('AUC') +
  scale_x_discrete('# of SEC datasets') +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.25, 'lines'))
p1
ggsave("fig/analysis/orthogonal_fractionation/pct-sec.pdf", p1,
       width = 14.5, height = 4.5, units = "cm", useDingbats = FALSE)

# plot only four or eight for the main text
p2 = auc %>%
  filter(n_datasets %in% c(4, 8)) %>%
  mutate(pct_sec = n_sec / n_datasets,
         n_datasets = ifelse(n_datasets == 2, 2, 
                             paste(n_datasets, 'datasets'))) %>%
  ggplot(aes(x = factor(n_sec), y = auc, fill = pct_sec,
             color = pct_sec)) +
  facet_grid(~ n_datasets, scales = 'free', space = 'free') +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
               outlier.size = 0.3, size = 0.4) + 
  scale_fill_distiller(palette = 'Blues', direction = 1,
                       name = '% SEC', limits = c(-0.2, 1),
                       breaks = c(-0.2, 1), labels = c(0, 100)) +
  scale_color_distiller(palette = 'Blues', direction = 1,
                        name = '% SEC', limits = c(-0.2, 1),
                        breaks = c(-0.2, 1), labels = c(0, 100)) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  scale_y_continuous('AUC') +
  scale_x_discrete('# of SEC datasets') +
  boxed_theme() +
  theme(legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.25, 'lines'))
p2
ggsave("fig/analysis/orthogonal_fractionation/pct-sec-2-4.pdf", p2,
       width = 7, height = 4.64, units = "cm", useDingbats = FALSE)
