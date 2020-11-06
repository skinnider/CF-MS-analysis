# Plot the results of the main analysis grid shown in Figure 2.
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# 1. Complexes ####
dat1 = readRDS("data/analysis/analysis_grid/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))

# get the best analysis combination per ...
## metric
bests1 = dat1 %>%
  group_by(metric, transform, missing) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup() %>%
  group_by(metric) %>%
  filter(median == max(median)) %>%
  ungroup()
## missing values
bests2 = dat1 %>%
  group_by(missing, transform, metric) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup() %>%
  group_by(missing) %>%
  filter(median == max(median)) %>%
  ungroup()
## transformation
bests3 = dat1 %>%
  group_by(transform, missing, metric) %>%
  summarise(mean = mean(auroc), median = median(auroc), n = n()) %>%
  ungroup() %>%
  group_by(transform) %>%
  filter(median == max(median)) %>%
  ungroup()

# colour scheme
pal = brewer.pal(8, 'Set2')
col = pal[3]

# plot
## calculate medians
medians1a = dat1 %>%
  mutate(metric = clean_metric(metric)) %>%
  group_by(metric) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p1a = dat1 %>%
  # clean up metrics
  mutate(metric = clean_metric(metric),
         title = 'Measure of association') %>%
  # plot
  ggplot(aes(y = reorder(metric, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ title) +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians1a, aes(x = 0.345, y = metric, label = label), 
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.3, 0.9, 0.1)) +
  coord_cartesian(xlim = c(0.35, 0.72)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1a
ggsave("fig/analysis/analysis_grid/complex-AUC-metrics-CORUM.pdf", p1a,
       width = 7, height = 8, units = "cm", useDingbats = F)

## impact of missing values
medians1b = dat1 %>%
  group_by(missing) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p1b = dat1 %>%
  ggplot(aes(y = reorder(missing, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ 'Missing values') +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians1b, aes(x = 0.345, y = missing, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.3, 0.9, 0.1)) +
  coord_cartesian(xlim = c(0.35, 0.72)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1b
ggsave("fig/analysis/analysis_grid/complex-AUC-missing-CORUM.pdf", p1b,
       width = 5, height = 2.6, units = "cm", useDingbats = F)

## impact of log-transformation
medians1c = dat1 %>%
  group_by(transform) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p1c = dat1 %>%
  ggplot(aes(y = reorder(transform, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ 'Normalization') +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians1c, aes(x = 0.345, y = transform, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.3, 0.9, 0.1)) +
  coord_cartesian(xlim = c(0.35, 0.72)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1c
ggsave("fig/analysis/analysis_grid/complex-AUC-normalize-CORUM.pdf", p1c,
       width = 5, height = 2.6, units = "cm", useDingbats = F)

# combine
p1 = (p1a + theme(axis.title.x = element_blank(),
                  plot.margin = margin(rep(0, 4)))) / 
  (p1b + theme(axis.title.x = element_blank(),
               plot.margin = margin(rep(0, 4)))) / 
  (p1c + theme(plot.margin = margin(rep(0, 4)))) +
  plot_layout(heights = c(24, 3, 3))
p1
ggsave("fig/analysis/analysis_grid/complex-AUCs-CORUM.pdf", p1, 
       width = 7.5, height = 11.5, units = "cm", useDingbats = F)

# 2. GO ####
dat2 = readRDS("data/analysis/analysis_grid/GO.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # recode transformation
  mutate(transform = fct_recode(transform, 
                                'log-transform' = 'log',
                                'quantile normalization' = 'quantile'))
# for each combination, filter by breadth, and get median AUROC across GO terms
dat2 %<>%
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup()

# colour scheme
pal = brewer.pal(8, 'Set2')
col = pal[7]

# metric
medians2a = dat2 %>%
  mutate(metric = clean_metric(metric)) %>%
  group_by(metric) %>%
  summarise(median = median(auroc),
            q1 = quantile(auroc, probs = 0.25),
            q3 = quantile(auroc, probs = 0.75),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p2a = dat2 %>%
  mutate(metric = clean_metric(metric),
         title = 'Measure of association') %>%
  ggplot(aes(y = reorder(metric, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ title) +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.6, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2a, aes(x = 0.448, y = metric, label = label), 
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.48, 0.9, 0.02)) +
  coord_cartesian(xlim = c(0.45, 0.565)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2a
ggsave("fig/analysis/analysis_grid/GO-AUC-metrics.pdf", p2a,
       width = 7, height = 8, units = "cm", useDingbats = F)

# missing values
medians2b = dat2 %>%
  group_by(missing) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p2b = means %>%
  ggplot(aes(y = reorder(missing, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ 'Missing values') +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2b, aes(x = 0.448, y = missing, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.48, 0.9, 0.02)) +
  coord_cartesian(xlim = c(0.45, 0.565)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2b
ggsave("fig/analysis/analysis_grid/GO-AUC-missing.pdf", p2b,
       width = 5, height = 2.6, units = "cm", useDingbats = F)

# normalization
medians2c = dat2 %>%
  group_by(transform) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
p2c = dat2 %>%
  ggplot(aes(y = reorder(transform, auroc, stats::median),
             x = auroc, fill = '1', color = '1')) +
  facet_grid(~ 'Normalization') +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(width = 0.6, alpha = 0.5, outlier.shape = NA, coef = 0) + 
  geom_text(data = medians2c, aes(x = 0.448, y = transform, label = label),
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.48, 0.9, 0.02)) +
  coord_cartesian(xlim = c(0.45, 0.565)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2c
ggsave("fig/analysis/analysis_grid/GO-AUC-normalize.pdf", p2c,
       width = 5, height = 2.25, units = "cm", useDingbats = F)

# combine
p2 = (p2a + theme(axis.title.x = element_blank(),
                  plot.margin = margin(rep(0, 4)))) / 
  (p2b + theme(axis.title.x = element_blank(),
               plot.margin = margin(rep(0, 4)))) / 
  (p2c + theme(plot.margin = margin(rep(0, 4)))) +
  plot_layout(heights = c(24, 3, 3))
p2
ggsave("fig/analysis/analysis_grid/GO-AUCs.pdf", p2, 
       width = 7, height = 11.5, units = "cm", useDingbats = F)

# combine complexes (CORUM) and GO
p3 = p1 | p2
ggsave("fig/analysis/analysis_grid/AUCs.pdf", p3, 
       width = 12.8, height = 12, units = "cm", useDingbats = F)
