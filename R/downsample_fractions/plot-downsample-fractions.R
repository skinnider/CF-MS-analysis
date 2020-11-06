setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# complexes
dat1 = readRDS("data/analysis/downsample_fractions/complexes.rds") %>%
  filter(quant_mode == 'iBAQ') %>%
  # summarise over ten iterations
  group_by(accession, experiment, analysis, metric, transform, missing, 
           n_fractions, quant_mode) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), mean) %>%
  ungroup()

# filter to experiments with at least 50 fractions
cdf = readRDS("data/QC/protein-groups-CDF.rds")
n_fractions = cdf %>%
  group_by(accession, experiment, quant_mode, search) %>%
  summarise(fraction_count = max(n_fractions)) %>%
  ungroup()
expts = filter(n_fractions, fraction_count >= 50)

# calculate global mean over # of fractions
## complexes
mean1 = dat1 %>%
  inner_join(expts, by = c('accession', 'experiment', 'quant_mode')) %>%
  group_by(n_fractions, quant_mode, metric) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), list(mean = mean, sd = sd)) %>%
  ungroup()

# plot
pal = colours.cafe322[c(1, 3)]
p1a = mean1 %>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  mutate(color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_line(aes(x = n_fractions, y = auroc_mean, color = color)) +
  geom_point(aes(x = n_fractions, y = auroc_mean, color = color),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', limits = c(0.58, 0.67),
                     breaks = seq(0.58, 0.68, 0.02)) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('AUC', 'Change in AUC')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p1a
p1b = mean1 %>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  group_by(quant_mode, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc_mean - lag(auroc_mean)) %>%
  ungroup() %>%
  ggplot() + 
  geom_line(aes(x = n_fractions, y = delta_auroc, color = quant_mode)) +
  geom_point(aes(x = n_fractions, y = delta_auroc, color = quant_mode),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', position = 'right') +
  scale_color_manual('', values = pal[-1]) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1)
p1b
# superimpose in cowplot
aligned_plots = align_plots(p1a, p1b, align = "hv", axis = "tblr")
p1 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p1
ggsave("fig/analysis/downsample_fractions/AUC-deltaAUC-CORUM-pearson.pdf",
       p1, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

# repeat for MI
p2a = mean1 %>%
  filter(metric == 'MI') %>%
  filter(n_fractions <= 50) %>%
  mutate(color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_line(aes(x = n_fractions, y = auroc_mean, color = color)) +
  geom_point(aes(x = n_fractions, y = auroc_mean, color = color),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', limits = c(0.58, 0.67),
                     breaks = seq(0.58, 0.68, 0.02)) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('AUC', 'Change in AUC')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p2b = mean1 %>%
  filter(metric == 'MI') %>%
  filter(n_fractions <= 50) %>%
  group_by(quant_mode, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc_mean - lag(auroc_mean)) %>%
  ungroup() %>%
  ggplot() + 
  geom_line(aes(x = n_fractions, y = delta_auroc, color = quant_mode)) +
  geom_point(aes(x = n_fractions, y = delta_auroc, color = quant_mode),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', position = 'right',
                     limits = c(0, 0.03)) +
  scale_color_manual('', values = pal[-1]) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1)
aligned_plots = align_plots(p2a, p2b, align = "hv", axis = "tblr")
p2 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p2
ggsave("fig/analysis/downsample_fractions/AUC-deltaAUC-CORUM-MI.pdf",
       p2, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

# GO
dat2 = readRDS("data/analysis/downsample_fractions/GO.rds") %>% 
  filter(quant_mode == 'iBAQ') %>%
  # filter by breadth, and get median AUROC across GO terms
  filter(between(n_chromatograms, 10, 100)) %>%
  group_by_at(vars(-go_term, -n_proteins, -n_chromatograms, -auroc)) %>%
  summarise(auroc = median(auroc)) %>%
  ungroup() %>%
  # summarise over ten iterations
  group_by_at(vars(-iteration, -n_obs, -n_pairs, -auroc)) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), mean) %>%
  ungroup()

# finally, filter to expts with at least 50 fractions, and calculate global mean
mean2 = dat2 %>%
  inner_join(expts, by = c('accession', 'experiment', 'quant_mode')) %>%
  group_by(n_fractions, quant_mode, metric) %>%
  summarise_at(vars(auroc, n_obs, n_pairs), list(mean = mean, sd = sd)) %>%
  ungroup()

# plot
pal = colours.cafe322[c(1, 3)]
p3a = mean2 %>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  mutate(color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_line(aes(x = n_fractions, y = auroc_mean, color = color)) +
  geom_point(aes(x = n_fractions, y = auroc_mean, color = color),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', limits = c(0.52, 0.545),
                     breaks = seq(0.52, 0.54, 0.01)) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('AUC', 'Change in AUC')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p3a
p3b = mean2 %>%
  filter(metric == 'pearson') %>%
  filter(n_fractions <= 50) %>%
  group_by(quant_mode, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc_mean - lag(auroc_mean)) %>%
  ungroup() %>%
  ggplot() + 
  geom_line(aes(x = n_fractions, y = delta_auroc, color = quant_mode)) +
  geom_point(aes(x = n_fractions, y = delta_auroc, color = quant_mode),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', position = 'right',
                     breaks = seq(0, 0.01, 0.005),
                     limits = c(NA, 0.011)) +
  scale_color_manual('', values = pal[-1]) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1)
p3b
# superimpose in cowplot
aligned_plots = align_plots(p3a, p3b, align = "hv", axis = "tblr")
p3 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p3
ggsave("fig/analysis/downsample_fractions/AUC-deltaAUC-GO-pearson.pdf",
       p3, width = 5.5, height = 4.75, units = "cm", useDingbats = F)

# repeat for MI
p4a = mean2 %>%
  filter(metric == 'MI') %>%
  filter(n_fractions <= 50) %>%
  mutate(color = factor(quant_mode, levels = c('iBAQ', '2'))) %>%
  ggplot() + 
  geom_rect(data = data.frame(),
            aes(xmin = 37.5, xmax = 42.5, ymin = -Inf, ymax = +Inf), 
            fill = 'grey93', color = NA, alpha = 0.5) +
  geom_line(aes(x = n_fractions, y = auroc_mean, color = color)) +
  geom_point(aes(x = n_fractions, y = auroc_mean, color = color),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('AUC', limits = c(0.51, 0.54)) +
  scale_color_manual('', values = pal, drop = FALSE,
                     labels = c('AUC', 'Change in AUC')) +
  boxed_theme() +
  theme(legend.position = 'top',
        aspect.ratio = 1)
p4a
p4b = mean2 %>%
  filter(quant_mode == 'iBAQ', metric == 'MI') %>%
  filter(n_fractions <= 50) %>%
  group_by(quant_mode, metric) %>%
  arrange(n_fractions) %>%
  mutate(delta_auroc = auroc_mean - lag(auroc_mean)) %>%
  ungroup() %>%
  ggplot() + 
  geom_line(aes(x = n_fractions, y = delta_auroc, color = quant_mode)) +
  geom_point(aes(x = n_fractions, y = delta_auroc, color = quant_mode),
             size = 0.8) +
  scale_x_continuous('# of fractions', breaks = seq(10, 50, 10)) +
  scale_y_continuous('Change in AUC', position = 'right',
                     breaks = seq(0, 0.01, 0.005)) +
  scale_color_manual('', values = pal[-1]) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1)
aligned_plots = align_plots(p4a, p4b, align = "hv", axis = "tblr")
p4 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
p4
ggsave("fig/analysis/downsample_fractions/AUC-deltaAUC-GO-MI.pdf",
       p4, width = 5.5, height = 4.75, units = "cm", useDingbats = F)
