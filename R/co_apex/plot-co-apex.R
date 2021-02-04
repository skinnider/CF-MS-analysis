setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(ggstance)
source("R/theme.R")

# complexes
dat1 = readRDS("data/analysis/co_apex/complexes.rds") %>%
  filter(quant_mode == 'iBAQ')

pal = brewer.pal(8, 'Set2')
col = pal[3]
medians1 = dat1 %>%
  filter(quant_mode == 'iBAQ') %>%
  mutate(metric = clean_metric(metric)) %>%
  group_by(metric) %>%
  summarise(median = median(auroc),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup() %>% 
  arrange(median) %>% 
  mutate(metric = factor(metric))
rect = medians1 %>% filter(metric == 'Co-apex score') 
p1 = dat1 %>%
  filter(quant_mode == 'iBAQ') %>%
  # clean up metrics
  mutate(metric = clean_metric(metric),
         metric = reorder(metric, auroc, stats::median),
         title = 'Measure of association') %>%
  # plot
  ggplot(aes(y = metric, fill = '1', color = '1')) +
  facet_grid(~ title) +
  geom_blank() +
  geom_crossbarh(data = rect, aes(xmin = -Inf, xmax = +Inf, x = 0.5), 
                 # fill = 'grey90', fill = pal,
                 alpha = 0.2, color = NA, width = 1.05) +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(aes(x = auroc), width = 0.6, alpha = 0.5, outlier.shape = NA,
                coef = 0) + 
  geom_text(data = medians1, aes(x = 0.38, y = metric, label = label), 
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.5, 0.75, 0.1)) +
  coord_cartesian(xlim = c(0.39, 0.77)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p1
ggsave("fig/analysis/co_apex/complex-AUC-metrics-CORUM.pdf", p1,
       width = 6.3, height = 7.05, units = "cm", useDingbats = F)

# GO
dat2 = readRDS("data/analysis/co_apex/GO.rds") %>%
  filter(quant_mode == 'iBAQ') 
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
medians2 = dat2 %>%
  mutate(metric = clean_metric(metric)) %>%
  group_by(metric) %>%
  summarise(median = median(auroc),
            q1 = quantile(auroc, probs = 0.25),
            q3 = quantile(auroc, probs = 0.75),
            label = formatC(median, format = 'f', digits = 3)) %>%
  ungroup()
rect = medians2 %>% filter(metric == 'Co-apex score') 
p2 = dat2 %>%
  mutate(metric = clean_metric(metric),
         metric = reorder(metric, auroc, stats::median),
         title = 'Measure of association') %>%
  ggplot(aes(y = metric, fill = '1', color = '1')) +
  facet_grid(~ title) +
  geom_blank() +
  geom_crossbarh(data = rect, aes(xmin = -Inf, xmax = +Inf, x = 0.5), 
                 # fill = 'grey90', fill = pal,
                 alpha = 0.3, color = NA, width = 1.05) +
  geom_vline(aes(xintercept = 0.5), color = 'black', size = 0.3, 
             linetype = 'dotted') +
  geom_boxploth(aes(x = auroc), width = 0.6, alpha = 0.6, outlier.shape = NA,
                coef = 0) + 
  geom_text(data = medians2, aes(x = 0.395, y = metric, label = label), 
            size = 2, color = 'grey20', hjust = 0) +
  scale_x_continuous('AUC', breaks = seq(0.45, 0.65, 0.05)) +
  coord_cartesian(xlim = c(0.4, 0.65)) +
  scale_fill_manual('', values = col, guide = F) +
  scale_color_manual('', values = col, guide = F) +
  boxed_theme() +
  theme(axis.title.y = element_blank())
p2
ggsave("fig/analysis/co_apex/GO-AUC-metrics-CORUM.pdf", p2,
       width = 6.3, height = 7.05, units = "cm", useDingbats = F)
