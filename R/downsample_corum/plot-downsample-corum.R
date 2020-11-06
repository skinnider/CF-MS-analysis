# Plot the effect of downsampling protein complexes on classifier performance
# to get a sense of the minimum number of complexes needed. 
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(minpack.lm)
source("R/theme.R")

# read the dataset
auc = readRDS("data/analysis/downsample_corum/AUC.rds") %>%
  # filter to the final grid
  filter(split_by == 'proteins', feature_select == 'best_first') %>%
  # filter to held-out complexes only
  filter(complex_set == 'held-out')

# print all variables
auc %>%
  extract(, -seq(which(colnames(.) == 'auc'), ncol(.))) %>%
  map(unique) %>%
  map(sort) %>%
  extract(lengths(.) > 1)

# convert to long
long = auc %>%
  mutate(corum_pairs = corum_TPs2 + corum_TNs2) %>%
  dplyr::select(split_by, classifier, n_features, n_datasets, downsample_pct,
                sample_idx, auc, corum_proteins2, corum_TPs2, corum_pairs) %>%
  gather('outcome', 'value', corum_proteins2:corum_pairs) %>%
  mutate(outcome = fct_recode(outcome,
                              '# of proteins' = 'corum_proteins2',
                              '# of true-positive PPIs' = 'corum_TPs2',
                              '# of total PPIs' = 'corum_pairs'))

# strong sample effect, so scale to % of max
scaled = long %>%
  group_by_at(vars(-auc, -value, -downsample_pct)) %>%
  mutate(max = max(auc), 
         q99 = quantile(auc, probs = 0.99),
         q95 = quantile(auc, probs = 0.95),
         pct_of_max = auc / max,
         pct_of_q99 = auc / q99,
         pct_of_q95 = auc / q95,
         n = n()) %>%
  ungroup()

# plot number of proteins
pal = pals::stepped3() %>% extract(c(10, 14))
## highlight low-data regime for total PPIs
cutoff = scaled %>%
  filter(pct_of_max < 0.98,
         classifier == 'NB',
         n_datasets == 4,
         outcome == '# of proteins') %>%
  group_by(classifier) %>%
  arrange(desc(value)) %>%
  dplyr::slice(1)
print(cutoff)
p1 = scaled %>%
  filter(n_datasets == 4, 
         outcome == '# of proteins') %>%
  mutate(classifier = fct_recode(classifier,
                                 'naive Bayes' = 'NB',
                                 'random forest' = 'RF')) %>%
  ggplot(aes(x = value, y = pct_of_max, color = classifier)) +
  # facet_grid(~ n_datasets) +
  geom_rect(data = cutoff, 
            aes(xmin = -Inf, xmax = value, ymin = -Inf, ymax = +Inf), 
            alpha = 0.3, fill = 'grey90', color = NA) +
  geom_vline(data = cutoff, aes(xintercept = value), linetype = 'dotted',
             color = 'grey20', size = 0.2) +
  geom_point(size = 0.2) + 
  scale_color_manual('', values = pal) +
  scale_x_continuous('# of proteins') +
  scale_y_continuous('AUC (% of maximum)', labels = function(x) x * 100) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.5, 'lines'),
        legend.position = c(0.97, 0.03), 
        legend.justification = c(1, 0))
p1
ggsave("fig/analysis/downsample_corum/n_proteins-n_datasets=4.pdf", p1,
       width = 4.75, height = 4.75, units = "cm", useDingbats = FALSE)

# also plot # of datasets
pal = brewer.pal(4, 'Blues')[-1]
p2 = scaled %>%
  filter(outcome == '# of proteins') %>%
  mutate(classifier = fct_recode(classifier,
                                 'naive Bayes' = 'NB',
                                 'random forest' = 'RF'),
         n_datasets = fct_recode(as.character(n_datasets),
                                 '2 datasets' = '2',
                                 '3 datasets' = '3',
                                 '4 datasets' = '4')) %>%
  ggplot(aes(x = value, y = pct_of_max, color = n_datasets)) +
  facet_grid(~ classifier, scales = 'free') +
  geom_point(size = 0.2) + 
  scale_color_manual('', values = pal) +
  scale_x_continuous('# of total PPIs (thousands)',
                     labels = function(x) x / 1e3) +
  scale_y_continuous('AUC (% of maximum)', labels = function(x) x * 100) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.5, 'lines'),
        legend.position = c(0.98, 0.03), 
        legend.justification = c(1, 0))
p2
ggsave("fig/analysis/downsample_corum/n_proteins-n_datasets.pdf", p2,
       width = 9, height = 4.75, units = "cm", useDingbats = FALSE)
