# Plot Supplementary Figure 7:
#' - n_datasets
#' - downsample_complexes
setwd("~/git/CF-MS-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

################################################################################
###### a. n_datasets
################################################################################

# read AUCs from the replicate_integration experiment
auc = readRDS("data/analysis/replicate_integration/AUC.rds") %>%
  # drop some (constant) columns
  dplyr::select(-min_fractions, -n_folds) %>%
  # pre-filter 1 dataset
  filter(n_datasets > 1)

df1 = filter(auc,
             feature_select == 'best_first',
             complex_set == 'held-out',
             split_by == 'proteins', 
             n_features <= 10)
pal = brewer.pal(12, 'Set3')
cols = pal[c(1, 10)]

pA = df1 %>%
  mutate(n_datasets = ifelse(n_datasets == 1, '1 dataset',
                             ifelse(n_datasets <= 10,
                                    paste(n_datasets, 'datasets'),
                                    n_datasets)),
         n_datasets = fct_relevel(n_datasets,
                                  '1 dataset',
                                  paste(2:10, 'datasets'),
                                  '15', '20', '30', '40'),
         classifier = fct_recode(classifier, 'Naive Bayes' = 'NB',
                                 'Random forest' = 'RF')) %>%
  ggplot(aes(x = factor(n_features), y = auc, fill = classifier, 
             color = classifier)) +
  facet_grid(~ n_datasets, space = 'free', scales = 'free') +
  geom_hline(aes(yintercept = 0.5), color = 'grey86') +
  geom_boxplot(width = 0.7, alpha = 0.6, outlier.shape = NA) + 
  scale_y_continuous('AUC') +
  scale_x_discrete('# of features') +
  scale_color_manual('', values = cols) +
  scale_fill_manual('', values = cols) +
  boxed_theme() +
  theme(legend.position = 'top',
        legend.key.size = unit(0.45, 'lines'))
pA
ggsave("fig/final/supp7/a_n_datasets.pdf", pA,
       width = 18.5, height = 5, units = "cm", useDingbats = FALSE)

################################################################################
###### b. downsample_complexes
################################################################################

# read the dataset
downs = readRDS("data/analysis/downsample_corum/AUC.rds") %>%
  # filter to the final grid
  filter(split_by == 'proteins', feature_select == 'best_first') %>%
  # filter to held-out complexes only
  filter(complex_set == 'held-out')

# convert outcomes to long
long = downs %>%
  mutate(corum_pairs = corum_TPs2 + corum_TNs2) %>%
  dplyr::select(split_by, classifier, n_features, n_datasets, downsample_pct,
                sample_idx, auc, corum_proteins2, corum_TPs2, corum_pairs) %>%
  gather('outcome', 'value', corum_proteins2:corum_pairs) %>%
  mutate(outcome = fct_recode(outcome,
                              '# of proteins' = 'corum_proteins2',
                              '# of true-positive PPIs' = 'corum_TPs2',
                              '# of total PPIs' = 'corum_pairs'))

# scale within sample_idx to % of maximum AUC
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

# plot
cols = brewer.pal(4, 'Blues')[-1]
pB = scaled %>%
  filter(outcome == '# of proteins') %>%
  mutate(classifier = fct_recode(classifier,
                                 'naive Bayes' = 'NB',
                                 'random forest' = 'RF'),
         n_datasets = fct_recode(as.character(n_datasets),
                                 '2 datasets' = '2',
                                 '3 datasets' = '3',
                                 '4 datasets' = '4')) %>%
  ggplot(aes(x = value, y = pct_of_max, color = factor(n_datasets))) +
  facet_grid(~ classifier) +
  geom_point(size = 0.2) + 
  scale_color_manual('', values = cols) +
  scale_x_continuous('# of complex proteins') +
  scale_y_continuous('% of maximum AUC', labels = function(x) x * 100) +
  guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                              override.aes = list(size = 0.8))) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.size = unit(0.5, 'lines'),
        legend.position = c(0.98, 0.03), 
        legend.justification = c(1, 0))
pB
ggsave("fig/final/supp7/b_downsample_corum.pdf", pB,
       width = 9, height = 4.75, units = "cm", useDingbats = FALSE)
